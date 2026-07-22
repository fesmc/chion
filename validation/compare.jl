"""
Field-by-field comparison of a chion output file against a Chion.jl one.

THE TWO FILE-FORMAT HAZARDS (docs/porting_notes.md D14, D15)
------------------------------------------------------------
D14 -- axis order differs of necessity. Chion.jl declares its variables
("t","x","y") and ("t","layer","x","y"); chion writes (time,yc,xc) and
(time,layer,yc,xc), because netCDF requires the unlimited dimension to be
slowest-varying. Read back through NCDatasets the two therefore arrive in
DIFFERENT Julia orders -- (t,x,y) from Chion.jl, (x,y,time) from chion.

This module never hard-codes either. `read_canonical` permutes BY DIMENSION
NAME into one canonical order, so a future change to either writer's dimension
order cannot silently transpose a comparison into nonsense.

D15 -- chion fills unmapped cells with MV = -9999, Chion.jl with NaN. Both are
mapped to `missing`, and the comparison asserts the two files agree about WHICH
cells are missing before comparing any values. A field that is missing in one
and present in the other is a failure, not something to be quietly skipped.

RECORD ALIGNMENT
----------------
chion writes an initial record before any step, then one record per step;
Chion.jl writes only post-step records. So chion record k+1 and Chion.jl record
k are the same model state, and both carry the same time stamp -- which is
asserted rather than assumed.
"""

using NCDatasets
using Printf
using Statistics

const MV = -9999.0

"""The 18 variables that Chion.jl and chion both write for BESSI."""
const BESSI_VARS = ["thickness", "wet_mass", "bulk_density", "liquid_water",
                    "mass_base", "smb_ice", "runoff", "melt", "refreezing",
                    "sublimation", "latent_heat_flux_sum", "Tsrf", "albedo",
                    "N", "mass", "mass_w", "density", "temperature"]

"""The 4 variables both write for PDD."""
const PDD_VARS = ["snowpack_swe", "smb_ice", "runoff", "pdd_sum"]

const TIME_NAMES = ("t", "time")
const X_NAMES = ("x", "xc")
const Y_NAMES = ("y", "yc")
const LAYER_NAMES = ("layer",)

_find(dn, names) = findfirst(d -> d in names, dn)

"""
    read_canonical(ds, name) -> (data, has_layer)

Read a variable and permute it to canonical order: (time, y, x), or
(time, layer, y, x) when it is layer-resolved. Missing data (NaN or MV) becomes
`missing`.
"""
function read_canonical(ds, name::AbstractString)
    haskey(ds, name) || error("variable '$name' not present in $(path(ds))")
    v = ds[name]
    dn = collect(NCDatasets.dimnames(v))

    it = _find(dn, TIME_NAMES)
    ix = _find(dn, X_NAMES)
    iy = _find(dn, Y_NAMES)
    il = _find(dn, LAYER_NAMES)
    (it === nothing || ix === nothing || iy === nothing) &&
        error("variable '$name' has unrecognised dimensions $dn")

    a = Array(v)
    perm = il === nothing ? [it, iy, ix] : [it, il, iy, ix]
    a = permutedims(a, perm)

    out = Array{Union{Missing,Float64}}(undef, size(a))
    @inbounds for i in eachindex(a)
        x = Float64(a[i])
        out[i] = (isnan(x) || x == MV) ? missing : x
    end
    return out, il !== nothing
end

"""Per-field comparison result. `scale` is the reference field's own magnitude."""
struct FieldDiff
    name::String
    scale::Float64
    maxabs::Float64
    relscale::Float64
    ulps::Float64
    onset::Union{Nothing,Int}   # first record where relscale exceeds 1e-9
    npoints::Int
end

"""
    compare_files(chion_path, julia_path, vars; eps_wp)

Compare every variable in `vars`. `eps_wp` is the machine epsilon of the chion
build being tested, so `ulps` expresses each difference in units of that build's
own resolution rather than an absolute number that means different things in the
sp and dp builds.
"""
function compare_files(chion_path::AbstractString, julia_path::AbstractString,
                       vars::Vector{String}; eps_wp::Float64, drop_first::Bool=true)
    diffs = FieldDiff[]
    NCDataset(chion_path) do dc
        NCDataset(julia_path) do dj
            # --- record alignment, asserted ---------------------------------
            # Chion.jl writes NO time coordinate variable -- only a bare `t`
            # dimension (upstream defect: its output carries no record timing,
            # so a file cannot be interpreted without knowing the forcing that
            # produced it). Alignment is therefore checked structurally, on
            # counts, and chion's own axis is checked for uniformity.
            tc = Array(dc["time"])
            njl = dj.dim[haskey(dj.dim, "t") ? "t" : "time"]
            # drop_first=false compares two chion files, which both carry the
            # initial pre-step record, so the counts match exactly.
            length(tc) == njl + (drop_first ? 1 : 0) || error(
                "record-count mismatch: chion has $(length(tc)), Chion.jl has " *
                "$njl; expected chion = Chion.jl + 1 (chion's initial record). " *
                "Check dt_out == dt in the chion namelist.")
            if length(tc) > 2
                steps = diff(tc)
                maximum(steps) - minimum(steps) <= 1e-6 * maximum(abs.(steps)) || error(
                    "chion's time axis is not uniform: steps range " *
                    "$(minimum(steps))..$(maximum(steps)). dt_out must equal dt " *
                    "for a per-step comparison.")
            end

            for name in vars
                ac, layc = read_canonical(dc, name)
                aj, layj = read_canonical(dj, name)
                layc == layj || error("'$name' is layer-resolved in one file only")

                if drop_first
                    ac = layc ? ac[2:end, :, :, :] : ac[2:end, :, :]
                end

                size(ac) == size(aj) || error(
                    "'$name' shape mismatch after alignment: chion $(size(ac)) " *
                    "vs Chion.jl $(size(aj))")

                # D15: the two files must agree about which cells are missing.
                mc = ismissing.(ac); mj = ismissing.(aj)
                mc == mj || error(
                    "'$name': the files disagree about which cells are missing " *
                    "($(count(mc)) in chion vs $(count(mj)) in Chion.jl). MV and " *
                    "NaN are treated as equivalent, so this is a real mask difference.")

                idx = findall(.!mc)
                isempty(idx) && continue
                vc = Float64[ac[i] for i in idx]
                vj = Float64[aj[i] for i in idx]

                d = abs.(vc .- vj)
                scale = maximum(abs.(vj))
                scale = scale == 0 ? 1.0 : scale
                maxabs = maximum(d)
                rel = maxabs / scale

                # First record at which the fields visibly part company. This is
                # the diagnostic that matters when a difference is a trajectory
                # divergence rather than a uniform offset.
                onset = nothing
                nrec = size(ac, 1)
                for k in 1:nrec
                    sl = layc ? (ac[k, :, :, :], aj[k, :, :, :]) : (ac[k, :, :], aj[k, :, :])
                    kk = findall(.!ismissing.(sl[1]))
                    isempty(kk) && continue
                    dk = maximum(abs(Float64(sl[1][i]) - Float64(sl[2][i])) for i in kk)
                    if dk / scale > 1e-9
                        onset = k
                        break
                    end
                end

                push!(diffs, FieldDiff(name, scale, maxabs, rel,
                                       rel / eps_wp, onset, length(idx)))
            end
        end
    end
    return diffs
end

"""
    coverage(chion_path, names) -> Vector{NamedTuple}

Summarize, per column, which code paths a run actually reached.

This exists because a scenario that does not fire is a weak test that looks
like a strong one. The shipped column example had exactly that problem before
the D18 retuning: peak accumulation stayed below `mass_max`, so the layer
machinery never ran (peak layer count 1) and nothing refroze, while the run
still "passed". The harness therefore asserts coverage rather than assuming it.
"""
function coverage(chion_path::AbstractString, names::Vector{String})
    out = NamedTuple[]
    NCDataset(chion_path) do dc
        N, _ = read_canonical(dc, "N")
        melt, _ = read_canonical(dc, "melt")
        refr, _ = read_canonical(dc, "refreezing")
        runoff, _ = read_canonical(dc, "runoff")
        mbase, _ = read_canonical(dc, "mass_base")
        dens, _ = read_canonical(dc, "density")

        ncol = size(N, 3)
        for i in 1:ncol
            # Record 1 is the initial state, written before any step, when
            # every column is empty. Including it would make "this column went
            # bare" trivially true for every scenario.
            col(a) = collect(skipmissing(a[2:end, 1, i]))
            n = col(N)
            d = collect(skipmissing(vec(dens[2:end, :, 1, i])))
            dpos = filter(>(0), d)
            push!(out, (name = i <= length(names) ? names[i] : "col$i",
                        nmax = isempty(n) ? 0 : Int(maximum(n)),
                        nmin = isempty(n) ? 0 : Int(minimum(n)),
                        melt = last(col(melt)),
                        refr = last(col(refr)),
                        runoff = last(col(runoff)),
                        mass_base = last(col(mbase)),
                        dmin = isempty(dpos) ? 0.0 : minimum(dpos),
                        dmax = isempty(dpos) ? 0.0 : maximum(dpos)))
        end
    end
    return out
end

function report_coverage(cov::Vector{NamedTuple})
    println()
    println("--- scenario coverage (from the chion dp run) ---")
    @printf("%-16s %5s %5s %12s %12s %12s %12s %14s\n",
            "column", "Nmax", "Nmin", "melt", "refreezing", "runoff", "mass_base",
            "density range")
    for c in cov
        @printf("%-16s %5d %5d %12.1f %12.1f %12.1f %12.1f %6.0f -> %-6.0f\n",
                c.name, c.nmax, c.nmin, c.melt, c.refr, c.runoff, c.mass_base,
                c.dmin, c.dmax)
    end
end

function report(diffs::Vector{FieldDiff}, label::AbstractString)
    println()
    println("--- $label ---")
    @printf("%-22s %14s %14s %14s %10s %8s\n",
            "field", "scale", "max abs diff", "rel(scale)", "ulp(wp)", "onset")
    for d in diffs
        @printf("%-22s %14.4E %14.4E %14.4E %10.1f %8s\n",
                d.name, d.scale, d.maxabs, d.relscale, d.ulps,
                d.onset === nothing ? "-" : string(d.onset))
    end
end

"""
    pdd_closure(chion_path, forcing_path) -> (residual, input_total, terms)

End-to-end mass closure for chion's PDD, read straight off the output file.

    snowfall + rainfall == d(snowpack_swe) + d(smb_ice) + d(runoff)

This REPLACES the Chion.jl comparison as PDD's gate. chion's PDD deliberately
implements a different budget from Chion.jl's (docs/porting_notes.md D23,
Chion.jl issue #19), so agreement with Chion.jl is no longer the property worth
asserting -- and Chion.jl's own PDD cannot satisfy this identity, because it
credits `smb_ice` with `d(snowpack_swe)` as well and therefore counts the
reservoir twice.

Ice melt cancels between `smb_ice` (negative) and `runoff` (positive), which is
correct: it is mass drawn from the ice body below, not from anything the column
received.
"""
function pdd_closure(chion_path::AbstractString, forcing_path::AbstractString)
    local swe, smb, runoff, sf, rf, times
    NCDataset(chion_path) do dc
        swe, _ = read_canonical(dc, "snowpack_swe")
        smb, _ = read_canonical(dc, "smb_ice")
        runoff, _ = read_canonical(dc, "runoff")
        times = Array(dc["time"])
    end
    NCDataset(forcing_path) do df
        sf = Array(df["SF"]);  rf = Array(df["RF"])   # (x, y, time)
    end

    dt_sec = (times[2] - times[1]) * 86400.0
    ncol = size(swe, 3)
    worst = 0.0; worst_in = 0.0
    for i in 1:ncol
        input = sum(Float64(sf[i, 1, k]) + Float64(rf[i, 1, k])
                    for k in 1:size(sf, 3)) * dt_sec
        dswe = Float64(swe[end, 1, i]) - Float64(swe[1, 1, i])
        dsmb = Float64(smb[end, 1, i]) - Float64(smb[1, 1, i])
        drun = Float64(runoff[end, 1, i]) - Float64(runoff[1, 1, i])
        res = abs(dswe + dsmb + drun - input)
        if res / max(input, 1.0) > worst / max(worst_in, 1.0)
            worst = res; worst_in = input
        end
    end
    return worst, worst_in
end
