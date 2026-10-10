#!/usr/bin/env julia
"""
WP16 validation harness -- entry point.

    julia --project=validation validation/validate.jl [--quick]

Generates one forcing file, runs Chion.jl and both chion builds (sp and dp) on
it, and reports per-field differences. Exit status is non-zero if any gated
field exceeds its tolerance.

WHY BOTH CHION BUILDS
---------------------
The dp build is the one that is GATED. chion's production precision is sp, but
sp and Chion.jl's Float64 do not merely differ by a bounded per-field offset:
the two share discrete branch points (the mass_max split, the mass_min merge,
surface_has_snow, the melt gate) and not the round-off that decides which side
of them a given step falls on. Once any branch fires a step apart the layer
structure diverges by O(1) even though the physics is identical, so an sp gate
would fail on a correct port. Validating at dp removes that confound: a residual
is then attributable to the port.

The sp build is run too, and REPORTED rather than gated. That turns the sp-vs-dp
difference into a measured number -- PLAN.md section 3.1 measured candidate
expressions in isolation but never end-to-end -- and the `onset` column shows
where sp first parts company with the reference, which is the honest way to
present a trajectory divergence.

See README.md for the tolerance derivations.
"""

using Pkg
# CHION_VALIDATION_PROJECT selects another environment, e.g. one that
# develops a Chion.jl checkout on a different machine or at another commit
# (the reference is whatever Chion.jl that environment develops), so the
# committed Manifest.toml is left untouched.
Pkg.activate(get(ENV, "CHION_VALIDATION_PROJECT", @__DIR__))

include("forcing.jl")
include("compare.jl")
include("runners.jl")

const QUICK = "--quick" in ARGS
const WORKDIR = joinpath(@__DIR__, "work")

"""
Storage resolution of the REFERENCE files, and the floor on this comparison.

Chion.jl computes in Float64 but writes Float32 (`src/io.jl`: the NetCDF output
buffers are `Matrix{Float32}` / `Array{Float32,3}`). Every reference value is
therefore quantized to Float32 before we ever see it, so a difference of up to
half a Float32 ulp exists between the two models even if their arithmetic is
bit-identical. No comparison made through these files can resolve below this,
whatever precision either model runs at.

This is why the tolerances below are expressed in ulp of Float32 rather than of
the chion build's own `wp`: the file, not the model, sets the resolution.
"""
const EPS_FILE = eps(Float32)

"""
Gated tolerances for the dp build, as a multiple of EPS_FILE relative to each
field's own scale.

DERIVATION -- not fitted numbers.

At wp = dp, chion's own arithmetic contributes ~1e-16, four orders below the
Float32 storage floor, so it drops out entirely. What remains is:

  * the +-0.5 ulp quantization of each reference value on write, and
  * accumulation of that quantization through fields that are built from other
    fields (bulk_density is a mass sum over a thickness sum; Tsrf and albedo
    are carried forward step to step).

A single quantization is 0.5 ulp; a field formed from a handful of quantized
inputs stays within a small multiple. The gate is set at 4 ulp of Float32
(4.8e-07), which bounds that with margin while remaining far below anything a
real divergence would produce: a shifted layer split or a flipped melt branch
moves a field by O(1) relative, four to seven orders of magnitude above this.

The same gate applies to the layer-resolved fields. That is a deliberately
strong claim -- it asserts the layer STRUCTURE is identical, not merely similar
-- and it is testable at dp precisely because dp removes the chion-side
round-off that would otherwise shift a split or merge by a step.

A tolerance that has to be loosened is a signal to re-read the code, not to
loosen it further (PLAN.md WP16).
"""
const TOL_ULP_DP = 4.0

const GATED_VARS = vcat(BESSI_AGING_VARS, PDD_VARS, ITM_VARS)

function gate(diffs::Vector{FieldDiff}, label::AbstractString)
    nfail = 0
    tol = TOL_ULP_DP * EPS_FILE
    println()
    println("--- gate: $label  (tol = $(TOL_ULP_DP) ulp of Float32 = $(@sprintf("%.2E", tol))) ---")
    for d in diffs
        d.name in GATED_VARS || continue
        ulp_file = d.relscale / EPS_FILE
        if d.relscale <= tol
            @printf("  ok   : %-22s rel = %10.3E  = %5.2f ulp\n", d.name, d.relscale, ulp_file)
        else
            @printf("  FAIL : %-22s rel = %10.3E  = %5.2f ulp%s\n", d.name, d.relscale, ulp_file,
                    d.onset === nothing ? "" : "   (diverges at record $(d.onset))")
            nfail += 1
        end
    end
    return nfail
end

"""
Assert that each scenario reached the code paths it exists to exercise.

An agreement result is only as strong as the coverage behind it: two models can
agree perfectly on a column where almost nothing happened. These are the five
coverage requirements of PLAN.md WP16, stated as assertions rather than as
comments in the scenario table.
"""
function check_coverage(cov::Vector{NamedTuple})
    byname = Dict(c.name => c for c in cov)
    checks = [
        ("cold_dry builds layers",            () -> byname["cold_dry"].nmax > 1),
        ("cold_dry densifies",                () -> byname["cold_dry"].dmax > 350.0),
        ("cold_dry never melts",              () -> byname["cold_dry"].melt == 0.0),
        ("melting melts",                     () -> byname["melting"].melt > 0.0),
        ("melting refreezes",                 () -> byname["melting"].refr > 0.0),
        ("melting runs off",                  () -> byname["melting"].runoff > 0.0),
        ("bare_recover goes bare (N -> 0)",   () -> byname["bare_recover"].nmin == 0),
        ("bare_recover recovers (N > 0 later)", () -> byname["bare_recover"].nmax > 0),
        ("ntot_capacity reaches Ntot = 15",   () -> byname["ntot_capacity"].nmax == 15),
        ("ntot_capacity exports at the base", () -> byname["ntot_capacity"].mass_base > 0.0),
    ]
    nfail = 0
    println()
    println("--- coverage assertions ---")
    for (label, f) in checks
        ok = try f() catch; false end
        if ok
            println("  ok   : $label")
        else
            println("  FAIL : $label")
            nfail += 1
        end
    end
    return nfail
end

"""
ITM coverage, read off the chion dp+legacy output: each column must reach the
background-albedo branch and the part of the budget it exists for.
"""
function check_itm_coverage(chion_path::AbstractString, forcing_path::AbstractString,
                            names::Vector{String})
    local hs, alb, smbi, melt, refr, tsrf, tt
    NCDataset(forcing_path) do df
        tt = Array(df["TT"])                                  # (x, y, time)
    end
    NCDataset(chion_path) do dc
        tsrf, _ = read_canonical(dc, "Tsrf")
        hs, _ = read_canonical(dc, "H_snow")
        alb, _ = read_canonical(dc, "alb_s")
        smbi, _ = read_canonical(dc, "smb_ice")
        melt, _ = read_canonical(dc, "melt_cum")
        refr, _ = read_canonical(dc, "refreezing_cum")
    end
    # Record 1 is chion's initial state; columns are the scenarios in order.
    col(a, name) = Float64.(a[2:end, 1, findfirst(==(name), names)])
    steps(a, name) = diff(Float64.(a[:, 1, findfirst(==(name), names)]))
    p = ITM_PARAMS
    alb_land = p.alb_land * 0.5 + p.alb_forest * 0.5      # PDDA = 500
    checks = [
        ("ice_ablation melts out (H_snow -> 0)",
         () -> minimum(col(hs, "ice_ablation")) == 0.0),
        ("ice_ablation melts ice with no snow",
         () -> any((col(hs, "ice_ablation") .== 0.0) .& (steps(smbi, "ice_ablation") .< 0.0))),
        ("ice_ablation melts snow",
         () -> any((col(hs, "ice_ablation") .> 0.0) .& (steps(melt, "ice_ablation") .> 0.0))),
        ("ice_ablation refreezes",
         () -> last(col(refr, "ice_ablation")) > 0.0),
        ("ice_ablation firn warming fires (Tsrf > t2m, D27 path)",
         () -> (i = findfirst(==("ice_ablation"), names);
                any(col(tsrf, "ice_ablation") .> Float64.(tt[i, 1, :])))),
        ("ice_ablation snow-free albedo is alb_ice",
         () -> any((col(hs, "ice_ablation") .== 0.0) .& (col(alb, "ice_ablation") .== p.alb_ice))),
        ("ice_accum_cap holds H_snow at H_snow_max",
         () -> all(col(hs, "ice_accum_cap") .== p.H_snow_max)),
        ("ice_accum_cap exports the excess to ice",
         () -> last(col(smbi, "ice_accum_cap")) > 0.0),
        ("land_seasonal snow-free albedo is the land background",
         () -> any((col(hs, "land_seasonal") .== 0.0) .&
                   (abs.(col(alb, "land_seasonal") .- alb_land) .<= 1e-12))),
        ("ocean snow-free albedo is alb_ocean",
         () -> any((col(hs, "ocean") .== 0.0) .& (col(alb, "ocean") .== p.alb_ocean))),
    ]
    nfail = 0
    println()
    println("--- ITM coverage assertions ---")
    for (label, f) in checks
        ok = try f() catch; false end
        println(ok ? "  ok   : $label" : "  FAIL : $label")
        ok || (nfail += 1)
    end
    return nfail
end

"""
BESSI with `albedo = :aging` (Chion.jl 6d06af6), gated like the default
configuration: dp+legacy, same forcing, plus the `snow_age_days` output. The
coverage check asserts that the snow actually aged, under both timescales.
"""
function run_bessi_aging(fbessi::AbstractString)
    jl = run_julia_bessi(; forcing=fbessi, outfile="julia_bessi_aging.nc",
                         workdir=WORKDIR, ntot=15, years=1, albedo=:aging,
                         vars=BESSI_AGING_VARS)
    ch = run_chion(; precision=:dp, legacy=true, forcing=fbessi,
                   outfile="chion_bessi_aging_dp_legacy.nc", workdir=WORKDIR,
                   model="bessi", dt_out=1.0, dt=1.0, consts=AGING_CONSTS)
    d = compare_files(ch, jl, BESSI_AGING_VARS; eps_wp=eps_of(:dp))
    report(d, "BESSI port fidelity, albedo = aging: chion dp+legacy vs Chion.jl")
    nfail = gate(d, "BESSI port fidelity (albedo = aging)")

    age, alb, ts = NCDataset(ch) do dc
        (read_canonical(dc, "snow_age_days")[1], read_canonical(dc, "albedo")[1],
         read_canonical(dc, "Tsrf")[1])
    end
    aged = [a > 0 for a in skipmissing(age)]
    melting_aged = [a > 0 && t > 273.149 for (a, t) in zip(age, ts)
                    if !ismissing(a) && !ismissing(t)]
    checks = [("aging: snow ages (snow_age_days > 0)", any(aged)),
              ("aging: albedo relaxes below alpha_dry",
               any(a -> 0.70 < a < 0.81, skipmissing(alb))),
              ("aging: snow ages on a melting surface (tau_melt)", any(melting_aged))]
    println()
    println("--- coverage assertions (albedo = aging) ---")
    for (label, ok) in checks
        println(ok ? "  ok   : $label" : "  FAIL : $label")
        ok || (nfail += 1)
    end
    return nfail
end

"""
Relative humidity of the humidity-on BESSI configuration [1]: dry enough that
every column sublimates and evaporates, with the latent flux a few W m-2.
"""
const RH_HUMID = 0.7

"""
BESSI with humidity on (WP11): the default configuration's columns with uniform
RHZ and sea-level pressure, so the parameterized turbulent latent flux, its
phase-dependent latent heat (Lv at a melting surface, Lv+Lm below T0 and on
bare ice) and the gradient-based vapour mass (Chion.jl d0146e1) are exercised.
Gated like the default configuration: dp+legacy, every BESSI field.
"""
function run_bessi_humid(nstep::Int)
    fh = joinpath(WORKDIR, "forcing_bessi_rh.nc")
    write_forcing(fh, BESSI_SCENARIOS; nstep=nstep, dt_days=1.0,
                  relative_humidity=RH_HUMID)
    jl = run_julia_bessi(; forcing=fh, outfile="julia_bessi_rh.nc", workdir=WORKDIR,
                         ntot=15, years=1, humidity=true)
    ch = run_chion(; precision=:dp, legacy=true, forcing=fh,
                   outfile="chion_bessi_rh_dp_legacy.nc", workdir=WORKDIR,
                   model="bessi", dt_out=1.0, dt=1.0, rh_default=RH_HUMID)
    d = compare_files(ch, jl, BESSI_VARS; eps_wp=eps_of(:dp))
    report(d, "BESSI port fidelity, humidity on: chion dp+legacy vs Chion.jl")
    nfail = gate(d, "BESSI port fidelity (humidity on, rh = $(RH_HUMID))")

    subl, lhf = NCDataset(ch) do dc
        (read_canonical(dc, "sublimation")[1], read_canonical(dc, "latent_heat_flux_sum")[1])
    end
    names = [s.name for s in BESSI_SCENARIOS]
    final(a, name) = a[end, 1, findfirst(==(name), names)]
    checks = [("humid: the melting column sublimates",
               final(subl, "melting") > 0),
              ("humid: latent heat flux is non-zero",
               any(x -> x != 0, skipmissing(lhf)))]
    println()
    println("--- coverage assertions (humidity on) ---")
    for (label, ok) in checks
        println(ok ? "  ok   : $label" : "  FAIL : $label")
        ok || (nfail += 1)
    end
    return nfail
end

"""Thermal ice substrate layers of the substrate configuration (Chion.jl's default)."""
const ICE_SUBSTRATE_LAYERS = 5

"""
The substrate configuration's extra column. `bare_recover` goes bare only in
the melt season, after the melting snow has brought the ice beneath it to T0
(measured: every bare step melts, Tsrf = T0), so it never shows the cold
content of bare ice. This column never holds snow: bare ice through the year,
cooled below T0 in winter, melting in summer only once re-warmed.
"""
const BARE_ICE_SCENARIO = Scenario("bare_ice",
    "dry bare ice over the substrate: cold in winter, melts only once re-warmed",
    263.0, 15.0, 0.0, 150.0, 130.0)

"""
BESSI with the thermal ice substrate (Chion.jl 03bb445, plan C3): 5 layers of
ice below the snow/firn on both sides, gated like the default configuration
(dp+legacy, every BESSI field), on the default columns plus `bare_ice`. The
forcing carries a uniform HI = 1000 m, because chion puts the substrate only
under ice (porting_notes D34); Chion.jl does not read HI and has it under every
column. Coverage: `bare_recover` goes bare over the substrate; `bare_ice` sits
below T0 without melting while cold, and melts only later -- the cold content
of bare ice, which the substrate-free path (surface held at T0) cannot produce.
"""
function run_bessi_substrate(nstep::Int)
    scen = vcat(BESSI_SCENARIOS, [BARE_ICE_SCENARIO])
    fs = joinpath(WORKDIR, "forcing_bessi_ice.nc")
    write_forcing(fs, scen; nstep=nstep, dt_days=1.0, h_ice=1000.0)
    jl = run_julia_bessi(; forcing=fs, outfile="julia_bessi_ice.nc", workdir=WORKDIR,
                         ntot=15, years=1, ice_substrate_layers=ICE_SUBSTRATE_LAYERS)
    ch = run_chion(; precision=:dp, legacy=true, forcing=fs,
                   outfile="chion_bessi_ice_dp_legacy.nc", workdir=WORKDIR,
                   model="bessi", dt_out=1.0, dt=1.0, name_hice="HI",
                   bessi=(ice_substrate_layers=ICE_SUBSTRATE_LAYERS,))
    d = compare_files(ch, jl, BESSI_VARS; eps_wp=eps_of(:dp))
    report(d, "BESSI port fidelity, ice substrate: chion dp+legacy vs Chion.jl")
    nfail = gate(d, "BESSI port fidelity (ice substrate, $(ICE_SUBSTRATE_LAYERS) layers)")

    N, ts, melt = NCDataset(ch) do dc
        (read_canonical(dc, "N")[1], read_canonical(dc, "Tsrf")[1],
         read_canonical(dc, "melt")[1])
    end
    # Record 1 is the initial state; record k+1 follows step k.
    names = [s.name for s in scen]
    function bare_stats(name)
        i = findfirst(==(name), names)
        bare = Float64.(N[2:end, 1, i]) .== 0.0
        tsrf = Float64.(ts[2:end, 1, i])
        dm = diff(Float64.(melt[:, 1, i]))
        cold = findall(bare .& (tsrf .< 273.0) .& (dm .== 0.0))
        melting = findall(bare .& (dm .> 0.0))
        println("      $(rpad(name, 13)) $(count(bare)) bare steps, $(length(cold)) below " *
                "T0 without melt, $(length(melting)) melting; min bare Tsrf = " *
                "$(any(bare) ? minimum(tsrf[bare]) : NaN)")
        return bare, cold, melting
    end
    println()
    bare_r, _, _ = bare_stats("bare_recover")
    bare_i, cold_i, melt_i = bare_stats("bare_ice")
    checks = [("ice substrate: bare_recover goes bare", any(bare_r)),
              ("ice substrate: bare_ice stays bare", all(bare_i)),
              ("ice substrate: bare ice below T0 without melt (cold content)",
               !isempty(cold_i)),
              ("ice substrate: bare ice melts only after re-warming",
               !isempty(cold_i) && !isempty(melt_i) && first(melt_i) > first(cold_i))]
    println()
    println("--- coverage assertions (ice substrate) ---")
    for (label, ok) in checks
        println(ok ? "  ok   : $label" : "  FAIL : $label")
        ok || (nfail += 1)
    end
    return nfail
end

"""
Fine near-surface layer thicknesses [m] of the fine-layer configurations
(Chion.jl's default `near_surface_layer_max_thicknesses_m`).
"""
const NEAR_SURFACE_THICKNESSES = (0.02, 0.05, 0.10, 0.30)

"""chion's `&bessi` entry for `NEAR_SURFACE_THICKNESSES`."""
const FINE_BESSI = (near_surface_layer_max_thicknesses=NEAR_SURFACE_THICKNESSES,)

"""
BESSI with fine near-surface layers (Chion.jl 03bb445, plan C4): the top four
layers held at `NEAR_SURFACE_THICKNESSES` by the conservative remesh, the
100 kg m-2 surface merge off. Gated like the default configuration (dp+legacy,
every BESSI field, same forcing). Coverage, on chion's output: the limited
layers sit at their target thickness at the end of every step with a layer
below them; cap-down fires (snowfall pushed into layer 5 of `cold_dry`, no
melt); fill-up fires (layer 5 gives mass while `melting` melts); and layer 5,
which Chion.jl never splits or merges (N7; chion's C4b split/merge is reverted
under legacy), outgrows `mass_max` while the columns stay at <= 5 layers.
"""
function run_bessi_fine(fbessi::AbstractString)
    jl = run_julia_bessi(; forcing=fbessi, outfile="julia_bessi_fine.nc", workdir=WORKDIR,
                         ntot=15, years=1, near_surface=NEAR_SURFACE_THICKNESSES)
    ch = run_chion(; precision=:dp, legacy=true, forcing=fbessi,
                   outfile="chion_bessi_fine_dp_legacy.nc", workdir=WORKDIR,
                   model="bessi", dt_out=1.0, dt=1.0, bessi=FINE_BESSI)
    d = compare_files(ch, jl, BESSI_VARS; eps_wp=eps_of(:dp))
    report(d, "BESSI port fidelity, fine near-surface layers: chion dp+legacy vs Chion.jl")
    nfail = gate(d, "BESSI port fidelity (fine near-surface layers)")

    N, m, rho, melt = NCDataset(ch) do dc
        (read_canonical(dc, "N")[1], read_canonical(dc, "mass")[1],
         read_canonical(dc, "density")[1], read_canonical(dc, "melt")[1])
    end
    names = [s.name for s in BESSI_SCENARIOS]
    h = collect(NEAR_SURFACE_THICKNESSES)
    nk = length(h)
    # Record 1 is the initial state; record t follows step t-1.
    nrec = size(N, 1)
    nchecked, noff = 0, 0
    for i in eachindex(names), t in 2:nrec
        n = Int(N[t, 1, i])
        n >= 2 || continue
        nchecked += 1
        any(abs(m[t, k, 1, i] / rho[t, k, 1, i] - h[k]) > 1e-9 * h[k]
            for k in 1:min(nk, n - 1)) && (noff += 1)
    end
    col(name) = findfirst(==(name), names)
    i = col("cold_dry")
    cap = count(t -> N[t, 1, i] >= 5 && N[t-1, 1, i] >= 5 && m[t, 5, 1, i] > m[t-1, 5, 1, i],
                3:nrec)
    i = col("melting")
    fill = count(t -> N[t, 1, i] >= 5 && N[t-1, 1, i] >= 5 && melt[t, 1, i] > melt[t-1, 1, i] &&
                      m[t, 5, 1, i] < m[t-1, 5, 1, i], 3:nrec)
    nmax = Int(maximum(skipmissing(N)))
    m5max = maximum(m[t, 5, 1, i] for i in eachindex(names), t in 1:nrec if N[t, 1, i] >= 5)
    println()
    println("      limited layers checked on $(nchecked) column-steps, $(noff) off target; " *
            "cap-down into layer 5 on $(cap) cold_dry steps; fill-up from layer 5 on " *
            "$(fill) melting steps; max N = $(nmax), max layer-5 mass = $(round(m5max; digits=1))")
    checks = [("fine: limited layers at their target thickness after every step",
               nchecked > 0 && noff == 0),
              ("fine: cap-down fires (snowfall pushed into layer 5, cold_dry)", cap > 0),
              ("fine: fill-up fires (layer 5 refills the top while melting)", fill > 0),
              ("fine: layer 5 is never split (N <= 5, layer 5 beyond mass_max; N7)",
               nmax <= 5 && m5max > 500.0)]
    println()
    println("--- coverage assertions (fine near-surface layers) ---")
    for (label, ok) in checks
        println(ok ? "  ok   : $label" : "  FAIL : $label")
        ok || (nfail += 1)
    end
    return nfail
end

"""
Fine near-surface layers over the thermal ice substrate (both Chion.jl
defaults since 03bb445), on the substrate configuration's forcing and columns
(uniform HI = 1000 m for chion, D34; plus `bare_ice`). Gated only: the coverage
of each part is asserted in its own configuration.
"""
function run_bessi_fine_substrate(nstep::Int)
    scen = vcat(BESSI_SCENARIOS, [BARE_ICE_SCENARIO])
    fs = joinpath(WORKDIR, "forcing_bessi_ice.nc")
    write_forcing(fs, scen; nstep=nstep, dt_days=1.0, h_ice=1000.0)
    jl = run_julia_bessi(; forcing=fs, outfile="julia_bessi_fine_ice.nc", workdir=WORKDIR,
                         ntot=15, years=1, ice_substrate_layers=ICE_SUBSTRATE_LAYERS,
                         near_surface=NEAR_SURFACE_THICKNESSES)
    ch = run_chion(; precision=:dp, legacy=true, forcing=fs,
                   outfile="chion_bessi_fine_ice_dp_legacy.nc", workdir=WORKDIR,
                   model="bessi", dt_out=1.0, dt=1.0, name_hice="HI",
                   bessi=(ice_substrate_layers=ICE_SUBSTRATE_LAYERS, FINE_BESSI...))
    d = compare_files(ch, jl, BESSI_VARS; eps_wp=eps_of(:dp))
    report(d, "BESSI port fidelity, fine layers + ice substrate: chion dp+legacy vs Chion.jl")
    return gate(d, "BESSI port fidelity (fine near-surface layers + ice substrate)")
end

"""
BESSI with the cloud-proxy longwave (Chion.jl 03bb445, plan C5): the
downwelling longwave from an emissivity in air temperature and a shortwave
cloudiness proxy, `n = 1 - SWD/(TOA*tau_clear(z))`, with chion's internal
daily TOA (the host TOA of D33 is not set). Gated like the default
configuration (dp+legacy, every BESSI field, same forcing). The TOA needs the
same day of year and solar longitude on both sides (see runners.jl). Coverage:
at LAT = 70 N the year has days on the shortwave proxy and polar-night days on
the night cloudiness, the proxy is unclamped on some days, and chion's result
differs from the graybody run `ch_graybody`.
"""
function run_bessi_cloud_proxy(fbessi::AbstractString, ch_graybody::AbstractString)
    jl = run_julia_bessi(; forcing=fbessi, outfile="julia_bessi_lwcp.nc", workdir=WORKDIR,
                         ntot=15, years=1, overrides=(longwave_scheme=:cloud_proxy,))
    ch = run_chion(; precision=:dp, legacy=true, forcing=fbessi,
                   outfile="chion_bessi_lwcp_dp_legacy.nc", workdir=WORKDIR,
                   model="bessi", dt_out=1.0, dt=1.0,
                   consts=(longwave_scheme="cloud_proxy",))
    d = compare_files(ch, jl, BESSI_VARS; eps_wp=eps_of(:dp))
    report(d, "BESSI port fidelity, cloud-proxy longwave: chion dp+legacy vs Chion.jl")
    nfail = gate(d, "BESSI port fidelity (cloud-proxy longwave)")

    # Which cloudiness branch each day takes, from Chion.jl's own TOA and
    # calendar (record k <-> day of year k on the harness axis).
    lat, zs, sw = NCDataset(fbessi) do ds
        (Float64(ds["LAT"][1, 1]), Float64(ds["SH"][1, 1]), Array(ds["SWD"]))   # SWD (x, y, time)
    end
    nt = size(sw, 3)
    toa = [Chion._daily_toa_shortwave(lat, Chion._solar_longitude_deg_from_calendar_day(Float64(k)),
                                      Float64(k)) for k in 1:nt]
    daily = toa .> 50.0
    tau = 0.85 + 0.075 * zs / 1000
    inner = count(k -> daily[k] && any(0.0 < 1 - sw[i, 1, k] / (toa[k] * tau) < 1.0
                                       for i in axes(sw, 1)), 1:nt)
    println()
    println("      $(count(daily)) days on the shortwave proxy, $(nt - count(daily)) on the " *
            "night cloudiness; $(inner) days with an unclamped cloudiness")
    checks = [("cloud proxy: days on the shortwave cloudiness proxy", count(daily) > 0),
              ("cloud proxy: polar-night days on the night cloudiness", count(.!daily) > 0),
              ("cloud proxy: cloudiness unclamped on some days", inner > 0)]
    println()
    println("--- coverage assertions (cloud-proxy longwave) ---")
    for (label, ok) in checks
        println(ok ? "  ok   : $label" : "  FAIL : $label")
        ok || (nfail += 1)
    end
    return nfail + check_moved(ch, ch_graybody, "cloud-proxy longwave", "the graybody run")
end

"""
BESSI with `seb_scheme = :semix` and `turbulent_flux_scheme = :bessi` (Chion.jl
d0146e1 split, plan C6): the longwave absorbed with the surface emissivity
(`eps_snow`, `eps_ice` on bare ice), BESSI's turbulence. Gated like the default
configuration (dp+legacy, every BESSI field, same forcing); coverage: the
result differs from the `:bessi` longwave run `ch_bessi`.
"""
function run_bessi_seb_semix(fbessi::AbstractString, ch_bessi::AbstractString)
    jl = run_julia_bessi(; forcing=fbessi, outfile="julia_bessi_seb_semix.nc", workdir=WORKDIR,
                         ntot=15, years=1, overrides=(seb_scheme=:semix,))
    ch = run_chion(; precision=:dp, legacy=true, forcing=fbessi,
                   outfile="chion_bessi_seb_semix_dp_legacy.nc", workdir=WORKDIR,
                   model="bessi", dt_out=1.0, dt=1.0, consts=(seb_scheme="semix",))
    d = compare_files(ch, jl, BESSI_VARS; eps_wp=eps_of(:dp))
    report(d, "BESSI port fidelity, seb_scheme = semix, turbulence bessi: chion dp+legacy vs Chion.jl")
    nfail = gate(d, "BESSI port fidelity (seb_scheme = semix, turbulent_flux_scheme = bessi)")
    return nfail + check_moved(ch, ch_bessi, "seb_scheme = semix", "the bessi longwave run")
end

"""
BESSI with Chion.jl's calibrated surface scheme, `seb_scheme = :semix` and
`turbulent_flux_scheme = :semix` (Julia's bulk turbulence at its 03bb445
defaults: sensible exchange factor 2.5, stable coefficient 40; plan C7), with
humidity on so the latent exchange and its vapour mass run too. Two
configurations: the default columns, and the ice-substrate columns (5 layers,
plus `bare_ice`), where the bare-ice latent heat of D35 and the ice roughness
apply. Both gated like the default configuration (dp+legacy: Julia's 287.05 and
bare-ice Lv, D35/D38); coverage: each differs from chion's run with BESSI
turbulence on the same forcing, so the switch is seen to act.
"""
function run_bessi_turb_semix(nstep::Int; substrate::Bool=false)
    scen = substrate ? vcat(BESSI_SCENARIOS, [BARE_ICE_SCENARIO]) : BESSI_SCENARIOS
    tag = substrate ? "turb_semix_ice" : "turb_semix"
    ft = joinpath(WORKDIR, "forcing_bessi_$(tag).nc")
    write_forcing(ft, scen; nstep=nstep, dt_days=1.0, relative_humidity=RH_HUMID,
                  h_ice=substrate ? 1000.0 : nothing)
    jl = run_julia_bessi(; forcing=ft, outfile="julia_bessi_$(tag).nc", workdir=WORKDIR,
                         ntot=15, years=1, humidity=true,
                         ice_substrate_layers=substrate ? ICE_SUBSTRATE_LAYERS : 0,
                         overrides=(seb_scheme=:semix, turbulent_flux_scheme=:semix))
    chion_run(turb, out) = run_chion(; precision=:dp, legacy=true, forcing=ft,
        outfile=out, workdir=WORKDIR, model="bessi", dt_out=1.0, dt=1.0,
        rh_default=RH_HUMID, name_hice=substrate ? "HI" : "None",
        bessi=substrate ? (ice_substrate_layers=ICE_SUBSTRATE_LAYERS,) : (;),
        consts=(seb_scheme="semix", turbulent_flux_scheme=turb))
    ch = chion_run("semix", "chion_bessi_$(tag)_dp_legacy.nc")
    ch_ref = chion_run("bessi", "chion_bessi_$(tag)_ref_dp_legacy.nc")
    label = "seb_scheme = turbulent_flux_scheme = semix, humidity on" *
            (substrate ? ", ice substrate" : "")
    d = compare_files(ch, jl, BESSI_VARS; eps_wp=eps_of(:dp))
    report(d, "BESSI port fidelity, $label: chion dp+legacy vs Chion.jl")
    nfail = gate(d, "BESSI port fidelity ($label)")
    return nfail + check_moved(ch, ch_ref, "semix turbulence" * (substrate ? ", ice substrate" : ""),
                               "the BESSI-turbulence run")
end

"""
BESSI with diurnal substepping (plan C8), Chion.jl's calibrated 03bb445 set on
both sides (`DIURNAL_CALIBRATED`: up to 8 substeps, a 1 K temperature cycle,
none below -8 C), on the default forcing and columns with dynamic albedo, the
graybody longwave and BESSI's surface scheme. Day of year and solar longitude
agree by construction (runners.jl). A day that is not split keeps its forcing
in chion but takes the interval average in Chion.jl, which zeroes the
synthetic forcing's polar-night shortwave; `legacy_chion` does the same (D39).
Gated like the default configuration (dp+legacy, every BESSI field).
Coverage: substepped days (counted with Chion.jl's own criterion), unsplit days
with shortwave in the polar night, and a result that differs from the run
without substeps, `ch_daily`.
"""
function run_bessi_diurnal(fbessi::AbstractString, ch_daily::AbstractString)
    jl = run_julia_bessi(; forcing=fbessi, outfile="julia_bessi_diurnal.nc", workdir=WORKDIR,
                         ntot=15, years=1, diurnal=DIURNAL_CALIBRATED)
    ch = run_chion(; precision=:dp, legacy=true, forcing=fbessi,
                   outfile="chion_bessi_diurnal_dp_legacy.nc", workdir=WORKDIR,
                   model="bessi", dt_out=1.0, dt=1.0, diurnal=DIURNAL_CALIBRATED)
    d = compare_files(ch, jl, BESSI_VARS; eps_wp=eps_of(:dp))
    report(d, "BESSI port fidelity, diurnal substeps (8, 1 K): chion dp+legacy vs Chion.jl")
    nfail = gate(d, "BESSI port fidelity (diurnal substeps, 8 substeps, 1 K cycle)")

    # Days split, with Chion.jl's criterion on the file's forcing and calendar.
    lat, tt, sw = NCDataset(fbessi) do ds
        (Float64(ds["LAT"][1, 1]), Array(ds["TT"]), Array(ds["SWD"]))   # (x, y, time)
    end
    dc = DIURNAL_CALIBRATED
    nsub(i, k) = Chion._diurnal_shortwave_substep_count(
        1.0, sw[i, 1, k], tt[i, 1, k], dc.min_air_temperature_c + 273.15, lat,
        Chion._solar_longitude_deg_from_calendar_day(Float64(k)), dc.threshold, dc.max_substeps)
    nt = size(sw, 3)
    split = count(nsub(i, k) == dc.max_substeps for i in axes(sw, 1), k in 1:nt)
    night = count(sw[i, 1, k] > 0 &&
                  Chion._daily_toa_shortwave(lat, Chion._solar_longitude_deg_from_calendar_day(Float64(k)),
                                             Float64(k)) == 0
                  for i in axes(sw, 1), k in 1:nt)
    println()
    println("      $(split) column-days split into $(dc.max_substeps) substeps; $(night) " *
            "column-days with shortwave in the polar night (unsplit)")
    checks = [("diurnal: days split into substeps", split > 0),
              ("diurnal: unsplit polar-night days with shortwave (D39 path)", night > 0)]
    println()
    println("--- coverage assertions (diurnal substeps) ---")
    for (label, ok) in checks
        println(ok ? "  ok   : $label" : "  FAIL : $label")
        ok || (nfail += 1)
    end
    return nfail + check_moved(ch, ch_daily, "diurnal substeps", "the daily run")
end

"""
BESSI with each model's own defaults (plan C11): Chion.jl's `BESSIModel(grid)`
and chion's `input/chion_defaults.nml`, nothing pinned (`defaults = true`), so
the gate covers Chion.jl's calibrated 03bb445 set as one package: alpha_ice
0.40, the cloud-proxy longwave, `seb_scheme = turbulent_flux_scheme = semix`
(2.5, 40), the 5-layer ice substrate, the fine near-surface layers and 8
diurnal substeps with a 1 K cycle. Forcing: the semix-turbulence substrate
configuration's (default columns plus `bare_ice`, humidity on, HI = 1000 m),
only what the physics needs: chion puts the substrate only under ice (D34),
and the latent exchange needs a humidity. `legacy_chion` reverts chion's
deviations on this path (C4b's split below the fine layers, D35, D38, D39).
Gated like the default configuration (dp+legacy, every BESSI field). Coverage
on chion's output: layer 1 holds its fine thickness, bare ice cools below T0
over the substrate, and the result differs from the semix-turbulence substrate
run `ch_turb_ice` (same forcing, the other options pinned).
"""
function run_bessi_defaults(nstep::Int)
    scen = vcat(BESSI_SCENARIOS, [BARE_ICE_SCENARIO])
    fd = joinpath(WORKDIR, "forcing_bessi_defaults.nc")
    write_forcing(fd, scen; nstep=nstep, dt_days=1.0, relative_humidity=RH_HUMID,
                  h_ice=1000.0)
    jl = run_julia_bessi(; forcing=fd, outfile="julia_bessi_defaults.nc", workdir=WORKDIR,
                         years=1, humidity=true, defaults=true)
    ch = run_chion(; precision=:dp, legacy=true, forcing=fd,
                   outfile="chion_bessi_defaults_dp_legacy.nc", workdir=WORKDIR,
                   model="bessi", dt_out=1.0, dt=1.0, rh_default=RH_HUMID,
                   name_hice="HI", defaults=true)
    d = compare_files(ch, jl, BESSI_VARS; eps_wp=eps_of(:dp))
    report(d, "BESSI port fidelity, Chion.jl defaults: chion dp+legacy vs Chion.jl")
    nfail = gate(d, "BESSI port fidelity (Chion.jl defaults, calibrated 03bb445 set)")

    N, m, rho, ts = NCDataset(ch) do dc
        (read_canonical(dc, "N")[1], read_canonical(dc, "mass")[1],
         read_canonical(dc, "density")[1], read_canonical(dc, "Tsrf")[1])
    end
    names = [s.name for s in scen]
    h1 = first(NEAR_SURFACE_THICKNESSES)
    fine = [abs(m[t, 1, 1, i] / rho[t, 1, 1, i] - h1) <= 1e-9 * h1
            for i in eachindex(names), t in 2:size(N, 1) if N[t, 1, i] >= 2]
    i = findfirst(==("bare_ice"), names)
    cold_bare = any(N[t, 1, i] == 0 && ts[t, 1, i] < 273.0 for t in 2:size(N, 1))
    checks = [("defaults: layer 1 at its fine thickness on every layered step",
               !isempty(fine) && all(fine)),
              ("defaults: bare ice below T0 over the substrate", cold_bare)]
    println()
    println("--- coverage assertions (Chion.jl defaults) ---")
    for (label, ok) in checks
        println(ok ? "  ok   : $label" : "  FAIL : $label")
        ok || (nfail += 1)
    end
    ch_turb_ice = joinpath(WORKDIR, "chion_bessi_turb_semix_ice_dp_legacy.nc")
    return nfail + check_moved(ch, ch_turb_ice, "Chion.jl defaults",
                               "the semix-turbulence substrate run")
end

"""
Coverage helper: assert that a configuration's chion output `ch` differs from
the run `ch_ref` without its switch (max |Tsrf difference| > 0.1 K), so a gate
pass cannot come from a switch that silently did nothing on either side.
"""
function check_moved(ch::AbstractString, ch_ref::AbstractString, label::AbstractString,
                     reflabel::AbstractString)
    tsrf(path) = NCDataset(dc -> read_canonical(dc, "Tsrf")[1], path)
    moved = maximum(abs(Float64(x)) for x in skipmissing(tsrf(ch) .- tsrf(ch_ref)))
    ok = moved > 0.1
    println()
    println("--- coverage assertions ($label) ---")
    println((ok ? "  ok   : " : "  FAIL : ") * "$label: result differs from $reflabel " *
            "(max |dTsrf| = $(round(moved; digits=2)) K)")
    return ok ? 0 : 1
end

"""
ITM against Chion.jl's ITMModel (ported from chion, 29eb867): dp+legacy, so
chion's tsrf uses the daily melt_net as Chion.jl does (D27), and all eight
written fields are gated.
"""
function run_itm_julia(nstep::Int)
    println("\n[3/3] ITM: $(length(ITM_SCENARIOS)) columns x $nstep daily steps")
    for (sc, g) in ITM_SCENARIOS
        println("      - $(rpad(sc.name, 16)) $(sc.what)")
    end
    fitm = joinpath(WORKDIR, "forcing_itm.nc")
    # Rain below T0 (t_snow_max < T0): the only way ITM's melt_net turns
    # positive at t2m < T0, which is where firn warming (D27) shows in tsrf.
    write_forcing(fitm, first.(ITM_SCENARIOS); nstep=nstep, dt_days=1.0,
                  t_snow_max=271.15, geometry=last.(ITM_SCENARIOS))

    jl_itm = run_julia_itm(; forcing=fitm, outfile="julia_itm.nc", workdir=WORKDIR)
    ch_itm = run_chion(; precision=:dp, legacy=true, forcing=fitm,
                       outfile="chion_itm_dp_legacy.nc", workdir=WORKDIR,
                       model="itm", dt_out=1.0, dt=1.0, nml_extra=itm_nml(),
                       name_hice="HI", name_pdds="PDDA")
    d = compare_files(ch_itm, jl_itm, ITM_VARS; eps_wp=eps_of(:dp))
    report(d, "ITM port fidelity: chion dp+legacy vs Chion.jl")
    nfail = gate(d, "ITM port fidelity")
    nfail += check_itm_coverage(ch_itm, fitm, [sc.name for (sc, g) in ITM_SCENARIOS])
    return nfail
end

"""
ITM's production reference is smbpal, and that comparison already exists as
the WP12 acceptance test. It is run here rather than reimplemented, so there is
one statement of smbpal equivalence rather than two that can drift apart.
"""
function run_itm_check()
    println("\n      ITM: smbpal equivalence (tests/test_itm.x, both precisions)")
    nfail = 0
    for prec in (:dp, :sp)
        exe = joinpath(CHION_ROOT, bindir(prec), "test_itm.x")
        if !isfile(exe)
            println("  FAIL : $exe not built ($(make_cmd("itm", prec)))")
            nfail += 1
            continue
        end
        # ignorestatus: a failing test exits non-zero; report it rather than abort.
        out = read(ignorestatus(Cmd(`$exe`; dir=CHION_ROOT)), String)
        worst = ""
        for ln in split(out, '\n')
            occursin("equivalence,", ln) && occursin("rel =", ln) && (worst = strip(ln))
        end
        if occursin("ALL CHECKS PASSED", out)
            println("  ok   : ITM vs smbpal, precision=$prec  ($(worst))")
        else
            println("  FAIL : ITM vs smbpal, precision=$prec")
            for ln in split(out, '\n')
                occursin("FAIL", ln) && println("         $(strip(ln))")
            end
            nfail += 1
        end
    end
    return nfail
end

function main()
    mkpath(WORKDIR)
    nfail = 0

    println("="^72)
    println(" chion WP16 validation: chion vs Chion.jl")
    println(" reference: ", reference_id())
    isempty(BESSI_SCHEME_PINS) || println(" BESSI pins: ", BESSI_SCHEME_PINS)
    println("="^72)

    # =================================================================
    # BESSI -- Chion.jl is authoritative, tight tolerances.
    # =================================================================
    nstep = QUICK ? 90 : 365
    println("\n[1/3] BESSI: $(length(BESSI_SCENARIOS)) columns x $nstep daily steps")
    for s in BESSI_SCENARIOS
        println("      - $(rpad(s.name, 16)) $(s.what)")
    end

    fbessi = joinpath(WORKDIR, "forcing_bessi.nc")
    write_forcing(fbessi, BESSI_SCENARIOS; nstep=nstep, dt_days=1.0)

    jl_bessi = run_julia_bessi(; forcing=fbessi, outfile="julia_bessi.nc",
                               workdir=WORKDIR, ntot=15, years=1)

    # PORT FIDELITY (gated): dp + legacy, so chion runs Chion.jl's own
    # constants and any residual is attributable to the port itself.
    ch_legacy = run_chion(; precision=:dp, legacy=true, forcing=fbessi,
                          outfile="chion_bessi_dp_legacy.nc", workdir=WORKDIR,
                          model="bessi", dt_out=1.0, dt=1.0)
    d = compare_files(ch_legacy, jl_bessi, BESSI_VARS; eps_wp=eps_of(:dp))
    report(d, "BESSI port fidelity: chion dp+legacy vs Chion.jl")
    nfail += gate(d, "BESSI port fidelity")

    cov = coverage(ch_legacy, [s.name for s in BESSI_SCENARIOS])
    report_coverage(cov)
    nfail += check_coverage(cov)

    nfail += run_bessi_aging(fbessi)
    nfail += run_bessi_humid(nstep)
    nfail += run_bessi_substrate(nstep)
    nfail += run_bessi_fine(fbessi)
    nfail += run_bessi_fine_substrate(nstep)
    nfail += run_bessi_cloud_proxy(fbessi, ch_legacy)
    nfail += run_bessi_seb_semix(fbessi, ch_legacy)
    nfail += run_bessi_turb_semix(nstep)
    nfail += run_bessi_turb_semix(nstep; substrate=true)
    nfail += run_bessi_diurnal(fbessi, ch_legacy)
    nfail += run_bessi_defaults(nstep)

    # PRECISION COST (reported): sp vs dp, chion against itself, so the number
    # is the cost of wp = sp alone with no reference-model effects mixed in.
    ch_sp = run_chion(; precision=:sp, forcing=fbessi,
                      outfile="chion_bessi_sp.nc", workdir=WORKDIR,
                      model="bessi", dt_out=1.0, dt=1.0)
    ch_dp = run_chion(; precision=:dp, forcing=fbessi,
                      outfile="chion_bessi_dp.nc", workdir=WORKDIR,
                      model="bessi", dt_out=1.0, dt=1.0)
    report(compare_files(ch_sp, ch_dp, BESSI_VARS; eps_wp=eps_of(:sp),
                         drop_first=false),
           "BESSI precision cost: chion sp vs chion dp (reported)")

    # PHYSICS CORRECTION (reported): the measured effect of chion's
    # full-precision gas constant vs Chion.jl's 8.314, plus gravity (D22, D25).
    report(compare_files(ch_dp, ch_legacy, BESSI_VARS; eps_wp=eps_of(:dp),
                         drop_first=false),
           "BESSI legacy constants: R, g vs Chion.jl's 8.314, 9.81 (reported)")

    # =================================================================
    # PDD -- Chion.jl is authoritative since it adopted chion's budget
    # (D23; Chion.jl ce6a68d). Both integrals, on monthly steps: `pism` is
    # the Calov-Greve expectation integral, `simple` the clipped mean.
    # =================================================================
    nmonth = QUICK ? 12 : 60
    println("\n[2/3] PDD: $nmonth monthly steps (dt = 30 d), pdd_method = simple and pism")
    fpdd = joinpath(WORKDIR, "forcing_pdd.nc")
    write_forcing(fpdd, BESSI_SCENARIOS; nstep=nmonth, dt_days=30.0)

    for method in (:simple, :pism)
        jl_pdd = run_julia_pdd(; forcing=fpdd, outfile="julia_pdd_$(method).nc",
                               workdir=WORKDIR, years=1, pdd_method=method)

        # PORT FIDELITY (gated) at dp. No legacy build is needed: none of
        # chion's deliberate corrections touches PDD.
        for prec in (:dp, :sp)
            ch = run_chion(; precision=prec, forcing=fpdd,
                           outfile="chion_pdd_$(method)_$(prec).nc", workdir=WORKDIR,
                           model="pdd", dt_out=30.0, dt=30.0,
                           pdd_method=string(method))
            d = compare_files(ch, jl_pdd, PDD_VARS; eps_wp=eps_of(prec))
            if prec === :dp
                report(d, "PDD port fidelity ($method): chion dp vs Chion.jl")
                nfail += gate(d, "PDD port fidelity ($method)")
            else
                report(d, "PDD ($method): chion sp vs Chion.jl (reported)")
            end

            # The full three-reservoir closure, read off chion's output. It
            # does not depend on the reference, so it stays gated in both
            # precisions.
            res, inp = pdd_closure(ch, fpdd)
            tol = prec === :dp ? 1.0e-9 : 4.0 * eps(Float32)
            println()
            println("--- gate: PDD mass closure ($method), precision=$prec ---")
            @printf("  d(swe)+d(smb_ice)+d(runoff) - (snowfall+rainfall)\n")
            @printf("    worst column: residual = %.4E on input %.4E  -> rel %.3E\n",
                    res, inp, res / max(inp, 1.0))
            if res / max(inp, 1.0) <= tol
                @printf("  ok   : mass closes to %.1E relative\n", tol)
            else
                @printf("  FAIL : mass closure exceeds %.1E relative\n", tol)
                nfail += 1
            end
        end
    end

    nfail += run_itm_julia(QUICK ? 365 : 3 * 365)
    nfail += run_itm_check()

    println()
    println("="^72)
    if nfail == 0
        println(" WP16: ALL GATED FIELDS PASSED")
    else
        println(" WP16: $nfail GATED FIELD(S) FAILED")
    end
    println("="^72)
    return nfail
end

exit(main() == 0 ? 0 : 1)
