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
Pkg.activate(@__DIR__)

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

const GATED_VARS = vcat(BESSI_VARS, PDD_VARS)

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
ITM has no Chion.jl counterpart -- `build_model(:itm, ...)` deliberately errors
there -- so its reference is smbpal, and that comparison already exists as the
WP12 acceptance test. It is run here rather than reimplemented, so there is one
statement of ITM equivalence rather than two that can drift apart.
"""
function run_itm_check()
    println("\n[3/3] ITM: smbpal equivalence (tests/test_itm.x, both precisions)")
    nfail = 0
    for prec in (:dp, :sp)
        exe = joinpath(CHION_ROOT, bindir(prec), "test_itm.x")
        if !isfile(exe)
            println("  FAIL : $exe not built (make itm" *
                    (prec === :dp ? " precision=dp" : "") * ")")
            nfail += 1
            continue
        end
        out = read(Cmd(`$exe`; dir=CHION_ROOT), String)
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

    # PHYSICS CORRECTION (reported): the measured effect of using the gas
    # constant in the densification Arrhenius terms (Chion.jl issue #18).
    report(compare_files(ch_dp, ch_legacy, BESSI_VARS; eps_wp=eps_of(:dp),
                         drop_first=false),
           "BESSI densification correction: gas constant vs 8.13 (reported)")

    # =================================================================
    # PDD -- Chion.jl for STRUCTURE only. Monthly steps exercise the
    # PISM/Calov-Greve expectation integral.
    # =================================================================
    nmonth = QUICK ? 12 : 60
    println("\n[2/3] PDD: $nmonth monthly steps (dt = 30 d, PISM branch)")
    fpdd = joinpath(WORKDIR, "forcing_pdd.nc")
    write_forcing(fpdd, BESSI_SCENARIOS; nstep=nmonth, dt_days=30.0)

    jl_pdd = run_julia_pdd(; forcing=fpdd, outfile="julia_pdd.nc",
                           workdir=WORKDIR, years=1)

    # PDD is NOT gated against Chion.jl. chion implements a different snowpack
    # budget on purpose (D23 / Chion.jl issue #19), so agreement with Chion.jl
    # is no longer the property worth asserting. It is gated on its own
    # mass-closure identity instead -- which is strictly stronger, and which
    # Chion.jl's PDD cannot satisfy at all.
    for prec in (:dp, :sp)
        ch = run_chion(; precision=prec, forcing=fpdd,
                       outfile="chion_pdd_$(prec).nc", workdir=WORKDIR,
                       model="pdd", dt_out=30.0, dt=30.0)
        report(compare_files(ch, jl_pdd, PDD_VARS; eps_wp=eps_of(prec)),
               "PDD, chion $prec vs Chion.jl (REPORTED, not gated -- " *
               "different budget by design)")

        res, inp = pdd_closure(ch, fpdd)
        tol = prec === :dp ? 1.0e-9 : 4.0 * eps(Float32)
        println()
        println("--- gate: PDD mass closure, precision=$prec ---")
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
