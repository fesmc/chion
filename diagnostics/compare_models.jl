# Cross-config surface-SMB skill table (chion vs MAR), GRL-16KM.
#
#   julia --project=diagnostics diagnostics/compare_models.jl [out_dir]
#
# Reads every output/cmp_<name>/chion_grl16.nc produced by
# scripts/run_layer_comparison.sh, computes the MAR-comparable surface SMB
# (accum - d(runoff) - d(sublimation), final year; identical to compare_smb.jl)
# and reports bias/RMSE/R2 globally and per elevation zone. Also pulls wall time
# and per-column-step from each run.log. Writes a markdown table to stdout and
# a CSV (cmp_smb_stats.csv) to out_dir.

using NCDatasets, Statistics, Printf

chion   = "/Users/alrobi001/models/chion"
out_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(chion, "output")
mar_file = "/Users/alrobi001/models/ice_data/Greenland/GRL-16KM/GRL-16KM_MARv3.11-ERA_monmean_1961-1990.nc"

# name => label (order preserved)
configs = [
    "bessi_n15" => "BESSI n=15",
    "bessi_n7"  => "BESSI n=7",
    "bessi_n5"  => "BESSI n=5",
    "bessi_n3"  => "BESSI n=3",
    "bessi_n2"  => "BESSI n=2",
    "bessi_n1"  => "BESSI n=1",
    "pdd"       => "PDD",
    "itm"       => "ITM",
]

# --- MAR reference, forcing, geometry ------------------------------------
m     = NCDataset(mar_file)
mask  = Array(m["mask"][:, :])
z     = Array(m["z_srf"][:, :])
ice   = mask .> 50
mar   = dropdims(mean(m["smb"][:, :, :], dims=3), dims=3) .* 365.0            # mm/yr
accum = dropdims(mean(m["sf"][:, :, :] .+ m["rf"][:, :, :], dims=3), dims=3) .* 365.0
close(m)

# Elevation zones (docs/steady_state_snowpack.md).
zones = [
    ("all ice",           ice),
    ("margin z<800",      ice .& (z .< 800)),
    ("lower 800-1500",    ice .& (z .>= 800)  .& (z .< 1500)),
    ("mid 1500-2200",     ice .& (z .>= 1500) .& (z .< 2200)),
    ("interior z>2200",   ice .& (z .>= 2200)),
]

stats(x, y) = (chion = mean(x), mar = mean(y),
               bias = mean(x .- y),
               rmse = sqrt(mean((x .- y).^2)),
               r2   = 1 - sum((x .- y).^2) / sum((y .- mean(y)).^2),
               n    = length(x))

# Parse "wall time    [s] :  49.1540" and "per column-step  :  3.79E-07".
function run_timing(logpath)
    wall = perc = NaN
    ncol = steps = 0
    if isfile(logpath)
        for ln in eachline(logpath)
            occursin("wall time", ln)      && (wall  = parse(Float64, split(ln, ":")[end]))
            occursin("per column-step", ln) && (perc  = parse(Float64, split(ln, ":")[end]))
            occursin("columns          :", ln) && (ncol = parse(Int, split(ln, ":")[end]))
            occursin("steps            :", ln) && (steps = parse(Int, split(ln, ":")[end]))
        end
    end
    (; wall, perc, ncol, steps)
end

# --- collect ------------------------------------------------------------
rows = []          # (label, timing, Dict(zone => stats))
for (name, label) in configs
    nc = joinpath(chion, "output", "cmp_$name", "chion_grl16.nc")
    lg = joinpath(chion, "output", "cmp_$name", "run.log")
    if !isfile(nc)
        @warn "missing $nc — skipping"
        continue
    end
    d  = NCDataset(nc)
    ru = d["runoff"][:, :, :]
    # PDD/ITM have no energy-balance sublimation term; treat as zero (their
    # surface SMB is precip - runoff, matching how they define it).
    su = haskey(d, "sublimation") ? d["sublimation"][:, :, :] : zero(ru)
    close(d)
    smb = accum .- (ru[:, :, end] .- ru[:, :, end-1]) .- (su[:, :, end] .- su[:, :, end-1])

    zst = Dict{String,Any}()
    for (zn, zmask) in zones
        zst[zn] = stats(smb[zmask], mar[zmask])
    end
    push!(rows, (label, run_timing(lg), zst))
end

# --- speed table --------------------------------------------------------
println("\n## Speed (GRL-16KM, 50 yr, serial / 1 core)\n")
println("| config | wall [s] | per col-step [s] | vs BESSI n=15 |")
println("|---|---|---|---|")
base = findfirst(r -> r[1] == "BESSI n=15", rows)
basewall = base === nothing ? NaN : rows[base][2].wall
for (label, t, _) in rows
    speedup = isnan(basewall) ? "" : @sprintf("%.2fx", basewall / t.wall)
    @printf("| %s | %.1f | %.3e | %s |\n", label, t.wall, t.perc, speedup)
end

# --- skill table: global then per zone ----------------------------------
for (zn, _) in zones
    println("\n## SMB skill — $zn\n")
    println("| config | n | chion | MAR | bias | RMSE | R² |")
    println("|---|---|---|---|---|---|---|")
    for (label, _, zst) in rows
        s = zst[zn]
        @printf("| %s | %d | %+.0f | %+.0f | %+.0f | %.0f | %.2f |\n",
                label, s.n, s.chion, s.mar, s.bias, s.rmse, s.r2)
    end
end

# --- CSV ----------------------------------------------------------------
csv = joinpath(out_dir, "cmp_smb_stats.csv")
open(csv, "w") do io
    println(io, "config,wall_s,per_col_step_s,zone,n,bias,rmse,r2")
    for (label, t, zst) in rows
        for (zn, _) in zones
            s = zst[zn]
            @printf(io, "%s,%.3f,%.4e,%s,%d,%.2f,%.2f,%.4f\n",
                    label, t.wall, t.perc, zn, s.n, s.bias, s.rmse, s.r2)
        end
    end
end
println("\nwrote $csv")
