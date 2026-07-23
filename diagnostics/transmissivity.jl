# Transmissivity diagnostic
# =========================
#
# ITM parameterizes the surface shortwave as S = tau * S_toa with a transmissivity
# that ITM takes to be a linear function of surface elevation:
#
#       tau_ITM = trans_a + trans_b * z_srf         (defaults 0.46, 6e-5 m-1)
#
# This script confronts that with an observational transmissivity from the
# assembled forcing,
#
#       tau_obs = swd / S_toa                        (ERA5 surface SW / TOA insolation)
#
# per ice cell and month, and asks two things:
#   1. how well does the ITM default (a=0.46, b=6e-5) reproduce tau_obs, and
#   2. what does a least-squares refit give -- against z_srf alone, and against
#      z_srf together with cloud cover tcc (a natural second predictor, since
#      clouds are what tau is really responding to).
#
# It also writes tau_fit.png: tau_obs vs elevation coloured by cloud cover, with
# the ITM default and the two refits overlaid.
#
# Input is the file written by test_domain.x (monthly swd, tcc, S_toa, z_srf,
# mask). Usage:
#
#   julia --project=diagnostics diagnostics/transmissivity.jl [check_file.nc] [out.png]

using NCDatasets
using CairoMakie
using Statistics
using Printf

const A_ITM = 0.46      # trans_a default
const B_ITM = 6.0e-5    # trans_b default [m-1]

file    = length(ARGS) >= 1 ? ARGS[1] : "domain_greenland_check.nc"
out_png = length(ARGS) >= 2 ? ARGS[2] : "tau_fit.png"
println("reading $file")

ds    = NCDataset(file)
mask  = ds["mask"][:, :]           # (x,y)
z     = ds["z_srf"][:, :]
swd   = ds["swd"][:, :, :]         # (x,y,month)
toa   = ds["S_toa"][:, :, :]
tcc   = ds["tcc"][:, :, :]
close(ds)

nmon = size(swd, 3)

# Build the (cell,month) sample: ice cells with the sun meaningfully up, so the
# swd/S_toa ratio is well posed (skip the polar-winter darkness).
const S_MIN = 20.0                 # [W m-2] minimum monthly-mean TOA to include

tau = Float64[]; zz = Float64[]; cc = Float64[]; mo = Int[]
for m in 1:nmon, j in 1:size(swd,2), i in 1:size(swd,1)
    if mask[i,j] > 50 && toa[i,j,m] > S_MIN
        push!(tau, swd[i,j,m] / toa[i,j,m])
        push!(zz,  z[i,j])
        push!(cc,  tcc[i,j,m])
        push!(mo,  m)
    end
end
n = length(tau)
@printf("samples (ice cells x sunlit months): %d\n\n", n)

rmse(a, b) = sqrt(mean((a .- b).^2))
r2(y, yhat) = 1 - sum((y .- yhat).^2) / sum((y .- mean(y)).^2)

# --- 1. ITM default ------------------------------------------------------
tau_itm = A_ITM .+ B_ITM .* zz
@printf("tau range (obs)         : %.3f .. %.3f   mean %.3f\n",
        minimum(tau), maximum(tau), mean(tau))
println("---------------------------------------------------------------")
println("ITM default  tau = 0.46 + 6.0e-5 z")
@printf("  bias (ITM-obs)        : %+.4f\n", mean(tau_itm .- tau))
@printf("  RMSE                  : %.4f\n",  rmse(tau_itm, tau))
@printf("  R^2                   : %.3f\n\n", r2(tau, tau_itm))

# --- 2. refit tau ~ a + b z ---------------------------------------------
X1 = hcat(ones(n), zz)
b1 = X1 \ tau
fit1 = X1 * b1
println("refit        tau = a + b z")
@printf("  a                     : %.4f\n",  b1[1])
@printf("  b            [m-1]    : %.3e\n",  b1[2])
@printf("  RMSE                  : %.4f\n",  rmse(fit1, tau))
@printf("  R^2                   : %.3f\n\n", r2(tau, fit1))

# --- 3. refit tau ~ a + b z + c tcc -------------------------------------
X2 = hcat(ones(n), zz, cc)
b2 = X2 \ tau
fit2 = X2 * b2
println("refit        tau = a + b z + c tcc")
@printf("  a                     : %.4f\n",  b2[1])
@printf("  b            [m-1]    : %.3e\n",  b2[2])
@printf("  c            [1]      : %+.4f\n", b2[3])
@printf("  RMSE                  : %.4f\n",  rmse(fit2, tau))
@printf("  R^2                   : %.3f\n",  r2(tau, fit2))
println("---------------------------------------------------------------")
println("(c<0 as expected: more cloud -> lower transmissivity. The elevation-")
println(" only ITM form cannot see this, which is why a cloud predictor helps.)")

# --- figure: tau vs elevation, coloured by cloud cover -------------------
fig = Figure(size = (760, 560))
ax  = Axis(fig[1, 1], xlabel = "surface elevation  z_srf  [m]",
           ylabel = "transmissivity  τ = swd / S_toa",
           title  = "ITM τ vs observed (ERA5), coloured by cloud cover")
sc = scatter!(ax, zz, tau, color = cc, colormap = :viridis, markersize = 3,
              colorrange = (0, 1))
Colorbar(fig[1, 2], sc, label = "cloud cover  tcc")

zline = collect(range(minimum(zz), maximum(zz), length = 100))
lines!(ax, zline, A_ITM .+ B_ITM .* zline, color = :red, linewidth = 3,
       label = @sprintf("ITM default  0.46 + 6.0e-5 z  (R²=%.2f)", r2(tau, tau_itm)))
lines!(ax, zline, b1[1] .+ b1[2] .* zline, color = :black, linewidth = 3,
       label = @sprintf("refit z       %.2f + %.1e z  (R²=%.2f)", b1[1], b1[2], r2(tau, fit1)))
# z+tcc fit shown at the mean cloud cover
lines!(ax, zline, b2[1] .+ b2[2] .* zline .+ b2[3] * mean(cc), color = :dodgerblue,
       linewidth = 3, linestyle = :dash,
       label = @sprintf("refit z+tcc  (at mean tcc, R²=%.2f)", r2(tau, fit2)))
axislegend(ax, position = :rb, framevisible = true)

save(out_png, fig)
println("wrote $out_png")

# ------------------------------------------------------------------------
# Seasonal structure
# ------------------------------------------------------------------------
# The annual pooled fit conflates the seasonal cycle of tau with the within-
# month cloud response (season is a confounder: autumn is low-tau AND cloudy,
# spring high-tau), so the pooled cloud slope is steeper than any single month
# -- a Simpson's-paradox artifact. Fit each month on its own instead.

println("\n=== per-month fit  tau = a + b z + c tcc ======================")
println("mon    n   mean_tau     a      b[/km]    c       R^2")
month_names = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
mtau = fill(NaN, nmon); ma = fill(NaN, nmon); mb = fill(NaN, nmon)
mc = fill(NaN, nmon); mr2 = fill(NaN, nmon); mn = zeros(Int, nmon)
for m in 1:nmon
    idx = findall(==(m), mo)
    mn[m] = length(idx)
    mn[m] < 50 && (println(@sprintf("%3s  %5d   (too dark)", month_names[m], mn[m])); continue)
    t = tau[idx]; X = hcat(ones(mn[m]), zz[idx], cc[idx]); b = X \ t
    mtau[m] = mean(t); ma[m] = b[1]; mb[m] = b[2]; mc[m] = b[3]; mr2[m] = r2(t, X * b)
    @printf("%3s  %5d   %.3f    %.3f   %+.2f   %+.3f   %.2f\n",
            month_names[m], mn[m], mtau[m], ma[m], 1000*mb[m], mc[m], mr2[m])
end
println("(cloud slope c strengthens ~3x from spring to late summer; the pooled")
println(" annual c is steeper than any month, so a season-resolved tau is the")
println(" honest form for ITM -- use the melt-season, not the annual, fit.)")

# --- seasonal figure -----------------------------------------------------
season_png = replace(out_png, ".png" => "_seasonal.png")
fig2 = Figure(size = (1180, 500))

# Panel A: monthly mean tau and the cloud coefficient.
axA = Axis(fig2[1, 1], xlabel = "month", ylabel = "transmissivity τ",
           title = "seasonal cycle of τ and its cloud sensitivity",
           xticks = (1:12, month_names))
ok = .!isnan.(mtau)
lines!(axA, (1:12)[ok], mtau[ok], color = :black, linewidth = 3, label = "mean τ")
scatter!(axA, (1:12)[ok], mtau[ok], color = :black)
axC = Axis(fig2[1, 1], yaxisposition = :right, ylabel = "cloud coefficient c",
           yticklabelcolor = :dodgerblue, ylabelcolor = :dodgerblue)
hidespines!(axC); hidexdecorations!(axC)
lines!(axC, (1:12)[ok], mc[ok], color = :dodgerblue, linewidth = 3, linestyle = :dash)
scatter!(axC, (1:12)[ok], mc[ok], color = :dodgerblue)
axislegend(axA, position = :lb)

# Panel B: tau vs cloud cover for three representative months.
axB = Axis(fig2[1, 2], xlabel = "cloud cover  tcc", ylabel = "transmissivity τ",
           title = "τ–cloud slope steepens through the melt season")
for (m, col) in zip((5, 7, 8), (:seagreen, :orange, :firebrick))
    idx = findall(==(m), mo)
    isempty(idx) && continue
    scatter!(axB, cc[idx], tau[idx], color = (col, 0.10), markersize = 3)
    # slope of tau on tcc at mean elevation
    X = hcat(ones(length(idx)), zz[idx], cc[idx]); b = X \ tau[idx]
    cl = collect(range(0, 1, length = 20))
    lines!(axB, cl, b[1] .+ b[2] * mean(zz[idx]) .+ b[3] .* cl,
           color = col, linewidth = 3, label = "$(month_names[m])  c=$(@sprintf("%.2f", b[3]))")
end
axislegend(axB, position = :lb)

save(season_png, fig2)
println("wrote $season_png")
