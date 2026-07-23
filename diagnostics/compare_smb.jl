# chion vs MAR surface-SMB comparison figure (CairoMakie).
#
#   julia --project=diagnostics diagnostics/compare_smb.jl <chion_output.nc> [mar_file.nc] [out.png]
#
# Surface SMB (MAR definition) = precip - runoff - sublimation, taken as the
# final-year value of a spun-up chion run. Writes a 4-panel PNG: chion map,
# MAR map, difference, and a scatter with bias / RMSE / R^2.

using NCDatasets
using CairoMakie
using Statistics
using Printf

chion_file = length(ARGS) >= 1 ? ARGS[1] : "chion_grl16.nc"
mar_file   = length(ARGS) >= 2 ? ARGS[2] :
    "/Users/alrobi001/models/ice_data/Greenland/GRL-16KM/GRL-16KM_MARv3.11-ERA_monmean_1961-1990.nc"
out_png    = length(ARGS) >= 3 ? ARGS[3] : "smb_chion_vs_mar.png"

# --- MAR reference and accumulation forcing ------------------------------
m     = NCDataset(mar_file)
mask  = Array(m["mask"][:, :])
xc    = Array(m["xc"][:]);  yc = Array(m["yc"][:])
ice   = mask .> 50
mar   = dropdims(mean(m["smb"][:, :, :], dims=3), dims=3) .* 365.0          # mm/yr
accum = dropdims(mean(m["sf"][:, :, :] .+ m["rf"][:, :, :], dims=3), dims=3) .* 365.0
close(m)

# --- chion surface SMB: accum - d(runoff) - d(subl), final year ----------
d  = NCDataset(chion_file)
ru = d["runoff"][:, :, :]
su = haskey(d, "sublimation") ? d["sublimation"][:, :, :] : zero(ru)   # PDD/ITM: no subl term
close(d)
chion = accum .- (ru[:, :, end] .- ru[:, :, end-1]) .- (su[:, :, end] .- su[:, :, end-1])

# --- statistics over ice -------------------------------------------------
x = chion[ice];  y = mar[ice]
bias = mean(x .- y)
rmse = sqrt(mean((x .- y).^2))
r2   = 1 - sum((x .- y).^2) / sum((y .- mean(y)).^2)
@printf("bias %+.1f  RMSE %.1f  R2 %.2f\n", bias, rmse, r2)

# masked fields for display (NaN off-ice -> transparent)
mask_off(v) = [ice[i, j] ? v[i, j] : NaN for i in axes(v, 1), j in axes(v, 2)]
chion_m = mask_off(chion);  mar_m = mask_off(mar);  diff_m = chion_m .- mar_m

# --- figure --------------------------------------------------------------
fig = Figure(size = (1500, 640))
smbrange = (-2000, 2000)
diffrange = (-800, 800)

function mappanel(pos, field, title, crange, cmap)
    ax = Axis(fig[1, pos], title = title, aspect = DataAspect(),
              xticksvisible = false, yticksvisible = false,
              xticklabelsvisible = false, yticklabelsvisible = false)
    hm = heatmap!(ax, xc, yc, field, colormap = cmap, colorrange = crange,
                  nan_color = :transparent)
    return hm
end

hm1 = mappanel(1, chion_m, "chion surface SMB", smbrange, :RdBu)
mappanel(2, mar_m, "MAR v3.11 SMB", smbrange, :RdBu)
hm3 = mappanel(3, diff_m, "chion − MAR", diffrange, :PuOr)
Colorbar(fig[2, 1:2], hm1, vertical = false, label = "mm w.e. yr⁻¹", flipaxis = false)
Colorbar(fig[2, 3], hm3, vertical = false, label = "mm w.e. yr⁻¹", flipaxis = false)

axs = Axis(fig[1, 4], title = @sprintf("bias %+.0f   RMSE %.0f   R² %.2f", bias, rmse, r2),
           xlabel = "MAR SMB", ylabel = "chion SMB", aspect = DataAspect())
scatter!(axs, y, x, markersize = 2, color = (:black, 0.12))
lines!(axs, [-4000, 3200], [-4000, 3200], color = :red)
limits!(axs, -4000, 3200, -4000, 3200)

Label(fig[0, :], basename(chion_file) * "  (ice cells, mm w.e. yr⁻¹)", fontsize = 16)

save(out_png, fig)
println("wrote $out_png")
