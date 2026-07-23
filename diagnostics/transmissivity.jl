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
# Input is the file written by test_domain.x (monthly swd, tcc, S_toa, z_srf,
# mask). Usage:
#
#   julia --project=validation diagnostics/transmissivity.jl [check_file.nc]

using NCDatasets
using Statistics
using Printf

const A_ITM = 0.46      # trans_a default
const B_ITM = 6.0e-5    # trans_b default [m-1]

file = length(ARGS) >= 1 ? ARGS[1] : "domain_greenland_check.nc"
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

tau = Float64[]; zz = Float64[]; cc = Float64[]
for m in 1:nmon, j in 1:size(swd,2), i in 1:size(swd,1)
    if mask[i,j] > 50 && toa[i,j,m] > S_MIN
        push!(tau, swd[i,j,m] / toa[i,j,m])
        push!(zz,  z[i,j])
        push!(cc,  tcc[i,j,m])
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
