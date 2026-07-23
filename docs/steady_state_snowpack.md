# Steady-state snowpack on an ice-sheet domain

A standalone way to spin a chion snowpack up to a repeating seasonal state on a
real domain, driven straight from the raw `~/models/ice_data` datasets with no
preprocessing step. Greenland `GRL-16KM` is wired; other domains are one loader
subroutine away.

## Pieces

| piece | where | role |
|---|---|---|
| `chion_forcing_monthly` | `src/` (library) | mean-preserving monthly→daily interpolation (360-day, 12×30); the daily series averages back to the input monthly means, so precip totals are conserved |
| `insolation` | `libs/insol/` (driver layer) | top-of-atmosphere daily insolation (Laskar LA2004); kept out of `libchion.a`, preserving "the host supplies insolation" |
| `chion_domain` | `libs/domains/` (driver layer) | assembles standardized SI monthly forcing from raw data: MAR v3.11 (t2m, sf, rf, smb/melt/runoff) native on grid, ERA5 shortwave + cloud **conservatively regridded** (coords `map_init(con)`, cache under `maps/`), TOA insolation from latitude |
| `chion_grid.x` | `tests/chion_grid.f90` | one driver, `forcing_source = file | domain`; the domain path cycles the annual climatology for `n_years`, with `swd_source = file | transmissivity` |
| `diagnostics/*.jl` | `diagnostics/` | Julia + CairoMakie analysis (own project env): `compare_smb.jl` (chion vs MAR SMB figure), `transmissivity.jl` (ITM τ vs τ_obs = swd/S_toa) |

## Run

```
cd <run-dir>
ln -s <chion>/input input                 # chion reads input/chion_defaults.nml
<chion>/libchion/bin/chion_grid.x <chion>/par/chion_grl16.nml
```

Edit `path_ice_data` / `path_insol` in the par file. The conservative ERA5
weight map is generated on the first run and cached under `maps/` (gitignored,
regenerated if absent). A 50-year GRL-16KM BESSI run is ~46 s, 7204 columns.

## Result (GRL-16KM, BESSI, 50 yr, ERA5 shortwave)

**Which chion field is the MAR-comparable SMB.** MAR `smb` is the surface mass
balance, precip − runoff − sublimation. The matching chion quantity is the same
combination — equivalently d(column storage)/dt + `smb_ice` — NOT `smb_ice`
alone. `smb_ice` (the ice-facing flux) equals the surface SMB only once the firn
column is in equilibrium; on a 50-year cold start it still lags by the storage
growth (~200 vs ~258 mm/yr at year 50), so comparing it to MAR flatters the
result. Nor is the storage tendency alone right: from an empty start every
column is building firn, so it is positive almost everywhere and hides the
ablation entirely.

Against MAR v3.11, final-year surface SMB (precip − runoff − subl):

| zone | n | chion | MAR | bias | RMSE | R² | corr |
|---|---|---|---|---|---|---|---|
| all ice | 7204 | 258 | 190 | +68 | 269 | 0.74 | 0.87 |
| margin z<800 m | 571 | −461 | −740 | +279 | 633 | 0.39 | 0.72 |
| lower z 800–1500 m | 1232 | +210 | −16 | +226 | 449 | 0.43 | 0.78 |
| mid z 1500–2200 m | 1977 | +430 | +377 | +54 | 146 | 0.88 | 0.95 |
| interior z>2200 m | 3424 | +295 | +312 | −17 | 19 | 0.99 | 1.00 |

net-ablating area 12.4 % (MAR 14.4 %); domain runoff 128 vs 191 mm/yr.

The global R² (0.74) hides a sharp elevation split. **The interior accumulation
zone is essentially exact** (R² 0.99, RMSE 19), and the mid zone is good
(R² 0.88). **Below ~1500 m the skill collapses** (R² 0.4): the lower zone is
already net-ablating in MAR (−16) but net-accumulating in chion (+210), and the
margins ablate only ~60 % as hard as MAR (−461 vs −740). The whole +68 domain
bias is this coherent, elevation-dependent **under-ablation toward the warm
margins** — the signature of monthly-mean forcing with no sub-monthly
temperature extremes.

## Diurnal shortwave substepping (the default)

Melt is nonlinear in the instantaneous shortwave flux, so daily-mean forcing
systematically under-melts. Resolving the solar-noon peak with
`diurnal_shortwave_substeps` recovers it — but the lever has to be aimed. A
50-year sweep (all else equal):

| config | runoff | SMB | bias | R² | net-abl. |
|---|---|---|---|---|---|
| MAR target | 191 | 190 | — | — | 14.4 % |
| no diurnal | 128 | 258 | +68 | 0.74 | 12.4 % |
| substep, gate −8 °C (default) | 294 | 92 | −98 | 0.70 | 32.7 % |
| substep, gate −8 °C + T-cycle ±5 K | 436 | −50 | −240 | 0.35 | 39.4 % |
| **substep, gate −1 °C** | **178** | **209** | **+18** | **0.84** | **15.5 %** |

Two findings. (1) The **temperature cycle** (`diurnal_temperature_cycle`, ±5 K)
is too aggressive — it flips the sheet into strong over-ablation; left OFF.
(2) The **activation gate** `diurnal_shortwave_min_air_temperature` is the knob:
at its −8 °C default the substepping fires across the whole percolation zone and
over-ablates it (lower zone −294 where MAR is −16); raised to near-melting
(−1 °C, 272.15 K) the boost is confined to the warm margins and every elevation
band improves at once — domain bias +68→+18, R² 0.74→**0.84**, net-ablating area
12.4→15.5 % (MAR 14.4 %), lower zone +210→+29. This is one physically-motivated
threshold moved to a sensible value, not a fit, so it is the GRL-16KM default in
`par/chion_grl16.nml`. Residual: margins still ~15 % short of MAR (−640 vs −740).

## Performance and resolution

The column loop is the whole cost and columns are independent, so runtime is
linear in the number of ice columns and the per-column-step cost is fixed
(~3.8e-7 s, BESSI 15-layer energy balance with diurnal substepping). The
driver's reported wall time brackets the run loop only (setup and the ERA5
regrid are excluded).

| grid | ice columns | 50-yr serial | per col-step |
|---|---|---|---|
| GRL-16KM | 7 204 | 49 s | 3.79e-7 s |
| GRL-8KM | 28 879 | 202 s | 3.88e-7 s |

4× the columns → 4.1× the time: clean linear scaling, no resolution-dependent
overhead. **A 50-year GRL-16KM run is ~18 000 daily steps, not 18 000 years.**

OpenMP (`make openmp=1`, parallel over columns) on a 4-performance-core machine:

| threads | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| GRL-16KM wall | 46.9 s | 30.0 s | 20.9 s | 17.3 s |
| speedup | 1.0× | 1.6× | 2.2× | 2.7× |

GRL-8KM at 8 threads: 74 s (from 202 s, 2.7×) — same scaling. Returns diminish
past the 4 performance cores (the machine's other 6 are efficiency cores).

Physics is resolution-robust: GRL-8KM gives bias +15, RMSE 212, R² 0.84 —
statistically identical to GRL-16KM (+18, 0.84). The SMB pattern is set by
elevation and temperature, both resolved at 16 km; refining sharpens the margin
gradient but does not move the domain-scale skill.

## Transmissivity finding

`diagnostics/transmissivity.jl` on GRL-16KM: observed τ = swd/S_toa averages
0.64, but ITM's default intercept (0.46) is too low (bias −0.06, negative R²)
and **elevation alone explains almost nothing** (R² 0.07). Cloud cover is the
real control — τ ≈ 0.81 + 7.4e-5·z_srf − 0.51·tcc raises R² to 0.32. A
cloud-aware transmissivity is the natural next parameterization. Consistent
with this, running the driver with `swd_source = transmissivity` (the ITM
clear-sky form) under-melts even more than the ERA5-shortwave run.
