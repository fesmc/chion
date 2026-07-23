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
| `diagnostics/transmissivity.jl` | `diagnostics/` | confronts ITM's τ = a + b·z_srf with τ_obs = swd/S_toa |

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

The spin-up converges: domain-mean surface SMB settles at ~252 mm/yr by year
~25, while the ice-facing export ramps from negative toward the surface SMB as
the firn column fills — the expected approach to mass-throughput balance.

Against MAR v3.11 (the source model), final-year surface SMB:

| zone | chion | MAR |
|---|---|---|
| accumulation, z 2000–3500 m | +307 | +324 mm/yr |
| mid, z 1000–2000 m | +411 | +255 mm/yr |
| ablation margin, z 0–1000 m | +112 | **−600** mm/yr |
| net-ablating area fraction | **0 %** | 14.4 % |

Spatial pattern correlation 0.71. **The accumulation zone is captured well; the
whole discrepancy is under-ablation at the warm margins** — chion produces no
net-ablating cells, where MAR loses up to −3880 mm/yr. This is the classic
melt deficit of a snowpack forced with monthly-mean climatology and no
sub-monthly temperature variability. Levers, untried here: the BESSI diurnal
temperature cycle (`diurnal_temperature_cycle`), albedo tuning, or a temperature
variability term at the margins.

## Transmissivity finding

`diagnostics/transmissivity.jl` on GRL-16KM: observed τ = swd/S_toa averages
0.64, but ITM's default intercept (0.46) is too low (bias −0.06, negative R²)
and **elevation alone explains almost nothing** (R² 0.07). Cloud cover is the
real control — τ ≈ 0.81 + 7.4e-5·z_srf − 0.51·tcc raises R² to 0.32. A
cloud-aware transmissivity is the natural next parameterization. Consistent
with this, running the driver with `swd_source = transmissivity` (the ITM
clear-sky form) under-melts even more than the ERA5-shortwave run.
