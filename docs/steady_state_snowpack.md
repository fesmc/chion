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

The spin-up converges: the domain-mean ice-facing SMB (`smb_ice`, "net mass
flux to the ice sheet") settles at ~200 mm/yr by year ~40 as the firn column
fills — the expected approach to mass-throughput balance.

**Which chion field is the MAR-comparable SMB.** MAR `smb` is the ice-facing
surface mass balance. The matching chion field is `smb_ice`, NOT the column
storage tendency d(snow+mass_base)/dt: on a cold start from an empty snowpack
every column is still building firn toward equilibrium, so the storage tendency
is positive almost everywhere and hides the ablation. `smb_ice` is what BESSI
drives negative through its bare-ice ablation branch once a column melts out.

Against MAR v3.11, final-year `smb_ice`:

| zone | chion | MAR |
|---|---|---|
| accumulation, z 2000–3500 m | +252 | +324 mm/yr |
| mid, z 1000–2000 m | +325 | +255 mm/yr |
| ablation margin, z 0–1000 m | −359 | −600 mm/yr |
| domain mean | +200 | +190 mm/yr |
| net-ablating area fraction | 12.4 % | 14.4 % |
| min (margin) | −2833 | −3881 mm/yr |

Spatial pattern correlation 0.86; margin runoff 945 vs 1178 mm/yr. chion tracks
MAR across the whole elevation range and reproduces the ablation zone (down to
−2833 mm/yr over 12.4 % of the ice area). The residual is a modest
**under-ablation at the very warmest margin cells** (−359 vs −600), the expected
melt deficit of monthly-mean climatological forcing with no sub-monthly
temperature variability. Levers, untried here: the BESSI diurnal temperature
cycle (`diurnal_temperature_cycle`), albedo tuning, or a margin temperature
variability term.

## Transmissivity finding

`diagnostics/transmissivity.jl` on GRL-16KM: observed τ = swd/S_toa averages
0.64, but ITM's default intercept (0.46) is too low (bias −0.06, negative R²)
and **elevation alone explains almost nothing** (R² 0.07). Cloud cover is the
real control — τ ≈ 0.81 + 7.4e-5·z_srf − 0.51·tcc raises R² to 0.32. A
cloud-aware transmissivity is the natural next parameterization. Consistent
with this, running the driver with `swd_source = transmissivity` (the ITM
clear-sky form) under-melts even more than the ERA5-shortwave run.
