# Steady-state snowpack on an ice-sheet domain

A standalone way to spin a chion snowpack up to a repeating seasonal state on a
real domain. Greenland (`GRL-16KM/8KM`, MAR v3.11 + ERA5, no preprocessing step)
and Antarctica (`ANT-32KM/16KM/8KM`, RACMO2.4/ANT-12 CORDEX, one CDO
preprocessing step) are both wired; other domains are one loader subroutine
away.

## Pieces

| piece | where | role |
|---|---|---|
| `chion_forcing_monthly` | `src/` (library) | mean-preserving monthly→daily interpolation (360-day, 12×30); the daily series averages back to the input monthly means, so precip totals are conserved |
| `insolation` | `libs/insol/` (driver layer) | top-of-atmosphere daily insolation (Laskar LA2004); kept out of `libchion.a`, preserving "the host supplies insolation" |
| `chion_domain` | `libs/domains/` (driver layer) | assembles standardized SI monthly forcing per domain: Greenland = MAR v3.11 (t2m, sf, rf, smb/melt/runoff) native on grid + ERA5 shortwave/cloud **conservatively regridded** (coords `map_init(con)`, cache under `maps/`); Antarctica = RACMO2.4/ANT-12 pre-regridded to the ANT grid (CDO, see below) + BedMachine geometry. TOA insolation from latitude |
| `chion_grid.x` | `tests/chion_grid.f90` | one driver, `forcing_source = file | domain`; the domain path cycles the annual climatology for `n_years`, with `swd_source = file | transmissivity | transmissivity_seasonal` |
| `diagnostics/*.jl` | `diagnostics/` | Julia + CairoMakie analysis (own project env): `compare_smb.jl` (chion vs MAR SMB figure), `transmissivity.jl` (ITM τ vs τ_obs = swd/S_toa) |

## Run

```
cd <run-dir>
ln -s <chion>/input input                 # chion reads input/chion_defaults.nml
<chion>/libchion/bin/chion_grid.x <chion>/par/chion_grl16.nml
```

Edit `path_ice_data` / `path_insol` in the par file. The conservative ERA5
weight map is generated on the first run and cached under `maps/` (gitignored,
regenerated if absent). A 50-year GRL-16KM BESSI run (7204 columns) is ~47 s on
one core; the build is OpenMP by default, so `OMP_NUM_THREADS=8` cuts it to
~17 s (see Performance).

### Antarctica (RACMO2.4 / ANT-12)

Antarctica needs a one-time preprocessing step because the source is a 30-year
CORDEX monthly *time series* on a *rotated-pole* grid, not a ready-made
climatology. Two scripts build it:

```
scripts/download_racmo_ant12.sh        # CORDEX tas,pr,rsds,clt 1981-2010 -> ~/data/racmo (~1 GB)
scripts/build_ant12_climatology.sh     # 12-month climatology, SI units, regridded to ANT-32/16/8KM
```

`build_ant12_climatology.sh` averages to 12 months, converts units, and
conservatively regrids (`cdo remapcon`) onto each ANT grid. The regrid is done
here rather than in-model because coords cannot remap from a rotated-pole
source ([fesmc/fesm-utils#8](https://github.com/fesmc/fesm-utils/issues/8)).
`load_antarctica` then reads the per-grid climatology plus BedMachine geometry;
set `path_racmo` (and `path_ice_data` for BedMachine) in the par file:

```
<chion>/libchion/bin/chion_grid.x <chion>/par/chion_ant32.nml
```

The CORDEX set is atmospheric only: no reference SMB (chion computes it), and
total precip is split into snow/rain at the freezing point.

## Result (GRL-16KM, BESSI, 50 yr, ERA5 shortwave)

> This section establishes how SMB is measured and the **baseline** (no diurnal
> substepping) behaviour, whose elevation-dependent under-ablation motivates the
> two sections that follow. The **current default** (diurnal substepping, gate
> −1 °C) improves the domain bias from +68 to +18 mm/yr and R² from 0.74 to
> 0.84 — see *Diurnal shortwave substepping*.

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

OpenMP is the default build (parallel over columns; `make openmp=0` for serial;
runtime threads via `OMP_NUM_THREADS`). On a 4-performance-core machine:

| threads | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| GRL-16KM wall | 46.9 s | 30.0 s | 20.9 s | 17.3 s |
| speedup | 1.0× | 1.6× | 2.2× | 2.7× |

GRL-8KM at 8 threads: 74 s (from 202 s, 2.7×) — same scaling. Returns diminish
past the 4 performance cores (the machine's other 6 are efficiency cores).

### Antarctica run times

The whole-Antarctica domain (grounded ice + shelves) carries far more ice
columns than Greenland at the same resolution, so absolute wall times are
larger even though per-column cost is unchanged. Measured on the same machine,
50-yr BESSI runs (wall time = run loop only; the CDO regrid is preprocessing,
not counted):

| grid | ice columns | serial | 8 threads | per col-step (8 thr) |
|---|---|---|---|---|
| ANT-32KM | 13 212 | 67 s | 30 s | 1.27e-7 s |
| ANT-16KM | 52 859 | — | 96 s | 1.01e-7 s |
| ANT-8KM | 211 070 | — | 393 s | 1.03e-7 s |

Per-column-step matches Greenland (ANT-32KM serial 2.84e-7 s vs GRL ~3.8e-7),
so throughput is not the issue — **ANT-8KM is slow because it has 211 k columns,
7× GRL-8KM's 29 k** (Antarctica's ice area is ~7× Greenland's). The 393 s is
domain size, not a per-cell regression; the same per-column engine drives both.
(Serial figures are shown only where measured; the others were run at the
default 8 threads.)

Physics is resolution-robust: GRL-8KM gives bias +15, RMSE 212, R² 0.84 —
statistically identical to GRL-16KM (+18, 0.84). The SMB pattern is set by
elevation and temperature, both resolved at 16 km; refining sharpens the margin
gradient but does not move the domain-scale skill.

## Model and layer-count comparison

How fast can BESSI be pushed, and where do the cheaper models break? A 50-year
GRL-16KM spin-up (identical ERA5 forcing, serial/1 core) for BESSI at
`Ntot = 15/7/5/3/2/1` and for PDD and ITM. Everything but `model` and `Ntot` is
the GRL-16KM default (diurnal substep, gate −1 °C). Reproduce with
`scripts/run_layer_comparison.sh` and
`diagnostics/compare_models.jl` (writes `output/cmp_smb_stats.csv`).

**Speed.**

| config | wall [s] | per col-step | vs BESSI 15 |
|---|---|---|---|
| BESSI n=15 | 48.2 | 3.72e-7 | 1.0× |
| BESSI n=7  | 36.6 | 2.82e-7 | 1.3× |
| BESSI n=5  | 32.7 | 2.52e-7 | 1.5× |
| BESSI n=3  | 28.4 | 2.19e-7 | 1.7× |
| BESSI n=2  | 25.4 | 1.96e-7 | 1.9× |
| BESSI n=1  | 22.5 | 1.73e-7 | 2.1× |
| PDD        |  3.7 | 2.88e-8 | 12.9× |
| ITM        |  8.4 | 6.47e-8 | 5.7× |

Cutting the BESSI layer cap barely helps: even all the way to a single layer is
only **2.1×**, not 15×. The per-column layer loop (energy balance, densification,
percolation) is not the whole cost — accumulation, surface fluxes and the diurnal
substep are fixed overhead per column, so dropping layers hits diminishing
returns fast. The real speed comes from dropping the layered column *and* the
energy balance: PDD (bulk degree-day) is **13×** faster and ITM
(insolation–temperature) **6×**.

**BESSI layer count barely touches the SMB — and a single layer is best.**
Domain skill is flat, in fact marginally *better* toward fewer layers:

| BESSI Ntot | bias | RMSE | R² |
|---|---|---|---|
| 15 | +18 | 214 | 0.84 |
| 7  | +18 | 213 | 0.84 |
| 5  | +16 | 210 | 0.84 |
| 3  | +14 | 208 | 0.85 |
| 2  | −0.5 | 197 | 0.86 |
| 1  | −2.3 | 197 | 0.86 |

Every elevation band is statistically identical across `Ntot` (interior
bit-for-bit: R² 0.99, RMSE 19; margins −640→−646). SMB is a column-*integrated*
quantity — storage tendency + ice export — and the depth cap that bounds the
firn column is the hard-coded 15-layer reference depth **regardless of `Ntot`**
(see `enforce_snow_depth_cap`), so even one layer holds the same total firn, as a
single thick layer. The vertical structure that extra layers resolve (thermal
gradient, refreeze profile) does not move the *annual* mass balance at domain
scale; the single evolving-density layer even lands closest to MAR (bias ≈ 0).
`Ntot=1` uses chion's dedicated single-layer closed-form energy solve, so it is a
genuine single-layer energy-balance column, not a degenerate multilayer run.
**On Greenland SMB, BESSI n=1 is the fastest config (2.1×) AND the most accurate;
the layers earn their keep only if you need the firn-column state itself**
(temperature, density, refreeze depth), not just SMB.

**PDD and ITM trade skill for speed, and exactly where the physics lives.** All
three models are near-perfect in the interior (pure accumulation, no melt: R²
0.99) and good in the mid zone. They diverge below ~1500 m, in the percolation
and margin zones where melt–refreeze detail decides the balance:

| config | all-ice R² | RMSE | margin R² | lower-zone R² |
|---|---|---|---|---|
| BESSI (any n) | 0.84–0.85 | 208–214 | 0.58 | 0.69–0.72 |
| PDD | 0.61 | 331 | −0.04 | 0.18 |
| ITM | 0.24 | 463 | −1.24 | −0.47 |

PDD keeps the domain bias small (+22) but loses the *pattern* in the melt zones
(lower-zone R² 0.69→0.18, margins go negative). ITM is worse on both counts:
with its default melt coefficients it massively **under-ablates the margins**
(chion +273 vs MAR −740; bias +1013 there), dragging the domain bias to +172 and
R² to 0.24 — the same margin miscalibration flagged in *Transmissivity finding*,
here isolated to the melt formulation since ITM is driven by the ERA5 shortwave.

**Bottom line.** For SMB skill, keep BESSI but drop to `Ntot = 5` (1.5× free,
no measurable cost). If ~13× throughput matters more than margin skill, PDD is
the honest cheap option — its bias is small and it is only wrong where every
bulk model is wrong. ITM needs recalibration before it competes. Figures:
`output/cmp_<config>/smb_chion_vs_mar.png`.

## Transmissivity finding

`diagnostics/transmissivity.jl` on GRL-16KM: observed τ = swd/S_toa averages
0.64, but ITM's default intercept (0.46) is too low (bias −0.06, negative R²)
and **elevation alone explains almost nothing** (R² 0.07). Cloud cover is the
real control. Running the driver with `swd_source = transmissivity` (the ITM
clear-sky form) under-melts even more than the ERA5-shortwave run.

**Seasonal structure — the annual fit is misleading.** Fitting each month on
its own (`tau_fit_seasonal.png`):

| | Mar | May | Jun | Jul | Aug | Sep |
|---|---|---|---|---|---|---|
| mean τ | 0.64 | 0.72 | 0.72 | 0.68 | 0.64 | 0.58 |
| cloud coef c | −0.04 | −0.15 | −0.24 | −0.29 | −0.31 | −0.11 |
| R² | 0.59 | 0.82 | 0.86 | 0.82 | 0.73 | 0.50 |

τ peaks in **May**, and the **cloud coefficient strengthens ~3× from spring to
late summer** (clouds suppress τ far more in Jul–Aug). The pooled annual fit
gives c=−0.51 — *steeper than any single month* — because season is a
confounder (autumn is low-τ and cloudy, spring high-τ), a Simpson's-paradox
artifact. So a season-resolved τ is the honest form, and ITM should be
calibrated on the **melt-season** relationship, not the annual pooled one.

### The seasonal transmissivity reproduces the ERA5 shortwave

`swd_source = "transmissivity_seasonal"` uses the JJA-pooled cloud-aware fit,
**τ = 0.768 + 6.99e-5·z_srf − 0.376·tcc** (clamped to [0,1], × S_toa; R²=0.67
within JJA). A 50-year GRL-16KM run driven this way — with NO ERA5 shortwave
data, only insolation, elevation and cloud cover — matches the ERA5-driven run:

| swd_source | runoff | SMB | bias | R² | net-abl. |
|---|---|---|---|---|---|
| MAR | 191 | 190 | — | — | 14.4 % |
| ERA5 file (default) | 178 | 209 | +18 | 0.84 | 15.5 % |
| transmissivity_seasonal | 170 | 216 | +26 | 0.82 | 15.2 % |

Band-by-band it tracks the ERA5 run (margins −595 vs −640, interior identical).
The cloud term is what makes it work: the clear-sky `transmissivity` form
(τ = 0.46 + 6e-5·z, no cloud) under-melts badly. This is the practical payoff of
the diagnostic — a shortwave parameterization good enough to stand in for the
ERA5 field, needing only cloud cover (more widely available, and estimable for
paleo domains). Still a single fit, not truly month-varying; a month-resolved
(a,b,c) is the next refinement.
