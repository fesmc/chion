# Changelog

All notable changes to chion are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[semantic](https://semver.org/spec/v2.0.0.html).

`v0.1.0` is the first tagged release. It covers the Chion.jl → Fortran port and
the CLIMBER-X SEMIX surface-scheme port on top of it; the history before this
tag is the port itself, summarised rather than enumerated.

## [Unreleased]

Synced with Chion.jl `main` `9ec6cc7` (= `03bb445`; `docs/PLAN_dev_nils.md`, deviations
`docs/porting_notes.md` D29-D43). Default BESSI results change; the validation gate is green
against `9ec6cc7` in twelve BESSI configurations, PDD and ITM.

### Breaking

- Namelists (re-sync host copies of `input/chion_defaults.nml`):
  - `Ki` removed from `&chion_const` (Calonne conductivity); a par file that sets it is
    not read for it.
  - Scheme aliases removed, canonical names only: albedo `bessi`/`legacy`, fresh-snow
    density `bessi`/`htessel`, PDD `calov_greve`, `semix_snow_albedo` `ww`/`warren` (now
    `warren_wiscombe`) are rejected.
  - `seb_scheme` selects the longwave only (`bessi` | `semix`); the turbulence is
    `turbulent_flux_scheme` (`bessi` | `semix` | `climberx`). chion's CLIMBER-X exchange is
    `"climberx"`: the former `seb_scheme = "semix"` is `seb_scheme = "semix"` +
    `turbulent_flux_scheme = "climberx"` (bit-identical); `"semix"` turbulence is now
    Chion.jl's bulk scheme. `semix_qsat` is renamed `climberx_qsat` (D37).
  - `diurnal_temperature_amplitude_max` (default 1 K) caps the amplitude; a larger
    `diurnal_temperature_amplitude` needs it raised.
  - `chion_grid.x` file mode reads `&ctrl name_hice`, `name_pdds` (`"None"` = 0).
- BESSI `H_ice = 0` is a land column (no ice substrate, no ice ablation, `alpha_land`):
  a host must fill `forc%H_ice` for BESSI, as for ITM (D34, D41).
- ITM output takes Chion.jl's names and units (D42): `alb_s`, rates `smb` (total SMB),
  `smbi`, `melt`, `runoff`, `refreezing`, `melt_net` [mmWE day-1], cumulative `smb_cum`,
  `smb_ice`, `melt_cum`, `runoff_cum`, `refreezing_cum` (were `albedo`, `smb_total`,
  cumulative `melt`/`runoff`/`refreezing`, ice-facing `smb` in kg m-2 s-1). Restarts
  unchanged.

### Added

- BESSI thermal ice substrate (`&bessi ice_substrate_layers`, `ice_substrate_top_thickness`;
  restart `ice_temperature`, older restarts start at `min(t_srf, T0)`; D34).
- Fine near-surface layers (`&bessi near_surface_layer_max_thicknesses`, 0 = no limit,
  D36); the first layer below them is split and merged by mass (D32).
- Cloud-proxy longwave (`&chion_const longwave_scheme`, `lw_*`); optional host
  `forc%toa_shortwave` (D33).
- `turbulent_flux_scheme = "semix"`: Chion.jl's bulk turbulence (`semix_*` parameters);
  bare ice with `Lv + Lm` (D35) and `R_dry` (D38).
- Albedo `"aging"` (`aging_*_timescale_days`; `snow_age_days` in output and restart), with
  a snowfall refresh `1 - exp(-S/aging_snowfall_ref)` (D30).
- Thin-snow albedo: snow albedo blended with the background by `min(1, SWE/swe_crit_albedo)`
  (CLIMBER-X `tanh` form under `albedo_scheme = "semix"`; `c_fsnow`, `c_fsnow_orog`); new
  state `albedo_snow`; land background `alpha_land` (D40, D41).
- Elevation-dependent diurnal T amplitude (`diurnal_temperature_amplitude_gradient`,
  `_reference_height`, `_max`); `surface_height` and `H_ice` packed into the step forcing.
- BESSI output `ice_temperature`; `calendar_solar_longitude_deg`;
  `chion_get_surface_flux_totals` (cumulative melt/runoff/refrz/subl).
- Build: `fpsafe=1` (value-safe `-O2`, `libchion/*-fpsafe`) for the tests and validation/;
  make creates every flavour's directories.
- validation/: PDD (both methods) and ITM gated against Chion.jl; BESSI configurations
  for each switchable `03bb445` option and for the defaults.

### Changed

- Defaults are Chion.jl `03bb445`'s calibrated GrIS set: `alpha_ice` 0.40, `seb_scheme =
  turbulent_flux_scheme = "semix"` (2.5, 40), `longwave_scheme = "cloud_proxy"`, 5 substrate
  layers, fine layers (0.02, 0.05, 0.10, 0.30 m), 8 diurnal substeps with a 1 K cycle; plus
  chion's thin-snow blend (`swe_crit_albedo` 10 kg m-2). PDD `pdd_method = "simple"`.
- BESSI physics from Chion.jl (unswitchable): Robin surface boundary (`Tsrf` = interface
  temperature), harmonic interface conductance, Calonne et al. (2019) conductivity,
  phase-dependent latent heat with a gradient-based vapour mass, wetness relaxation
  `(1-r)^dt`, depth cap 22.5 m, aging × `dt_days`.
- Default results (all of the above): 10-yr `chion_column` melt +27 %, refreezing +78 %;
  GRL-16KM vs MAR (Gt/yr; MAR SMB 348, melt 518, runoff 349): SMB 377 -> 279, melt 374 ->
  500, runoff 329 -> 427, R² 0.84 -> 0.77. Cost per step ~2x at 16 threads (GRL-16KM 50 yr
  20 -> 43 s, shared node), ~4x serial.
- Domain calibration: `par/chion_grl16.nml`, `chion_grl8.nml` take `alpha_ice = 0.50`
  (MAR grid 0.40-0.55 x `swe_crit_albedo` 5-40): GRL-16KM SMB 340, melt 440, runoff 367,
  R² 0.86 (baseline 0.84); GRL-8KM SMB 343 (MAR 358). Library defaults stay Chion.jl's.
- Shared physical constants from fesm-utils `phys_const_class` (`chion_init(..., cnst=)`);
  chion's own constants in `&chion_const`; `seconds_per_day` removed (D28).
- Performance: dynamic OpenMP schedule, solar geometry once per column-day
  (bit-identical); `chion_get_surface` and `chion_update` snapshots over active columns
  only, in parallel.
- Performance (PERF): BESSI ~32 % faster serial, ~28 % at 16 threads (GRL-16KM 50 yr,
  exclusive node: 549 -> 371 s, 40.4 -> 29.1 s): vectorised conductivities and no
  per-row cross-module calls in the energy solve, water content from the pore volume in
  percolation, no sub-ulp near-surface remesh transfers (D43), vectorised densification
  Arrhenius factors, solar declination once per step and substep sines once per
  column-day, semix neutral exchange coefficients as derived constants
  (`chion_const_derive`), and `src/physics` compiled as one translation unit
  (`chion_physics.f90` includes the modules; cross-module inlining without `-ipo`).
  Results move at the round-off noise floor (as for a recompile).
- Layer merges mix as `x1 + w2*(x2 - x1)` (exact for equal values; D31).
- `itm_par_load` takes optional `defaults_file`/`defaults_group` (sparse `&itm`).
- `legacy_chion=1` reverts the deliberate corrections to Chion.jl's values for validation
  (D24 list; gas constant 8.314, ITM `tsrf` scaling D27).

### Fixed

- Rain on a bare column runs off exactly once (Chion.jl dropped it, then double-counted
  it; D29).
- ITM per-step `tsrf`: `melt_net` scaled to the annual rate `firn_fac` is calibrated on
  (firn warming was ~360x too small; D27).
- BESSI vapour diagnostics report the mass actually removed (Chion.jl `03bb445`); fresh
  snow on a bare column takes the air temperature in every new layer.
- Diurnal substeps: a day not split keeps its forcing (no polar-night shortwave loss; D39).
- `chion_grid.x` (domain) drives ITM with the TOA insolation (it got the surface
  shortwave, attenuated twice: GRL-16KM SMB 663 -> 480, melt 62 -> 271 Gt/yr).
- Drivers use Chion.jl's calendar solar longitude (was 0 on 1 January).
- Build: objects depend on an `openmp` stamp, so switching it rebuilds them; a fresh
  `configme install` leaves an OpenMP `libchion.a`. `libchion.a` is a file target; missing
  module dependencies added (`make -j`).
- `chion_grid.x` reads `trans_a/b/c` only for the `swd_source` that uses them.

### Removed

- `Ki`, the scheme aliases, `semix_qsat` (see Breaking); `seconds_per_day`.

## [v0.1.0] — 2026-07-24

First tagged release: a complete Fortran port of Chion.jl, plus SEMIX's surface
energy balance and spectral/dust albedo as selectable, orthogonal options.

### Snowpack models

- **BESSI** — layered firn column: accumulation, layer splitting/merging,
  implicit conduction with a linearized surface energy balance and a
  melting-point re-solve, melt, percolation, refreezing, densification (BESSI
  and HTESSEL branches), and surface vapour exchange.
- **PDD** and **ITM** — bulk melt parameterizations, no energy balance.
- Precision policy: `wp = sp` for state and interfaces (`precision=dp` builds
  for reference comparison), `wp_acc = dp` mandatory for cumulative
  accumulators.

### SEMIX surface scheme (CLIMBER-X `src/smb/`; Willeit, Calov, Ganopolski)

Selectable on orthogonal axes — column structure (`Ntot`), surface energy
balance (`seb_scheme`) and albedo (`albedo_scheme`) vary independently. See
[docs/semix_port_scope.md](docs/semix_port_scope.md).

- **`albedo_scheme = "semix"`** — Warren & Wiscombe 1980 and Dang et al. 2015
  spectral snow albedo in four bands ({vis,nir}×{dir,dif}), grain-size aging,
  dust-in-snow darkening with seasonal-max-SWE melt amplification, optional
  host-supplied bare-ice albedo. Bands are collapsed to broadband using the
  incoming-SW spectral weights.
- **`seb_scheme = "semix"`** — CLIMBER-X's aerodynamic surface energy balance,
  dispatched at all three flux sites (the linearized surface row, the exact
  bare-ice fluxes, and the post-solve vapour mass):
  - `resistance`: snow-weighted roughness, log-law neutral exchange,
    bulk-Richardson stability
  - `ebal` sensible, latent and longwave coefficients, mapped onto chion's
    `q_const`/`q_lin` (coupling decision α — no skin node introduced)
  - saturation humidity selectable at runtime (`semix_qsat`), dew inhibition
    (`l_dew`), forced-neutral option (`l_neutral`)
  - emissivity-weighted downwelling longwave, with separate snow and ice
    emissivities

### Added

- `snow_vapor` module — the vapour-pressure parameterizations, extracted from
  `snow_surface_fluxes` so the SEMIX scheme can share them without a module
  cycle.
- `snow_seb_semix` module and `test_seb.x` acceptance test.
- Grid-driver knobs `rh_default`, `dust_dep_default` (both off at zero).
- `scripts/run_semix_matrix.sh` — the albedo × SEB × humidity configuration
  matrix on GRL-16KM.
- Domain loaders for Greenland (MAR/ERA5) and Antarctica (RACMO2.4).

### Validation

GRL-16KM, 50 yr, surface SMB vs MAR, `Ntot=1`, all ice:

| albedo | SEB | humidity | bias | RMSE | R² |
|---|---|---|---|---|---|
| dynamic | bessi | off | −2.3 | 197 | 0.86 |
| semix | bessi | off | −9.4 | 199 | 0.86 |
| dynamic | semix | off | +20.0 | 211 | 0.84 |
| semix | semix | off | +17.1 | 213 | 0.84 |
| semix | semix | 0.7 | +0.9 | 241 | 0.79 |

With comparable humidity forcing the full SEMIX configuration is
indistinguishable from BESSI (R² 0.83 both). The apparent gap at equal
`rh_default` is a water-vs-ice saturation-reference artifact, not physics — see
docs. 15/15 acceptance tests pass; `seb_scheme = "bessi"` output is
bit-identical to the pre-SEMIX baselines at every layer count.

### Not included

- **Rung 4**, SEMIX's `tstd`/Krapp-2016 statistical diurnal melt — decided
  against: chion's diurnal shortwave substepping already resolves sub-daily
  melt explicitly, the scheme is off by default in CLIMBER-X itself, and it is
  the one part of SEMIX that coupling α cannot map cleanly.
- Spectral net shortwave (currently a broadband collapse), SEMIX's continuous
  `f_snow` snow-cover blend, and a real humidity field. See
  [What is left](docs/semix_port_scope.md#what-is-left).
- CLIMBER-X integration — deliberately deferred; offline physics first.
