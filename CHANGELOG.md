# Changelog

All notable changes to chion are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[semantic](https://semver.org/spec/v2.0.0.html).

`v0.1.0` is the first tagged release. It covers the Chion.jl → Fortran port and
the CLIMBER-X SEMIX surface-scheme port on top of it; the history before this
tag is the port itself, summarised rather than enumerated.

## [Unreleased]

### Added

- BESSI `albedo_scheme = "aging"` (Chion.jl `6d06af6`): snowfall resets to
  `alpha_dry`, exponential relaxation to `alpha_wet` (`aging_cold_timescale_days`
  20 d, `aging_melting_timescale_days` 2 d); `snow_age_days` in output and
  restart. Default scheme unchanged.
- Aging albedo: snowfall rejuvenates in proportion to its mass,
  `f = min(1, S/aging_snowfall_ref)` (10 kg m-2), instead of Chion.jl's reset
  on any snowfall (D30; reverted under `legacy_chion`).
- `chion_get_surface_flux_totals`: cumulative melt/runoff/refrz/subl [kg m-2],
  so a host aggregating over many steps differences two calls.
- BESSI diurnal T amplitude rises with elevation, `clamp(A0 + γ·max(z − z_ref, 0),
  0, A_max)` (Chion.jl `d0146e1`): `diurnal_temperature_amplitude_gradient`
  [K km-1], `_reference_height` [m], `_max` [K]; neutral defaults (bit-identical).
  `surface_height` is now packed into the step forcing; a missing (NaN) height
  adds no excess.
- BESSI thermal ice substrate (Chion.jl `03bb445`): `&bessi ice_substrate_layers`
  (0 = none, the default here; upstream 5) and `ice_substrate_top_thickness`
  (0.05 m, doubling downward), insulated base, one implicit solve with the snow;
  bare ice carries cold content and must warm to T0 before melting. Only where
  `H_ice > 0` (now packed into the step forcing for BESSI); reset with the
  column; restart field `ice_temperature`, older restarts start at
  `min(t_srf, T0)` (D34). Default runs bit-identical.
- BESSI fine near-surface layers (Chion.jl `03bb445`): `&bessi
  near_surface_layer_max_thicknesses` [m] holds layers 1-4 at fixed thicknesses
  (upstream 0.02, 0.05, 0.10, 0.30) by a conservative remesh after accumulation
  and after refreezing; the `mass_min` surface merge is off with layer 1
  limited. 0 = no limit, the default here (D36). Default runs bit-identical.
- Fine near-surface layers: the first layer below them is split and merged by
  mass like the surface layer (D32; reverted under `legacy_chion`). Chion.jl
  keeps everything below the fine layers in one unsplit layer up to the 22.5 m
  depth cap. 10-yr `chion_column` example with fine layers, Chion.jl's
  behaviour -> chion: mean layer count 4.9 -> 10.0 (max 5 -> 15), runoff
  +7.7 %, refreezing -6.0 %, melt +0.6 %, final thickness -12.6 %, liquid water
  -40 %; runoff and refreezing within 2.5 % of the run without fine layers.

- BESSI cloud-proxy longwave (Chion.jl `03bb445`): `&chion_const longwave_scheme =
  "graybody" | "cloud_proxy"` (default `graybody` here; upstream `cloud_proxy`) and the
  six `lw_*` coefficients; emissivity from air temperature and a shortwave cloudiness
  proxy against a daily TOA, resolved once per step before the diurnal substeps.
  Optional host field `forc%toa_shortwave` (`has_toa_shortwave`) replaces chion's
  fixed-orbit TOA (D33). Default runs bit-identical.
- `calendar_solar_longitude_deg(day_of_year)` (re-exported by `chion`): Chion.jl's
  calendar-day solar longitude, for hosts and drivers without an orbital one.

### Changed

- `seb_scheme` selects the longwave only (`bessi` | `semix`); new
  `turbulent_flux_scheme` (`bessi` | `climberx`) selects the sensible and latent
  heat (Chion.jl `d0146e1`). CLIMBER-X SEMIX's aerodynamic exchange is
  `turbulent_flux_scheme = "climberx"`: the former `seb_scheme = "semix"` is
  `seb_scheme = "semix"` + `turbulent_flux_scheme = "climberx"` (bit-identical);
  `semix_qsat` is renamed `climberx_qsat` (`"climberx"` | `"bessi"`) (D37).

### Fixed

- `chion_column.x` and `chion_grid.x` derive the solar longitude with Chion.jl's
  calendar formula; it was `360*(doy-1)/year_length`, i.e. 0 (the March equinox)
  on 1 January. Affects only diurnal substeps, SEMIX `coszm` and the cloud proxy.

- BESSI dynamic albedo: aging scaled by `dt_days` (Chion.jl `6d077c5`), so
  diurnal substeps no longer age the albedo once each (upstream defect 19).
  Daily runs bit-identical.
- ITM per-step `tsrf`: `melt_net` scaled to the annual rate `firn_fac` is
  calibrated on (360-day year); firn warming was ~360x too small (D27).
- BESSI vapour exchange: `vapor_mass`/`sublimation` report the mass actually
  removed when sublimation exhausts the surface layer (Chion.jl `03bb445`;
  upstream defect 1). Diagnostics only.
- BESSI fresh snow on a bare column takes the air temperature in every new
  layer, not only layer 1 (Chion.jl `03bb445`).

### Changed

- BESSI surface boundary is Robin (Chion.jl `03bb445`): `Tsrf` is the snow-air
  interface temperature, eliminated through the top cell's half-thickness
  conductance `2K1/dz1`; no surface heat capacity, no one-layer closed form;
  melt energy `Q(T0) − Gs(T0 − T1)`; fluxes linearized and the vapour flux
  evaluated at `Tsrf`; melt sets `Tsrf = T0`. Changes default BESSI results
  (and `seb_scheme = "semix"`, which shares the solver).
- BESSI dynamic albedo wetness relaxation `α_wet + (α − α_wet)(1 − r)^dt`,
  `r = clamp(lwc/max_lwc_albedo, 0, 1)` (Chion.jl `03bb445`): substep-invariant;
  same as the linear law at daily steps up to round-off.
- BESSI depth cap is a constant 22.5 m (Chion.jl `03bb445`), no longer
  `15*mass_split*1.5/300`; identical at the default `mass_split`.
- Scheme names: aliases removed, canonical names only (Chion.jl `03bb445`):
  `semix_snow_albedo = "warren_wiscombe"` (was `"ww"`/`"warren"`); albedo
  `"bessi"`/`"legacy"`, fresh-snow density `"bessi"`/`"htessel"`, PDD
  `"calov_greve"`, `semix_qsat` `"climberx"`/`"chion"` are rejected.
- BESSI heat conduction: interface conductance is harmonic (half-layer
  resistances in series), `2KiKj/(Kj dzi + Ki dzj)`, with `beta = -dt/ci`
  (Chion.jl `81034fa`, `d0146e1`). Identical on uniform columns; changes
  default BESSI results where conductivity jumps between layers.
- BESSI snow thermal conductivity: Calonne et al. (2019), density- and
  temperature-dependent (Chion.jl `49990e6`), replacing `Ki*(rho/1000)^1.88`;
  `Ki` removed from `&chion_const`. Changes default BESSI results (K +14% at
  500, +23% at 917 kg m-3).
- BESSI turbulent latent heat (Chion.jl `d0146e1`): `Lv` at a melting snow
  surface, `Lv+Lm` below `T0` and on bare ice; the parameterized vapour mass is
  the humidity-gradient mass flux, independent of the latent heat (prescribed
  `q_lh` still converts with the phase's `L`). Changes default BESSI results
  when humidity forcing is on (latent flux ~12% smaller at melting surfaces).
- Layer merges mix density and temperature as `x1 + w2*(x2 - x1)`: equal
  values stay exact (two layers at T0 stay at T0). Round-off-level change (D31).
- Shared physical constants (`rho_i`, `rho_w`, `ci`, `cw`, `Lm`, `grav`, `T0`)
  come from fesm-utils `phys_const_class`: `chion_init(..., cnst=)` takes the
  host's record, otherwise chion loads `phys_const_file` (now the `phys_const`
  schema, Chion.jl values). chion's own constants moved to `&chion_const` in
  `input/chion_defaults.nml` (sparse overrides). `seconds_per_day` removed
  (`phys_constants:sec_day`). Standalone output unchanged (D28).
- `chion_get_surface`: outputs optional; parallel over columns (including the
  MV fill of inactive ones), no full-array copies, so it is cheap enough to
  call every step.
- `chion_update` snapshots (for `chion_get_smb`/`chion_get_surface_fluxes`)
  only active columns, in parallel (`chion_model_cum_active`); was a serial
  copy over all columns every step.
- `itm_par_load` takes optional `defaults_file`/`defaults_group`; `&itm` may now
  be sparse or absent, like `&bessi` and `&pdd` (needed to share smbpal's `&itm`
  group in yelmox).
- Build: `libchion.a` is a file target; objects depend on `libfesmutils.a`.
- Build: `fpsafe=1` (value-safe `-O2`, `libchion/*-fpsafe`) for the acceptance
  tests and validation/; make creates every flavour's build directories.
- PDD: default `pdd_method = "simple"` (was `"pism"`), matching Chion.jl;
  set `"pism"` explicitly for monthly steps or smbpal-like melt. validation/
  gates PDD against Chion.jl for both methods.
- ITM: `legacy_chion=1` reverts the D27 `tsrf` scaling (Chion.jl's daily
  rate); validation/ gates ITM against Chion.jl's `ITMModel`. `chion_grid.x`
  file forcing reads ITM ice thickness and annual PDDs (`name_hice`,
  `name_pdds`; `"None"` = 0). Production ITM unchanged.

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
