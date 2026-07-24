# Changelog

All notable changes to chion are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[semantic](https://semver.org/spec/v2.0.0.html).

`v0.1.0` is the first tagged release. It covers the Chion.jl → Fortran port and
the CLIMBER-X SEMIX surface-scheme port on top of it; the history before this
tag is the port itself, summarised rather than enumerated.

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
