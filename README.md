# chion

Fortran snowpack model interface — a FESM-style static library (`libchion.a`)
with a clean public API and pluggable snowpack models (BESSI, PDD, ITM), ported
from Chion.jl.

The energy-balance model (BESSI) is configurable on orthogonal axes: column
structure (`Ntot`), longwave (`seb_scheme = bessi | semix`), turbulent
exchange (`turbulent_flux_scheme = bessi | semix | climberx`) and albedo
(`albedo_scheme = constant | dynamic | prescribed | semix | aging`), where
`seb_scheme = semix`, `climberx` and the `semix` albedo are ports of
CLIMBER-X's SEMIX surface scheme.

## Docs

- [CHANGELOG.md](CHANGELOG.md) — release history
- [docs/PLAN.md](docs/PLAN.md) — implementation plan and architecture
- [docs/steady_state_snowpack.md](docs/steady_state_snowpack.md) — standalone
  domain-driven spin-up (Greenland GRL-16KM/8KM from raw `ice_data`, no
  preprocessing), MAR validation, performance, and the transmissivity study
- [docs/porting_notes.md](docs/porting_notes.md) — decisions made porting from Chion.jl
- [docs/pdd_defects.md](docs/pdd_defects.md) — PDD model defect notes
- [docs/semix_port_scope.md](docs/semix_port_scope.md) — the SEMIX port: design,
  per-rung results, the composed-configuration study, and
  [what is left](docs/semix_port_scope.md#what-is-left)

## Reference

Ported from Chion.jl `a9ec154`, synced to Chion.jl `main` `9ec6cc7`
([docs/PLAN_dev_nils.md](docs/PLAN_dev_nils.md)). `validation/` gates the port
against it field by field; chion's deliberate deviations are listed in
`docs/porting_notes.md` and reverted by the `legacy_chion=1` build.

## Defaults

`input/chion_defaults.nml` holds Chion.jl's calibrated GrIS set (`03bb445`):
BESSI with `alpha_ice` 0.40, dynamic albedo, `seb_scheme =
turbulent_flux_scheme = "semix"` (exchange factor 2.5, stable coefficient 40),
cloud-proxy longwave, 5-layer thermal ice substrate, fine near-surface layers
(0.02, 0.05, 0.10, 0.30 m), 8 diurnal substeps with a 1 K cycle; PDD
`pdd_method = "simple"`. chion adds a thin-snow albedo blend
(`swe_crit_albedo` 10 kg m-2) and land columns (`H_ice = 0`: `alpha_land`, no
ice ablation, no substrate), so a host must fill `forc%H_ice` for BESSI. The
Greenland domain par files (`par/chion_grl16.nml`, `chion_grl8.nml`) take
`alpha_ice = 0.50`, calibrated against MAR with the blend.

## Build

`make` builds `libchion.a` (OpenMP by default; `make openmp=0` for serial;
the objects depend on a stamp of the `openmp` setting, so switching it rebuilds
them and `configme install`, which builds serial then OpenMP, leaves an OpenMP
library). `make tests` builds the acceptance tests, `make drivers`
`chion_column.x` and `chion_grid.x`, `make all` both. Options (`make usage`):

- `precision=dp`: `wp = dp` (`libchion/*-dp`), for validation;
- `legacy_chion=1`: Chion.jl's values for chion's deliberate corrections
  (`*-legacy`), for validation only;
- `fpsafe=1`: value-safe `-O2` (`-fp-model precise` for ifx) in
  `libchion/*-fpsafe`. Build and run the acceptance tests and `validation/`
  with it (`make all fpsafe=1 [precision=dp [legacy_chion=1]]`); a machine's
  `-Ofast` breaks their one-ulp checks. Production keeps the machine's flags;
- `debug=1` (bounds and traps), `debug=2` (profiling).

`chion_grid.x` in domain mode needs a stack of ~40 MB for fesm-utils' ERA5
regridding (`ulimit -s unlimited`; see `docs/steady_state_snowpack.md`).
