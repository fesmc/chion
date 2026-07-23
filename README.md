# chion

Fortran snowpack model interface — a FESM-style static library (`libchion.a`)
with a clean public API and pluggable snowpack models (BESSI, PDD, ITM), ported
from Chion.jl.

## Docs

- [docs/PLAN.md](docs/PLAN.md) — implementation plan and architecture
- [docs/steady_state_snowpack.md](docs/steady_state_snowpack.md) — standalone
  domain-driven spin-up (Greenland GRL-16KM/8KM from raw `ice_data`, no
  preprocessing), MAR validation, performance, and the transmissivity study
- [docs/porting_notes.md](docs/porting_notes.md) — decisions made porting from Chion.jl
- [docs/pdd_defects.md](docs/pdd_defects.md) — PDD model defect notes

## Build

`make` builds `libchion.a` (OpenMP by default; `make openmp=0` for serial).
`make tests` builds the acceptance tests, `make grid` the gridded driver.
