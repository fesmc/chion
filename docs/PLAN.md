# chion — implementation plan

Fortran port of **Chion.jl v0.2.0** (`~/models/Chion.jl`, branch `main`), packaged as a
FESM-style static library `libchion.a` with a clean public API and pluggable snowpack models.

> **This plan supersedes an earlier draft** that was written against the stale `alex-dev`
> branch. Everything below is based on `main`.

## Reference material

| Source | Role |
|---|---|
| `~/models/Chion.jl` @ `main` (~8,400 LOC) | the model being ported — authoritative for all physics |
| `~/models/Chion.jl/docs/src/processes/*.md` | intended physics, prose form |
| `~/models/smbpal/src/smb_itm.f90` | the ITM scheme, which **does not exist in Chion.jl yet** |
| `~/models/yelmo` | Fortran conventions: `par`/`now` split, `select case` dispatch, nml, ncio, Makefile layout |
| `~/models/fesm-utils` | `precision`, `nml`, `ncio`, `variable_io`, `timestepping`, `timeout` |
| `~/models/configme` (`docs/DESIGN.md` §7) | build-configuration contract |
| `~/models/runme` (`src/runme/config.py`) | job-submission contract |
| `~/models/yelmox` (`libs/yelmox_domain.f90`) | eventual host; SMB wiring pattern |

---

## 1. What Chion.jl actually is

Facts that drive every decision below.

**Column-list discretisation, not a 2D grid.** State is `ncol` *independent* columns.
Layered arrays are `(Ntot, ncol)` — **layer-major**, layer 1 = surface, increasing downward
(`src/state.jl:147`, `src/column_state_utils.jl:41`). Julia and Fortran are both column-major,
so `real(wp) :: mass(Ntot,ncol)` is a byte-identical layout — no transposition anywhere.
`SnowpackGrid` carries `(x, y, js, is, mask)` to scatter columns back onto a 2D grid.

**Two models today, three eventually.**

| model | state | status |
|---|---|---|
| `BESSIModel` | `(Ntot,ncol)` layers + 15 per-column scalars | full physics |
| `PDDModel` | 4 per-column scalars, **no layers, no energy balance** | self-contained, but **not fully working** — see §3.2 |
| `ITMModel` | — | **does not exist in Chion.jl**; `build_model(:itm,…)` deliberately errors (`test/test_case_api.jl:538`). chion ports it from `smbpal/src/smb_itm.f90` instead (WP12). |

**Model choice and BESSI sub-options are different axes.** Within BESSI:
`albedo_scheme ∈ {constant, dynamic, prescribed}`, `fresh_snow_density_scheme ∈ {constant,
parameterized}`, `low_density_densification ∈ {bessi, htessel}` — stored as `UInt8` flags in
the constants struct (`src/constants.jl:43-49`).

**`SnowpackStepForcing` is already the backend contract** (`src/forcing.jl:234-257`): a
per-column, per-substep struct of 22 scalars with `has_*` flags marking which optional fluxes
are prescribed vs. internally parameterised. The Fortran port mirrors it 1:1.

**`column_step_core!` is the single-column kernel** (`src/step.jl:159-397`), device-safe and
called from a KernelAbstractions kernel over columns. Fixed order of operations:

```
1  accumulation (snow + rain, layer split/merge, depth cap)
2  if fresh snow onto a bare column: T(1) = T_air
3  no surface snow?  ->  bare-ice ablation, accumulate diagnostics, RETURN
4  albedo update (prescribed | dynamic | constant)
5  HTESSEL only: snapshot mass_w(1:n)
6  densification
7  latent-heat coefficients -> implicit energy solve (tridiagonal)
8  post-solve surface vapor mass flux
9  melt (if needs_melt): melt_mass = melt_energy_available / Lm
10 percolation -> runoff
11 HTESSEL only: liquid-water compaction
12 refreezing
13 final albedo fixup
```

**The only implicit solve is heat conduction** — backward Euler, Thomas algorithm, linearized
surface flux, zero-flux bottom, plus a two-pass melting-point re-solve. Everything else is
explicit forward Euler over `dt_seconds`.

**Coupling API already exists and is FastIsostasy-shaped**
(`src/integrators.jl`, `src/simulation.jl:96-113`):
`init_integrator(sim)` → `step!(integrator, Δt_days, force_dt=true)` → `finalize!(integrator)`,
with `sync_forcing!` to push externally-mutated forcing and `set_active_mask!` to
switch columns on/off at runtime. This maps directly onto `chion_init` / `chion_update` /
`chion_end`.

**No literature citations exist anywhere** in the process files or docs. Provenance is only
the informal names "BESSI", "HTESSEL", "PISM". Many bare magic numbers must be carried across
verbatim (see §5).

---

## 2. Design

### 2.1 Public API

`src/chion.f90` is the only module a host program needs.

```fortran
use chion

type(chion_class) :: chn

call chion_init(chn, filename, ncol, group="chion")        ! params + allocate; no physics
call chion_init_state(chn, time)                           ! cold start, or restart read
do while (...)
    ! host writes chn%forc%<field>(:) directly, then:
    call chion_update(chn, dt_days)
end do
call chion_end(chn)
```

Supporting public routines:

```fortran
chion_set_grid(chn, x, y, js, is, mask)      ! optional spatial mapping, for IO only
chion_set_active_mask(chn, active)           ! runtime column on/off (mirrors set_active_mask!)
chion_write_init(chn, filename, time)        chion_restart_write(chn, filename, time)
chion_write_step(chn, filename, time)        chion_restart_read (chn, filename, time)
```

Conventions taken from yelmo: object first, parameter file second, everything else
keyword/optional; `_init` loads params + allocates and computes no state; `_end` deallocates
only; the namelist group name is an argument so the package can be instantiated more than once.

**Deviation from yelmo, deliberate:** `chion_update` takes `dt_days` and advances exactly one
step, rather than taking a target `time` and sub-cycling. This matches Chion.jl's
`step!(integrator, Δt_days, force_dt)` and keeps the host in control of the time loop.
The diurnal sub-stepping inside `column_step` is the only internal sub-cycling.

### 2.2 Types

```fortran
type chion_class
    type(chion_param_class)   :: par      ! model choice, Ntot, mass thresholds, sub-scheme flags
    type(chion_const_class)   :: c        ! physical constants  (mirrors SnowpackPhysicalConstants)
    type(chion_grid_class)    :: grd      ! ncol, x, y, js, is, mask, active(:), active_idx(:)
    type(chion_forcing_class) :: forc     ! (ncol) arrays — the host writes these
    type(bessi_class)         :: bsi      ! allocated only when model=="bessi"
    type(pdd_class)           :: pdd
    type(itm_class)           :: itm
end type
```

`chion_step_forcing_class` — 22 scalars, one column, one substep, mirroring
`SnowpackStepForcing` field for field and in the same order:

```fortran
type chion_step_forcing_class
    real(wp) :: air_temperature, dt_days, snowfall_rate, rainfall_rate
    real(wp) :: shortwave_down, wind_speed
    real(wp) :: q_sw_net, q_lw_down, q_sh, q_lh
    logical  :: has_q_sw_net, has_q_lw_down, has_q_sh, has_q_lh
    real(wp) :: relative_humidity
    logical  :: has_relative_humidity
    real(wp) :: air_pressure, prescribed_albedo
    logical  :: has_prescribed_albedo
    real(wp) :: latitude_deg, day_of_year, solar_longitude_deg
end type
```

BESSI state, layer-major, matching `BESSIState` exactly:

```fortran
type bessi_state_class
    integer,  allocatable :: n_lay(:)          ! (ncol) active layer count
    real(wp), allocatable :: mass(:,:)         ! (Ntot,ncol) [kg m-2] solid
    real(wp), allocatable :: mass_w(:,:)       ! (Ntot,ncol) [kg m-2] liquid
    real(wp), allocatable :: density(:,:)      ! (Ntot,ncol) [kg m-3]
    real(wp), allocatable :: temperature(:,:)  ! (Ntot,ncol) [K]

    ! Cumulative accumulators -- MUST be wp_acc (dp). See section 3.1.
    real(wp_acc), allocatable :: mass_base(:), smb_ice(:), runoff(:)
    real(wp_acc), allocatable :: melt(:), refreezing(:)
    real(wp_acc), allocatable :: vapor_mass(:), sublimation(:), latent_heat_flux_sum(:)

    ! Instantaneous and diagnostic per-column scalars
    real(wp), allocatable :: t_srf(:), albedo(:)
    real(wp), allocatable :: thickness(:), wet_mass(:), bulk_density(:), liquid_water(:)
end type
```

PDD state is four `(ncol)` vectors: `snowpack_swe` (`wp`), and `smb_ice`, `runoff`, `pdd_sum`
(all `wp_acc`, being cumulative).

### 2.3 Model dispatch

`select case` on a character parameter, per yelmo — no abstract types, no procedure pointers.
The dispatcher `src/chion_model.f90` owns the column loop:

```fortran
select case(trim(chn%par%model))
    case("bessi")
        !$omp parallel do private(i,icol,forc_col)
        do i = 1, chn%grd%n_active
            icol = chn%grd%active_idx(i)
            call chion_pack_step_forcing(chn%forc, icol, dt_days, forc_col)
            call bessi_column_step(chn%bsi, icol, forc_col, chn%par, chn%c)
        end do
        !$omp end parallel do
    case("pdd")
        call pdd_step(chn%pdd, chn%forc, dt_days, chn%grd%active_idx, chn%par)
    case("itm")
        call itm_step(chn%itm, chn%forc, dt_days, chn%grd%active_idx, chn%par)
    case DEFAULT
        write(io_unit_err,*) "chion_update:: Error: model not recognized."
        write(io_unit_err,*) "model should be one of: ['bessi','pdd','itm']"
        write(io_unit_err,*) "model = ", trim(chn%par%model)
        stop
end select
```

The three sub-scheme flags dispatch the same way, as integer parameters inside `bessi_*`
routines (`CHION_ALBEDO_DYNAMIC` etc.), exactly mirroring the Julia `UInt8` flags.

**Energy workspace:** Julia allocates `(Ntot,ncol)` scratch (`EnergyWorkspace`). Since
`Ntot` is small (default 15), the Fortran port uses **automatic arrays of size `Ntot` local to
`bessi_energy_flux`**, which are stack-allocated per OpenMP thread and eliminate the workspace
type entirely. Note two Julia quirks that disappear with this: `layer_thickness` and
`thermal_conductivity` are allocated but never used, and `previous_temperature` is actually
reused as the diagonal scratch copy. Only 5 work arrays are genuinely needed.

### 2.4 Repository layout

```
chion/
  config/     Makefile (template), common.mk, Makefile_chion.mk
  .configme/  manifest.toml, .gitignore
  .runme/     info.json, config.default.toml
  src/
    chion.f90                 ! facade: re-export public API
    chion_defs.f90            ! precision, tolerances, constants, forcing/param/grid types
    chion_api.f90             ! chion_init / _init_state / _update / _end / _set_active_mask
    chion_model.f90           ! model dispatch + column loop + forcing pack/unpack
    chion_io.f90              ! ncio write / restart
    physics/
      snow_column_utils.f90   ! layer accessors, _surface_has_snow, bulk helpers
      snow_layers.f90         ! split / merge / shift / bottom-deplete / depth cap
      snow_accumulation.f90   ! snowfall+rain, fresh-snow density
      snow_albedo.f90         ! constant | dynamic | prescribed
      snow_densify.f90        ! BESSI | HTESSEL low-density + mid + high branches
      snow_energy.f90         ! tridiagonal implicit solve, Thomas, two-pass melt
      snow_surface_fluxes.f90 ! LW/SH/LH/rain heat, vapor mass, bare-ice ablation
      snow_melt.f90           ! melt application
      snow_percolation.f90    ! bucket cascade, irreducible water
      snow_refreezing.f90     ! cold content, refreeze, density update
      snow_diurnal.f90        ! solar geometry, interval averages, substep count
      snow_bessi.f90          ! assembles the above into bessi_column_step
      snow_pdd.f90            ! bulk PDD + PISM expectation integral
      snow_itm.f90            ! NEW — port from smbpal
      snow_diagnostics.f90    ! summarize_domain_state equivalents
  input/  chion_defaults.nml, chion_phys_const.nml, chion-variables-*.md
  par/    chion_Greenland.nml
  tests/  chion_column.f90, chion_grid.f90, test_layers.f90, test_energy.f90
  validation/  compare_julia.jl, forcing generators, reports
  libchion/{include,bin}/
  docs/   PLAN.md, porting_notes.md
```

---

## 3. Design decisions (settled)

1. **Discretisation: packed column list.** `ncol` independent columns with `(js,is)` scatter
   for IO, exactly as Chion.jl. Not `(nx,ny)`. The host passes only the points it wants
   computed.
2. **Precision: `wp = sp`, with `wp_acc = dp` for accumulators.** Matching yelmo and
   fesm-utils, so no conversion is needed at the yelmox boundary. Measured, not assumed —
   see §3.1.
3. **Cleanups allowed.** The port targets *physical* equivalence with Chion.jl, not
   bit-comparability. See §4.1 for the policy on what may be cleaned and what may not.
4. **ITM: port from `smbpal/src/smb_itm.f90` now** (WP12, promoted into batch 1). chion will
   ship with three models. Since chion is to replace smbpal in yelmox, chion-ITM must be able
   to reproduce smbpal-ITM — that is an explicit acceptance test.
5. **Parallelism: OpenMP over columns only.** No GPU path.
6. **chion replaces smbpal in yelmox**, eventually. WP19 is therefore a migration, not an
   addition: it must keep smbpal available during a transition period and demonstrate
   equivalence before smbpal is removed.

### 3.1 Precision policy

Chion.jl is `Float64` throughout. chion is `sp` for state and interfaces, which is a real
change and needed evidence rather than a guess. The expressions below were extracted from the
Chion.jl algorithms and evaluated in both precisions.

| expression | `sp` behaviour | verdict |
|---|---|---|
| cumulative diagnostics, 36,500 daily increments | 0.01 kg m⁻² increments onto a 1e5 kg m⁻² total are lost outright; **80 kg m⁻² drift (0.08%) over 100 yr** | **`dp` required** |
| pore volume `phi = m/rho - m/rho_i` | abs error ~3e-8 m; quantized at ~1e-5, so the `TOL_TINY` guard **can never fire**, and `phi` feeds the division `lwc = m_w/rho_w/phi` | **`dp` locally** |
| cold content `(T0-T)*ci*m_s` | 100% relative error below 3e-5 K from melting, but implied mass error only ~1e-5 kg m⁻² | `dp` locally (free) |
| densification gap `(rho_i - rho)` | error ~2e-5 kg m⁻³, absorbed by the monotone update and `rho_i` cap | `sp` fine |
| tridiagonal conduction solve | max 3.4e-5 K vs `dp`; matrix is strictly diagonally dominant | `sp` fine |
| melt energy residual | no error at representative values | `sp` fine |

Consequently:

- `wp = sp` — layered state (`mass`, `mass_w`, `density`, `temperature`), forcing, and all
  public interfaces.
- `wp_acc = dp` — the nine cumulative per-column accumulators: `smb_ice`, `runoff`, `melt`,
  `refreezing`, `vapor_mass`, `sublimation`, `latent_heat_flux_sum`, `mass_base`, `pdd_sum`.
  Per-column scalars, so the memory cost is negligible.
- `real(wp_acc)` **locals** in three named places: pore volume (percolation, albedo), cold
  content (refreezing), and surface-energy accumulation across diurnal substeps.
- Both tolerances are declared `dp` so comparisons promote correctly. A guard of the form
  `x <= TOL_TINY` is only meaningful when `x` itself was computed in `dp`.

Note the memory argument for `sp` is negligible at these sizes (Greenland 10 km ≈ 20,000
columns × 15 layers × 4 arrays = 4.8 MB vs 9.6 MB). The actual benefit is interface
consistency with yelmo/yelmox.

**This changes every conservation tolerance in the plan.** `sp` gives ~7 significant digits,
so acceptance thresholds are `1e-6` relative, not `1e-12`. Where a quantity is `wp_acc`, the
tighter `dp` threshold still applies and is stated per WP.

### 3.1b Deferred decisions — follow up before WP19

Open questions parked deliberately, each with the evidence already gathered. None blocks the
build; all should be resolved before the yelmox cutover.

| # | Question | Evidence so far | Decision |
|---|---|---|---|
| 1 | `8.13` vs the gas constant `8.314` in the densification Arrhenius denominators | The WP7 factors (~2.6, ~1400) were **arithmetically wrong**; the true factors are **1.11** and **1.86**. Provenance found: Herron & Langway (1980) stage 1 is `k0 = 11*exp(-10160/(R*T))` with `R = 8.314`, and the `10160` matches exactly | **RESOLVED (D22).** Corrected to the gas constant. 10-yr column run: <1% on every integrated quantity. Chion.jl issue #18. |
| 2 | `9.81` vs `9.80665` in the overburden | 3.6e-4 relative on overburden, ~1.1e-3 on the cubed mid/high tendencies (WP7) | **RESOLVED (D25).** Unified onto standard gravity, which is exact by definition. Covered by `legacy_chion=1`. |
| 3 | ITM's `L_m = 3.35e5` vs `chion_const_class%Lm = 3.34e5`, and its hard-coded `273.15` | `itm_c`/`itm_t` are calibrated against them; the true latent heat of fusion is 3.337e5, so 3.34e5 is the more accurate | **RESOLVED (D26).** ITM reads all four constants from `chion_const_class`. +0.30% on potential melt, measured. Two constant sets would let a host retuning `T0`/`Lm` leave ITM on different physics from BESSI. |
| 4 | PDD `smb_ice` convention: whole-column mass change (PDD) vs ice-only forcing (BESSI) | the two Chion.jl models disagree; `smb_ice` feeds `ice_sheet_net_forcing_yearly` (WP9). PDD has no firn representation, so a whole-column `smb_ice` implies a reservoir it does not have | **RESOLVED (D23).** BESSI's ice-facing convention adopted, with a capped one-layer reservoir as in smbpal. Full mass closure now holds. Chion.jl issue #19. |
| 5 | Whether `pdd_method` should default to `pism` rather than `simple` | the simple form loses 1.25 kg m-2 d-1 at -5 C and 5.98 at 0 C, concentrated at the ELA; smbpal always uses the integral (WP9) | Recommend `pism` at WP13. |

### 3.2 Caution: PDD is not fully working in Chion.jl

Chion.jl's `PDDModel` is known not to be fully functional. Consequences for the port:

- **Chion.jl is not authoritative for PDD.** Port the structure, but cross-check the physics
  against `smbpal/src/smb_pdd.f90` (`calc_ablation_pdd`, `calc_temp_effective`), which is the
  scheme actually in production use in yelmox.
- **WP9 must produce a defect list**, not just a port. Anything that looks wrong in
  `Chion.jl/src/processes/pdd.jl` gets written up in `docs/porting_notes.md` and reported back
  so it can be fixed upstream in Chion.jl.
- Known suspicious points already identified, to be resolved before implementing:
  - `refreezing_fraction` is applied to snow melt with **no cold-content or capacity limit** —
    a flat 0.6 of all snow melt refreezes regardless of temperature. smbpal's
    `calc_ablation_pdd` uses an `f_refrz_max` capacity limit instead.
  - Refrozen mass is simultaneously added back to `snowpack_swe`, credited to `smb_ice`, and
    subtracted from `runoff` — check this is not double-counting.
  - `snowpack_swe` has no cap, no aging and no densification, so it can grow without bound.
  - The "monthly" branch is auto-selected by `27 <= dt_days <= 32`, an implicit and fragile
    trigger for a change of physics.
  - `pdd_step!(model, state, forcing)` at `pdd.jl:366` applies monthly detection per step,
    while the loop at `pdd.jl:283` bypasses it entirely — two entry points with different
    physics.
- **PDD validation (WP16) therefore has two reference targets**: Chion.jl for structural
  agreement, and smbpal for physical plausibility. Where they disagree, smbpal wins and the
  divergence is documented.

---

## 4. Work packages

Each WP is self-contained, names its source-of-truth Julia file, its deliverable files, and an
acceptance test. Dependencies listed; same-level WPs run in parallel.

### 4.1 Cleanup policy

The port targets physical equivalence, not bit-comparability. To keep "cleanup" from turning
into "silent behaviour change", every WP applies this policy and records each cleanup in
`docs/porting_notes.md` with a one-line justification.

**Allowed without asking:**
- Replace numerical approximations with exact library functions where the difference is at
  round-off level — notably the Abramowitz–Stegun `_normal_cdf` polynomial (~1e-7 error) with
  `0.5_dp*erfc(-x/sqrt(2.0_dp))`.
- Promote magic numbers to named parameters, or to namelist parameters with the Julia value as
  the default. §5 item 3 lists the candidates. Prefer namelist parameters for anything that is
  arguably tunable physics (`max_lwc`, albedo aging coefficients, densification thresholds);
  named `parameter` constants for genuine physical constants.
- Reconcile `g = 9.81` in densification with the single gravity constant — **but** measure the
  impact on a decade-long column run first and record it. Measured in WP7: 3.6e-4 relative on
  overburden, entering the mid/high branches cubed for ~1.1e-3 on those tendencies. Small, but
  still untested end-to-end until WP8 exists.
- Delete the two unused `EnergyWorkspace` arrays and stop misusing `previous_temperature` as
  diagonal scratch (§5 item 12). The workspace becomes stack-local automatic arrays anyway.
- Remove the duplicate/inconsistent entry points that exist only for GPU dispatch in Julia.
- Fix outright bugs, provided they are reported back for upstream fixing rather than silently
  diverging.

**Not allowed without asking first** — these look like cleanups but change behaviour:
- Unifying the three empty-layer thresholds (§5 item 1). They gate different physics.
- Reconciling the linearized vs. exact surface-flux evaluations (§5 item 2).
- Adding `dt` to the albedo aging law (§5 item 5), even though it is dimensionally wrong.
- Making the two-pass energy re-solve a true Dirichlet row (§5 item 6).
- Enforcing volume conservation in the refreezing density cap (§5 item 7).
- Making the depth cap respect the configured `Ntot` (§5 item 11).
- **Replacing `8.13` with the gas constant `8.314` in the densification Arrhenius
  denominators.** This was originally listed as an allowed cleanup pending measurement. WP7
  measured it and it is **not** a cleanup: `exp(-10160/(8.13*T)) / exp(-10160/(8.314*T))` at
  260 K is a factor of **~2.6**, so substituting the gas constant more than doubles the
  low-density BESSI densification rate. The mid/high branches carry `-60000/(8.13*T)`, where
  the same substitution is a factor of **~1400** at 263 K. A number that load-bearing is a
  calibrated parameter, not a unit slip. It is ported verbatim as the named parameter
  `DENSIFY_R_GAS` so the experiment is a one-line change.
- **Reconciling ITM's `L_m = 3.35e5` with `chion_const_class%Lm = 3.34e5`,** or its hard-coded
  `273.15` with `T0`. `itm_c` and `itm_t` are calibrated against those values, so a host
  retuning `T0` would silently retune ITM's melt.

Each of these is a real modelling question, not a style issue. Raise them; do not decide them
inside a WP.

### Level 0 — foundations (sequential, must be first)

**WP1 — `chion_defs.f90` + repo skeleton**
Deliverables: the directory tree in §2.4 with placeholders, `.gitignore`, and
`src/chion_defs.f90` containing: `sp`/`dp`/`wp` re-exported from fesm-utils `precision`;
`TOL_TINY = 1.0e-12_wp_acc`, `TOL_EMPTY_LAYER = 1.0e-10_wp_acc` (both `dp`, see §3.1),
`io_unit_err`, `MV`;
`chion_const_class` (26 fields mirroring `SnowpackPhysicalConstants`, `src/constants.jl:61-88`,
defaults at `:173-200`); the scheme-flag integer parameters; `chion_step_forcing_class`;
`chion_forcing_class` (the `(ncol)` host-facing arrays); `chion_param_class`;
`chion_grid_class`; `chion_class`; and the yelmo-style helpers `chion_check_enum`,
`chion_check_file`, `chion_parse_path`, `chion_load_command_line_args`.
Source of truth: `Chion.jl/src/constants.jl`, `src/forcing.jl:234-257`, `src/domain.jl`.
Reference for style: `yelmo/src/yelmo_defs.f90`. **No physics.**
Acceptance: compiles standalone; a trivial program prints all default constants.
Blocks: everything.

**WP2 — build system (configme-compatible)**
Deliverables: `config/Makefile` template declaring `srcdir/objdir/bindir/libdir`,
`debug ?= 0`, `openmp ?= 0` **before** the single `<COMPILER_CONFIGURATION>` token, followed
immediately by `include config/common.mk`, the `DFLAGS` switch, and
`include config/Makefile_chion.mk`; `config/common.mk` wiring `FESMUTILSROOT`/`INC_FESMUTILS`/
`LIB_FESMUTILS`, `LIB_NC`, the `openmp=1` `include-serial`→`include-omp` swap,
`LFLAGS_EXTRA ?= -Wl,-zmuldefs`, and the final `LFLAGS`; `config/Makefile_chion.mk` with one
explicit rule per object plus object lists `chion_physics`, `chion_base`, `chion_tests`;
targets `chion-static` (→ `libchion/include/libchion.a`, `ranlib`,
`./check_githash > git_chion.txt`), `column`, `grid`, `clean`, `usage`;
`.configme/manifest.toml` (`package = "chion"`, `deps = ["fesm-utils:dev"]`) +
`.configme/.gitignore`; `check_githash` copied from yelmox; committed empty
`libchion/include/` and `libchion/bin/`.
Reference: `yelmo/config/*`, `yelmox/config/*`, `configme/docs/DESIGN.md` §7.
Acceptance: `configme -m macbook -c gfortran` generates a `Makefile`; `make chion-static`
builds the WP1 stub into `libchion.a`.
Depends: WP1.

**WP3 — `snow_column_utils.f90`**
Port `Chion.jl/src/column_state_utils.jl`. In Fortran the accessors collapse to plain array
indexing, so the real deliverables are the *predicates* and their exact thresholds:
`surface_has_snow` (`n>0 .and. mass(1,i) > EPS_EMPTY_LAYER`, from
`src/processes/energy_flux.jl:303`), `column_has_liquid_water` (`> EPS_TINY`, not
`EPS_EMPTY_LAYER` — the docstring is wrong, follow the code), `bulk_snow_density`,
`total_snow_water_mass` (note: clips negatives, unlike the diagnostics version).
**Three different "empty" thresholds are used deliberately** — `> 0`, `> EPS_TINY`,
`> EPS_EMPTY_LAYER`. Do not unify them.
Acceptance: unit test asserting each predicate flips at exactly the documented threshold.
Depends: WP1.

### Level 1 — physics kernels (4 agents in parallel)

All of these are pure single-column routines over plain arrays. None of them touch derived
types beyond `chion_const_class` and `chion_param_class`.

**WP4 — layer dynamics (`snow_layers.f90`)**
Port `Chion.jl/src/processes/layer_structure.jl` in full:
`reset_layer_at_index`, `split_surface_layer`, `merge_surface_layer` (both the partial-transfer
and full-merge branches), `merge_bottom_layer` (with the `combined_density > rho_i` export to
`mass_base`/`smb_ice`), `remove_surface_layer`, `remove_depleted_surface_and_route_water`,
`continuous_bottom_deplete`, `free_slot_for_surface_split` (including the `Ntot == 1` special
case), `enforce_snow_depth_cap`.
Note the depth cap uses hard-coded `BESSI_REFERENCE_LAYER_COUNT = 15` and
`BESSI_REFERENCE_DEPTH_DENSITY = 300` with a `1.5` factor, independent of the configured `Ntot`.
Deliverable also: `tests/test_layers.f90` driving accumulation/melt sequences and asserting
`sum(mass) + sum(mass_w) + mass_base + runoff` is conserved to `1e-6` relative (sp; see
§3.1). Note `mass_base` and `runoff` are `wp_acc`, so accumulate the check in `dp`.
This is the highest-risk WP — the split/merge/shift index arithmetic is fiddly.
Depends: WP1, WP3.

**WP5 — energy (`snow_energy.f90`, `snow_surface_fluxes.f90`)**
Port `Chion.jl/src/processes/energy_flux.jl` and `surface_fluxes.jl`. Contents:
- Thomas algorithm (forward/backward), with Julia's indexing convention: `lower(k)` is the
  sub-diagonal entry of **row k+1**. No pivoting.
- Matrix assembly: `beta_i = -2*dt/(ci*safe_pos(mass(i)))`;
  `interface(k) = (K_k*dz_k + K_{k+1}*dz_{k+1}) / safe_pos((dz_k+dz_{k+1})**2)`;
  `K = Ki*(rho*1e-3)**1.88`. Surface row carries `+lambda*Q_lin` on the diagonal and
  `+lambda*Q_const` on the rhs, `lambda = dt/(ci*m1)`. Bottom row is zero-flux.
- The **single-layer shortcut** (`n == 1`) is a separate closed-form branch, not the solver.
- The **two-pass melting-point re-solve**: if `T(1) > T0` after pass 1, rebuild the rhs from
  the original temperatures, set `rhs(1) = T0`, subtract `diag_surf` from `diag(1)`, re-solve,
  add `(T0 - T(1))*ci*m1` to `energy_to_melting`, then force `T(1) = T0` and clamp all layers.
  This is **not** a strict Dirichlet row — get this exactly right.
- `melt_energy_available = max((Q_const - Q_lin*T_srf)*dt - energy_to_melting, 0)`.
- Return a derived type with the 8 documented fields.
- `surface_fluxes.f90`: the **unlinearized** twin used for bare ice (evaluated at `T = T0`) and
  for post-solve vapor mass (evaluated at the new `T`). The two evaluations are deliberately
  inconsistent with the linearized ones in the solver — reproduce, do not "fix".
- Vapor-pressure parameterisations: `e_sat_w` over water for air (`611.2, 17.27, 243.12`),
  `e_sat_i` over ice for the surface (`611.2, 22.46, 272.62`), and its derivative;
  `D_lf = latent_heat_flux_ratio * D_sh/cp_air * 0.622 * (Lv+Lm)`.
- `bare_ice_ablation_mass` returning `(melt_mass, vapor_mass, sublimation_mass,
  latent_heat_flux, net_mass_change)`; sign convention **positive vapor_mass = deposition**.
- `diagnose_latent_heat_flux_coefficients` returning `(linear, constant)` — snowfall takes
  precedence over rainfall.
Note `shortwave_absorbed` has **no** `max(SWdn,0)` in `energy_flux.jl:99` but **does** in
`surface_fluxes.jl:74`. Preserve both.
Acceptance: `tests/test_energy.f90` — a steady-state analytic conduction case
(cf. `Chion.jl/docs/src/tests/test_energy_flux_analytical.md`), plus a check that the
tridiagonal solver reproduces a dense LU solve to `1e-14` on random SPD systems.
Depends: WP1, WP3.

**WP6 — water (`snow_percolation.f90`, `snow_refreezing.f90`, `snow_melt.f90`)**
Port `Chion.jl/src/processes/percolation.jl`, `refreezing.jl`, `melt.jl`.
- Percolation: single-pass top-down bucket cascade. Pore volume
  `phi = m_s/rho - m_s/rho_i`; retention `0.1 * rho_w * phi` (fixed volumetric fraction —
  **not** a density-dependent Coléou–Lesaffre form). Three cases per layer: no solid mass →
  water goes **straight to runoff** (note: the doc page says otherwise; follow the code);
  collapsed pore space (`phi <= EPS_TINY`) → all water pushed down or to runoff; normal
  retention → excess pushed down or to runoff. "Lowest layer" test is
  `k == n .or. mass(k+1,i) <= 0`, short-circuited so `mass(k+1,i)` is never evaluated when
  `k == n`. Percolation does **not** modify `mass`, `density` or `temperature`.
- Refreezing: per-layer, independent, no vertical coupling. Entry condition is the strict
  triple `m_s > 0 .and. m_w > 0 .and. T < T0` (strict zeros, no `EPS_TINY`).
  Cold content `(T0-T)*ci*m_s` vs. available latent `m_w*Lm`. Partial branch sets `T = T0`
  exactly; complete branch uses the mixing formula `(m_w*Lm/ci + m_w*T0 + T*m_s)/(m_w+m_s)`.
  Density scales by the mass ratio and is capped at `rho_i` **while mass still gains the full
  increment** — mass is conserved, volume is not. Preserve.
- Melt: fast path when `melt < mass(1) - EPS_EMPTY_LAYER`; otherwise a loop consuming layers
  from the top, with depleted-layer removal and water routing, then a surface merge when
  `mass(1) < mass_min`.
Acceptance: unit tests for (i) water conservation through percolation, (ii) refreezing energy
closure `refrozen*Lm == Σ consumed cold content` to `1e-6` relative, (iii) melt mass
conservation to `1e-6` relative. Cold content is computed in `dp` locally (§3.1), so the
energy closure may be checked more tightly than the `sp` mass fields.
Depends: WP1, WP3, WP4 (melt calls the layer routines).

**WP7 — accumulation, albedo, densification, diurnal**
Port `Chion.jl/src/processes/accumulation.jl`, `albedo.jl`, `densification.jl`,
`diurnal_shortwave.jl`.
- Accumulation: fresh-snow density (constant `rho_s`, or the parameterized
  `rho_s_a + rho_s_b*(T-T0) + rho_s_c*sqrt(max(wind,0))`, both clamped to `[50, rho_i]`);
  volume-weighted density mixing on the surface layer; the split `while` loop (with the
  `Ntot <= 2` special case calling `free_slot_for_surface_split` instead of
  `merge_bottom_layer`); the merge `while` loop; then the depth cap.
- Albedo, three schemes. Dynamic law, applied **once per call, with no `dt`** (so it is
  timestep-dependent — call exactly once per step, as Julia does):
  `a = min(a_prev, a_prev - (1.35e-3*(Ts-T0) + 0.0278))`, floored at `alpha_wet`; then the
  wetness relaxation `a - (a - alpha_wet)*lwc/max_lwc_albedo`; then clamp to
  `[alpha_wet, alpha_dry]`. Snowfall brightening:
  `min(alpha_dry, a + (alpha_dry-alpha_wet)*(1 - exp(-snowfall_mass/3)))` — e-folding **3 kg m-2**,
  and a single event can brighten by at most `alpha_dry - alpha_wet`.
  `ALBEDO_PRESCRIBED` takes the *dynamic* code path inside these routines and is overridden by
  the caller — reproduce that structure.
- Densification: overburden `sigma = 9.81*(M_above + m/2)` (note **9.81**, not 9.80665);
  three density regimes `<550`, `[550,800)`, `>=800`; the scheme flag selects **only** the
  `<550` branch (BESSI: `0.011*exp(-10160/(8.13*T))*(rho_i-rho)*max(A_t,0)` — note **8.13**,
  not the gas constant; HTESSEL: `rho*(sigma/eta + xi)` with the viscosity and thermal-
  metamorphism forms). Hard-coded `rho_e = 815`, `P_atm = 101325`. Monotone update
  `max(rho, rho + drho*dt)` capped at `rho_i`. Plus `apply_htessel_liquid_water_compaction`.
- Diurnal: declination from `solar_longitude_deg` with obliquity `23.439291`; sunset hour
  angle with polar day/night branches; `I_day = 2*(h0*A + B*sin(h0))`; interval average
  `S*I_ab/w` with `S = Qbar*2*pi/I_day`, divided by the **full** interval width so nocturnal
  intervals dilute correctly; temperature interval average
  `Tbar + amplitude*(sin(h_b)-sin(h_a))/w`. `substep_count` returns **only 1 or
  `max_substeps`**, gated by eight conditions including `dt_days ∈ [0.75, 1.25]`.
Acceptance: unit tests for (i) fresh-snow density clamps, (ii) albedo bounds under extreme
inputs, (iii) densification monotonicity and the `rho_i` cap, (iv) diurnal interval averages
summing back to the daily mean over a full `[-pi,pi]` tiling.
Depends: WP1, WP3, WP4.

### Level 2 — model assembly

**WP8 — BESSI (`snow_bessi.f90`)**
Assemble WP4–WP7 into `bessi_par_class` / `bessi_state_class` / `bessi_class`,
`bessi_par_load`, `bessi_alloc`/`_dealloc`, `bessi_init_state`, `bessi_reset_columns`, and
`bessi_column_step` reproducing `column_step_core!` (`src/step.jl:159-397`) **in the exact
order listed in §1**, plus the diurnal substep wrapper `column_step!` (`src/step.jl:106-157`).
Watch: the early-return bare-ice branch skips percolation and refreezing entirely; the
`temperature(1) = T_air` fresh-snow-on-bare-column rule; the HTESSEL snapshot taken *before*
the energy solve and used *after* percolation.
Acceptance: a single column driven by a synthetic annual cycle produces a mass-conserving,
physically plausible snowpack; profiles printed for eyeball comparison with Julia.
Depends: WP4, WP5, WP6, WP7.

**WP9 — PDD (`snow_pdd.f90`)** — *read §3.2 first; Chion.jl's PDD is not fully working*
Port `Chion.jl/src/processes/pdd.jl`, cross-checking every step against
`smbpal/src/smb_pdd.f90` (`calc_ablation_pdd`, `calc_temp_effective`). Where they disagree,
smbpal wins and the divergence is written up. **Deliverable includes a defect list** for
upstream Chion.jl, covering at minimum the five suspicious points in §3.2.
Bulk, no layers. Six-line core:
```
pdd_sum        += pdd
available_snow  = snowpack_swe + snowfall
snow_melt       = min(available_snow, ddf_snow*pdd)
remaining_pdd   = max(pdd - snow_melt/ddf_snow, 0)
ice_melt        = ddf_ice*remaining_pdd
refrozen        = refreezing_fraction*snow_melt     ! snow melt ONLY; rain and ice melt never refreeze
snowpack_swe    = available_snow - snow_melt + refrozen
smb_ice        += snowfall - snow_melt + refrozen - ice_melt
runoff         += rainfall + snow_melt - refrozen + ice_melt
```
Two PDD flavours, auto-selected by timestep: simple `max(T-273.15,0)*dt` normally, and the
PISM/Calov–Greve expectation integral
`dt*(sigma*phi(z) + Tbar*Phi(z))`, `z = Tbar/sigma`, when `27 <= dt_days <= 32`.
Per §4.1, use `0.5_dp*erfc(-z/sqrt(2.0_dp))` rather than the Abramowitz–Stegun polynomial.
Replace the implicit `27 <= dt_days <= 32` trigger with an explicit
`pdd_method = "simple"|"pism"` namelist parameter, defaulting to the Julia behaviour; note
this in the defect list. Use the constants struct rather than PDD's hard-coded `273.15`
and `86400.0`.
Acceptance: mass-balance identity
`Δsmb_ice + Δrunoff == snowfall + rainfall - Δsnowpack_swe` to `1e-10` (all three terms are
`wp_acc`, so this stays a `dp` check); PISM integral checked
against its limits (`Tbar→+inf → dt*Tbar`; `Tbar=0 → dt*sigma/sqrt(2pi)`); and a run against
the same forcing as smbpal's PDD giving physically comparable annual SMB.
Depends: WP1.

**WP10 — diagnostics (`snow_diagnostics.f90`)**
Port `Chion.jl/src/diagnostics.jl`: `summarize_domain_state` computing `thickness`,
`wet_mass`, `bulk_density`, `liquid_water` with the exact (and mutually inconsistent)
inclusion criteria — thickness sums only layers with `mass > 0 .and. density > EPS_TINY`,
while `bulk_density`'s numerator is the **unconditional** mass sum, and `wet_mass` /
`liquid_water` do **not** clip negatives. Plus `get_state` / `print_state` for a single column.
Acceptance: matches a hand-computed 3-layer example including the edge cases.
Depends: WP1, WP3.

### Level 3 — integration

**WP11 — dispatcher + public API (`chion_model.f90`, `chion_api.f90`, `chion.f90`)**
The `(ncol)`→per-column forcing pack, the OpenMP column loop, the model `select case`, the
active-mask machinery (`set_active_mask` + column reset on deactivation, mirroring
`_reset_model_columns!`), `chion_init`/`_init_state`/`_update`/`_end`, and the `chion.f90`
re-export facade.
Acceptance: a 100-column run gives bit-identical results to 100 single-column runs, and to the
same run with `omp_num_threads=1` vs `=8`.
Depends: WP8, WP9, WP10.

**WP12 — ITM (`snow_itm.f90`)** — *moved into batch 1; source is smbpal, not Chion.jl*
Port `smbpal/src/smb_itm.f90` to the `chion_step_forcing_class` contract: `itm_par_class`,
`itm_par_load`, `itm_alloc`, `itm_init_state`, `itm_step`. Keep the physics unchanged
(`calc_itm`, `calc_atmos_transmissivity`, `calc_albedo_planet`, `itm_c_lat`,
`calc_albedo_surface`). Insolation arrives as forcing (`forc%q_sw_net` or `shortwave_down`);
chion does not own an insolation module — this is a change from smbpal, which computes
insolation internally via `calc_insol_day`, so the host must now supply it.
Because chion is to replace smbpal in yelmox, **acceptance is equivalence with smbpal**: given
the same parameters, forcing and insolation, `chion` with `model="itm"` must reproduce
`smbpal_update_monthly` with `abl_method="itm"` to within `sp` round-off (~1e-6 relative).
Build that comparison as
part of this WP, not as an afterthought — it is the evidence that WP19 is safe.
Also record, for a future upstream contribution, what an `ITMModel` in Chion.jl would need.
Depends: WP1, WP3.

**WP13 — parameters (`input/chion_defaults.nml`, `par/chion_Greenland.nml`)**
Full defaults schema in the yelmo column layout with `! [unit] description` comments and
enumerated allowed values; one `nml_read` per parameter in each `*_par_load`; `nml_validate`
against the defaults file; `chion_check_enum` for `model`, `albedo_scheme`,
`fresh_snow_density_scheme`, `low_density_densification`. Also `input/chion_phys_const.nml`
with an `&Earth` group carrying the 26 `chion_const_class` values.
Acceptance: every parameter read by the code appears in the defaults file, verified by a script.
Depends: WP8, WP9, WP11.

**WP14 — IO (`chion_io.f90`, `input/chion-variables-*.md`)**
`chion_write_init`/`chion_write_step` via `ncio` and `variable_io` markdown metadata tables
(one per model); `chion_restart_write`/`_read` covering the full prognostic state including
the `(Ntot,ncol)` BESSI arrays and `n_lay`. Output uses `scatter_to_grid` semantics
(`Chion.jl/src/io.jl:45`) when `(js,is,mask)` are set, otherwise writes a bare column
dimension. Variable names and units must match Chion.jl's `NETCDF_METADATA`
(`Chion.jl/src/io.jl:3-32`) so the two models' output files are directly comparable.
Depends: WP11, WP13.

**WP15 — driver programs (`tests/chion_column.f90`, `tests/chion_grid.f90`)**
`chion_column.x`: `ncol=1`, `&ctrl` group specifying synthetic forcing — the direct comparison
target for Chion.jl. `chion_grid.x`: driven by a forcing NetCDF file, mirroring
`load_forcing_file`. Both use `chion_load_command_line_args` (exactly one argv) and write to CWD.
Depends: WP11, WP13, WP14.

### Level 4 — validation and tooling

**WP16 — validation harness (`validation/`)** — *the most important WP after WP4/WP5*
A script that (i) generates a forcing NetCDF, (ii) runs `Chion.jl` on it, (iii) runs
`chion_grid.x` on it, (iv) reports per-field max-abs and max-rel differences per timestep, and
(v) fails on a threshold. Must cover, at minimum: a cold dry column (layers + densification
only), a melting column (energy + melt + percolation + refreezing), a column that goes bare and
recovers, a column at `Ntot` capacity exercising bottom merge and the depth cap, and a
monthly-timestep PDD run exercising the PISM branch.

Because cleanups are allowed (§4.1), the target is **physical equivalence, not
bit-comparability**. Set per-field relative tolerances explicitly and justify each one; a
tolerance that has to be loosened is a signal to re-read the code, not to loosen it further.
Every field that cannot meet a tight tolerance must be traced to a specific documented cleanup.

Three reference targets, not one:
- **BESSI** — Chion.jl is authoritative. Tight tolerances.
- **PDD** — Chion.jl for structure only (§3.2); smbpal for physical plausibility. Loose
  tolerances, plus a written explanation of every divergence.
- **ITM** — smbpal is authoritative (there is no Chion.jl ITM). Round-off agreement expected;
  this is the WP12 comparison, run as part of the same harness.
Depends: WP12, WP15.

**WP17 — runme integration**
`.runme/info.json` with all 8 required keys (`exe_default = "column"`,
`exe_aliases = {"column":"libchion/bin/chion_column.x","grid":"libchion/bin/chion_grid.x"}`,
`par_path_as_argument: true`, `links: ["input"]`), `.runme/config.default.toml` with all 7
keys (`jobname = "chion"`), `.gitignore` entry for `.runme/config.toml`, `output/` convention.
Depends: WP15.

**WP18 — configme registry entry**
`configme/data/packages/chion.toml` in the configme repo: `[package]` name/org/repo/dir,
`config_style = "makefile-template"`, `[[package.links]] dep = "fesm-utils"`,
`[package.build] make_target = "chion-static"`, `variants = ["serial","omp"]`,
`[package.artifacts]` listing `libchion/include/libchion.a`.
Depends: WP2.

**WP19 — yelmox integration: migrate off smbpal** — *do last*
chion replaces smbpal, but not in one step. Two phases:

*Phase A — coexistence.* `CHIONROOT`/`INC_CHION`/`LIB_CHION` in `yelmox/config/common.mk`, a
`chion-static` recursion target, `LIB_CHION` on `LFLAGS`; a `"chion"` case in `ctl%smb_method`
(`yelmox/libs/yelmox_domain.f90:81`) **alongside** the existing `"smbpal"`, including the
column pack/unpack between yelmox's `(nx,ny)` fields and chion's column list (build it from
the ice mask so only ice-covered points are passed); `&chion`, `&bessi`, `&pdd`, `&itm` groups
appended to `yelmox/yelmox/yelmox_Greenland.nml`. Both paths buildable and runnable.

*Phase B — cutover, only after evidence.* Run a full Greenland simulation with
`smb_method="smbpal"` and with `smb_method="chion"` (`model="itm"`, matched parameters) and
compare SMB fields, integrated mass balance, and the resulting ice-sheet evolution. Only once
that is signed off: switch the default, deprecate `"smbpal"`, and remove
`yelmox/libs/smbpal/` in a separate commit.

Note chion does not compute insolation (WP12), so yelmox must supply it — decide whether that
comes from the vendored `libs/insol/` or from the climate forcing.
Depends: WP16, WP17.

**WP20 — docs (rolling)**
`README.md`, `docs/api.md`, and `docs/porting_notes.md` recording every intentional deviation
from Chion.jl and every Julia quirk deliberately preserved (see §5).

---

## 5. Traps carried over from Chion.jl

Collected here because they will silently break a "reasonable" port. `docs/porting_notes.md`
must restate each one next to the code that honours it.

1. **Three distinct empty-layer thresholds**, used deliberately and not interchangeably:
   `> 0` (percolation entry, `is_lowest_active_snow_layer`), `> EPS_TINY = 1e-12`
   (`column_has_liquid_water`, surface LWC, `bulk_snow_density`, `column_summary`),
   `> EPS_EMPTY_LAYER = 1e-10` (`surface_has_snow`, constant albedo, albedo update).
   The energy-solve early exit uses yet another: `mass(1) <= 0`.
2. **Two inconsistent surface-flux evaluations by design**: linearized about `T^n` in the
   energy solve; exact at a known `T` in `surface_fluxes.f90` (at `T0` for bare ice, at
   `T^{n+1}` for post-solve vapor mass). The energy and mass budgets therefore use slightly
   different latent fluxes. Do not reconcile them.
3. **Magic constants not in the constants struct**, all of which must be replicated:
   `g = 9.81` in densification (vs `DEFAULT_GRAVITY = 9.80665` elsewhere); `8.13` in both
   densification Arrhenius denominators (not the gas constant `8.314`); `rho_e = 815`;
   `P_atm = 101325`; regime thresholds `550`/`800`; vapor-pressure coefficients
   `611.2/17.27/243.12` (water) and `611.2/22.46/272.62` (ice); `0.622`; obliquity
   `23.439291`; `dt_days ∈ [0.75, 1.25]`; albedo `1.35e-3`, `0.0278`, snowfall e-folding `3`;
   percolation `max_lwc = 0.1`; PDD's hard-coded `273.15` and `86400.0`.
4. **`max_lwc` (percolation, 0.1) and `c%max_lwc_albedo` (albedo, 0.1) are different
   parameters** that happen to share a default. Keep them independently configurable.
5. **The albedo aging law has no `dt`** — it decays per *call*. Call it exactly once per step.
6. **The two-pass energy re-solve is not a Dirichlet row.** Row 1 keeps its conduction
   coupling; only the surface flux term is removed and the rhs replaced by `T0`.
7. **The refreezing density cap breaks volume conservation on purpose** — mass gains the full
   increment while density is capped at `rho_i`.
8. **`_uses_htessel_densification` is an `if/else`, not a 3-way select** — any flag value other
   than HTESSEL falls through to BESSI.
9. **`ALBEDO_PRESCRIBED` takes the dynamic code path** inside the albedo routines and is
   overridden by the caller; if `has_prescribed_albedo` is false it silently behaves as dynamic.
10. **Julia doc pages contradict the code in at least one place**: `percolation.md` says
    water in a massless layer is routed to the next layer down; the code sends it straight to
    runoff. Follow the code, and note the discrepancy upstream.
11. **The depth cap ignores the configured `Ntot`**, using the hard-coded reference count 15.
12. **`EnergyWorkspace` has two unused arrays** and reuses `previous_temperature` as the
    diagonal scratch. Only 5 work arrays are needed.

---

## 6. Execution order

```
WP1 ──► WP2 ─────────────────────────────────────────────────────────► WP18
  └──► WP3 ──┬──► WP4 ──┐
             ├──► WP5 ──┤
             ├──► WP6 ──┼──► WP8 ──┐
             ├──► WP7 ──┘          │
             ├──► WP10 ────────────┼──► WP11 ──► WP13 ──► WP14 ──► WP15 ──► WP16 ──► WP19
             ├──► WP12 ────────────┤                                    ┌──► WP17 ──┘
       WP9 ──────────────────────--┘                                    │
```

Parallel batches:

| batch | WPs | agents |
|---|---|---|
| 0 | WP1 → {WP2, WP3} | 1 sequential, then 2 parallel |
| 1 | WP4, WP5, WP6, WP7, WP9, WP10, WP12 | 7 parallel |
| 2 | WP8 | 1 |
| 3 | WP11 → WP13 → WP14 → WP15 | 1 sequential (tightly coupled) |
| 4 | WP16, WP17, WP18 | 3 parallel |
| 5 | WP19 | 1, gated on WP16 sign-off |
| — | WP20 | rolling, alongside everything |

WP4 and WP5 are the two that decide whether the port succeeds; if only two WPs get careful
review, make it those. WP9 is the one most likely to surface upstream bugs rather than port
bugs — treat a clean WP9 with suspicion.
