# Porting notes: Chion.jl -> chion

Every intentional deviation from `Chion.jl` (branch `main`), and every Julia quirk
deliberately preserved. Required reading before changing any physics module.

Policy for what may and may not be cleaned up: `docs/PLAN.md` section 4.1.
The full list of traps: `docs/PLAN.md` section 5.

Each entry states **what**, **why**, and **impact**.

---

## WP1 — `chion_defs.f90`

### D1. `wp = sp` for state, `wp_acc = dp` for accumulators
**What:** Chion.jl is `Float64` throughout. chion uses `wp = sp` for the layered state,
forcing and all public interfaces, and `wp_acc = dp` for the nine cumulative per-column
accumulators (`smb_ice`, `runoff`, `melt`, `refreezing`, `vapor_mass`, `sublimation`,
`latent_heat_flux_sum`, `mass_base`, `pdd_sum`). Three expressions are additionally evaluated
in `dp` locally: pore volume, cold content, and surface-energy accumulation across diurnal
substeps.

**Why:** `wp = sp` matches yelmo and fesm-utils, so nothing has to be converted at the yelmox
boundary. The split was measured rather than assumed — see `docs/PLAN.md` section 3.1 for the
full table. In summary:

- The tridiagonal conduction solve is *safe* in `sp` (max 3.4e-5 K vs `dp`; the matrix is
  strictly diagonally dominant). This was the expected risk and turned out not to be one.
- Cumulative accumulators are *not* safe: 0.01 kg m-2 increments onto a 1e5 kg m-2 total are
  lost outright in `sp`, drifting 80 kg m-2 (0.08%) over a 100-year daily run. This is the
  only failure with real physical impact, and it is why `wp_acc` exists.
- Pore volume `phi = m/rho - m/rho_i` is quantized at ~1e-5 m in `sp`, so the `TOL_TINY`
  guard can never fire, and `phi` feeds the division `lwc = m_w/rho_w/phi`.
- Cold content `(T0-T)*ci*m_s` is pure noise below 3e-5 K from melting, though the implied
  mass error is only ~1e-5 kg m-2. Computed in `dp` because it costs nothing.

**Impact:**
- No conversion needed at the yelmo/yelmox interface. `wp_chion` is still exported (mirroring
  yelmo's `wp_yelmo`) so host code can be explicit.
- **All conservation tolerances are `1e-6` relative, not `1e-12`.** Checks over `wp_acc`
  quantities may still use `dp` thresholds; every WP states which applies.
- Any new cumulative quantity added later must be `wp_acc`. This is the easiest thing in the
  port to get wrong, because it compiles and runs fine and only shows up as slow mass drift.

### D1b. Both tolerances are declared `dp`
**What:** `TOL_TINY` and `TOL_EMPTY_LAYER` are `real(wp_acc)`, not `real(wp)`.
**Why:** So comparisons promote correctly. `TOL_TINY = 1e-12` is below `sp` resolution for
anything of order 1 or larger, so `x <= TOL_TINY` is only a meaningful test when `x` itself
was computed in `dp`. The acceptance test asserts both halves of this explicitly.
**Impact:** Guards must be applied to `dp`-computed quantities. Where a guard protects an
`sp` quantity, it is effectively a test against zero — which is what Julia gets too, just for
a different reason.

### D2. Tolerances renamed `EPS_*` -> `TOL_*`
**What:** `EPS_TINY` -> `TOL_TINY`, `EPS_EMPTY_LAYER` -> `TOL_EMPTY_LAYER`.
**Why:** Consistency with yelmo, which uses `TOL` / `TOL_UNDERFLOW`. `EPS` collides
conceptually with machine epsilon, which these are not.
**Impact:** Naming only. Values and usage are unchanged, and the two remain distinct — see
trap 1 in `docs/PLAN.md` section 5.

### D3. Scheme flags kept inside the constants type
**What:** `albedo_scheme`, `fresh_snow_density_scheme` and `low_density_densification` live in
`chion_const_class`, not in `chion_param_class`, even though they are configuration rather
than physical constants.
**Why:** Mirrors Chion.jl's `SnowpackPhysicalConstants`. Every physics routine already
receives `c`; moving the flags would mean threading an extra argument through roughly ten
routines for no behavioural gain.
**Impact:** Cosmetic. Revisit if the parameter/constant split becomes confusing in WP13.

### D4. Scheme flags are integers; namelist input is strings
**What:** Julia stores `UInt8` flags and constructs them from `Symbol`s. chion stores plain
`integer` parameters (`CHION_ALBEDO_*` etc.) and converts from namelist strings via
`chion_albedo_scheme_flag` and friends.
**Why:** Integer branching in inner loops; readable, validatable namelist input.
**Impact:** None. The Chion.jl aliases are preserved: albedo `bessi`/`legacy` -> `constant`;
fresh-snow-density `bessi` -> `constant`, `htessel` -> `parameterized`.

### D5. Unrecognized densification scheme is an error, not a silent fallback
**What:** Chion.jl dispatches densification with `if _uses_htessel_densification(c) ... else`
(`src/processes/densification.jl:262`), so any value other than HTESSEL silently selects
BESSI. `chion_densify_scheme_flag` rejects unknown names with a message and `stop`.
**Why:** A typo in a namelist should not silently change the physics.
**Impact:** Behaviour differs only for inputs that were already invalid. Trap 8 in
`docs/PLAN.md` section 5.

### D6. Time metadata is scalar in `chion_forcing_class`
**What:** `day_of_year` and `solar_longitude_deg` are scalars, not `(ncol)` arrays. Julia
stores them per column and per timestep.
**Why:** They are properties of the timestep, not the column. Latitude, which genuinely varies
per column, remains an array.
**Impact:** None for any current use. If a host ever needs per-column time (it should not),
this becomes an array.

### D7. `chion_class` is not defined in `chion_defs`
**What:** The top-level `chion_class` will be defined in `chion_api.f90` (WP11), not here.
**Why:** It contains the model state types (`bessi_class`, `pdd_class`, `itm_class`), which
are defined in their own modules. yelmo puts everything in `yelmo_defs`; chion keeps model
state with the model.
**Impact:** `chion_defs` stays physics-free and compiles standalone, which is what makes the
Level 1 work packages independently testable.

---

## WP2 — build system

### D9. configme has no local package-override tier — WP18 was a blocker
**What:** `configme -m macbook -c gfortran` initially failed with *"package 'chion' is not
supported"*. configme resolves packages only from its own shipped `data/packages/`; unlike
machines and compilers, there is **no user-level or repo-level override tier**
(`src/configme/data.py:34`).
**Why it matters:** WP18 was therefore a *blocker* for WP2's acceptance test, not a
downstream task as the plan originally had it. Registering a new FESM package requires a
change to the configme repo plus a reinstall before it takes effect — worth knowing for any
future package.
**Resolution:** `src/configme/data/packages/chion.toml` added on branch `add-chion-package`
in `~/models/configme` (commit `a5656da`); full configme test suite passes (210 tests).
After `pip install -U`, `configme -m macbook -c gfortran` generates the root `Makefile`
correctly and `make chion-static` builds against it.
**Still open:** the branch is pushed (`fesmc/configme:add-chion-package`) but unmerged, and
chion is not yet listed in `data/orchestrators/yelmox.toml` `default_packages` — that belongs
to WP19 phase A.

### D11. `.gitignore` patterns must anchor `Makefile` to the repo root
**What:** `/Makefile`, not `Makefile`.
**Why:** An unanchored pattern matches at every level and silently excludes
`config/Makefile`, which is the configme *template* and must stay tracked. This was caught
only because the file was missing from a commit summary. yelmo and yelmox both use the
unanchored form and get away with it because their templates were tracked before the ignore
rule existed — do not copy that pattern into a new repo.

---

## WP3 — `snow_column_utils.f90`

### D8. Kernels take contiguous column slices, not (array, index)
**What:** Julia physics kernels take the full `(Ntot,ncol)` array plus a column index `idx`,
and reach elements through `_get_layer(mass, k, idx)`. chion passes a column slice —
`mass(:,icol)` — plus the active layer count `n`.
**Why:** The Julia indirection exists only so one kernel body can run over a `Matrix` on GPU
and a `Vector` on CPU; chion has no GPU path. Because the arrays are `(Ntot,ncol)` and
Fortran is column-major, a slice is contiguous and is passed by reference with no copy. The
resulting signatures are far more readable and make OpenMP privacy obvious by inspection.
**Impact:** Applies to every Level 1 physics module, so it is a convention, not a local
choice. The outer dispatcher still works in `(state, icol)` terms and slices at the call
site. Routines that change the layer count (split/merge) additionally take `n` as
`intent(INOUT)`.

### D10. Fortran does not short-circuit `.or.`
**What:** `is_lowest_active_snow_layer` is written as a nested `if`, not as the single
expression `k >= n .or. mass(k+1) <= 0`.
**Why:** Julia's `||` short-circuits, so `mass(k+1)` is never evaluated when `k == n`.
Fortran leaves evaluation order unspecified, and `k == n == Ntot` would read one element past
the array. The acceptance test covers exactly this case.
**Impact:** None behaviourally; it is a correctness requirement. The same hazard will recur
anywhere a Julia guard relies on `||`/`&&` short-circuiting — check for it in every ported
predicate.

---

## Preserved deliberately (do NOT "fix")

These are listed in full in `docs/PLAN.md` section 5. Restated here as they are encountered:

- **Three distinct empty-layer thresholds** (`> 0`, `> TOL_TINY`, `> TOL_EMPTY_LAYER`) gate
  different physics and are not interchangeable. `chion_defs` defines both tolerances and
  comments the hazard at the declaration.
- **`BESSI_REFERENCE_LAYER_COUNT = 15`** is used by the depth cap regardless of the configured
  `Ntot`. Preserved and flagged at the declaration.

---

## WP11 — public API

### D13. `chion_get_smb` returns one reconciled ice-facing rate in `[kg m-2 s-1]`
**What:** the three models disagree natively on what "SMB" means. BESSI's `smb_ice` is
ice-only; PDD's `smb_ice` is a whole-column change; ITM reports both `smb` (whole-column) and
`smbi` (ice-facing), as a per-step rate in mm/d. `chion_get_smb` returns a single quantity for
all three: **net mass flux to the ice sheet, `[kg m-2 s-1]`, positive = ice gains mass,
averaged over the step just completed.**

| model | source quantity | note |
|---|---|---|
| BESSI | `smb_ice` | already ice-facing (bottom export, bare-ice vapour, bare melt, ice melt) |
| PDD | `smb_ice - snowpack_swe` | identically `-ice_melt`; removes the whole-column reservoir term (defect D2) |
| ITM | `smbi_cum` | smbpal's `snow_to_ice + refrz - melted_ice`; **not** its `smb`, which is whole-column |

`chion_class` retains `smb_cum_prev` and `dt_last`, refreshed at the top of every
`chion_update`, so the host tracks nothing.

**Why:** the consumer is an ice-sheet model, which needs a flux to apply over its own
timestep, not a running total. SI matches the units `chion_forcing_class` already uses for
`snowfall_rate`/`rainfall_rate`, so every mass flux across the boundary is in one unit.
Returning `m i.e. yr-1` would bake chion's `sec_year`/`rho_ice` into the host's mass budget;
the header documents the conversion for a host that wants yelmo's convention.

**Impact:** all three mappings are exact identities on the models' own bookkeeping —
`sum(smb*dt_seconds)` reproduces the cumulative accumulator to `sp` round-off, asserted in
the test. Two consequences worth knowing before WP19:

- **PDD's ice-facing flux is never positive**, structurally. PDD has no densification and no
  snow-to-ice conversion, so the only ice-facing term is `-ice_melt`. Not a bug.
- **PDD's recovery differences an `sp`-stored reservoir**, so it carries unbiased round-off
  bounded by `eps_sp*snowpack_swe` per step (~4e-10 kg m-2 s-1 at 300 kg m-2). It errs in both
  directions and cancels over a run; the only visible effect is that a purely accumulating
  PDD column may report ~-1e-10 rather than exactly 0. Upstream fix: add a `wp_acc` cumulative
  `ice_melt` to `pdd_state_class`.
- BESSI and PDD legitimately return exactly 0 for a fresh accumulating column. An ice-facing
  flux is not the surface accumulation rate.

---

## WP14/WP15 — IO and drivers

### D14. Output dimension order differs from Chion.jl, of necessity
**What:** Chion.jl declares `("t","x","y")`, landing in the file as `var(y,x,t)`. chion writes
`var(time,yc,xc)`.
**Why:** `time` is the unlimited dimension, and netCDF requires the unlimited dimension to be
slowest-varying. This is also CF-standard and yelmo's convention.
**Impact:** **WP16's comparison harness must permute axes.** Variable names, units and
`long_name` are unaffected and match exactly for all 20 shared variables.

### D15. Unmapped grid cells use `MV = -9999`, not NaN
**What:** Chion.jl fills cells with no column at NaN; chion uses `chion_defs`' `MV`, written as
both `missing_value` and `_FillValue`.
**Why:** a NaN in an `sp` field is indistinguishable from one produced by an FPE-trapped
operation under `debug=1`, which would make a real bug look like ordinary masking.
**Impact:** WP16 must treat `MV` and NaN as equivalent when comparing.

### D16. Restart is host-driven, not folded into `chion_init_state`
**What:** `chion_restart_read` is a separate call. `par%restart` is loaded from the namelist
and consumed by the host or driver, not by `chion_init_state`.
**Why:** `chion_io` depends on `chion_api` for `chion_class`, so `chion_api` cannot call
`chion_io` without a circular dependency. This is the same pattern yelmox uses for smbpal.
**Impact:** the documented sequence is in the `chion_init_state` header. `use chion` exposes
both. If `chion_class` ever moves to a lower module, this can be folded in.

### D17. `chion_write_step` recomputes BESSI's four diagnostics
**What:** `thickness`, `wet_mass`, `bulk_density` and `liquid_water` are recomputed via
`summarize_domain_state` at write time rather than read from `chn%bsi%now`.
**Why:** `bessi_column_step` does not update them (Chion.jl produces them on the way to
output), so the stored values are stale by construction.
**Impact:** the *restart* still writes the stored values, because `bessi_reset_columns`
deliberately does not reset them (upstream defect 23) and that staleness is part of what an
exact restart must reproduce.

### D18. The shipped column example was retuned to exercise the model
**What:** `par/chion_column.nml` originally described a strongly ablating column: peak
accumulation stayed below `mass_max`, so the layer machinery never fired (peak layer count 1)
and there was no refreezing. Retuned to a percolation-zone firn site — `t2m_mean = 258.15`,
`t2m_amplitude = 16`, `pr_mean = 3.0e-5`, `sw_mean = 165`, `sw_amplitude = 145`.
**Why:** a shipped example that exercises one code path is a poor demonstration and a weak
WP16 comparison target.
**Impact:** BESSI over 10 years now reaches 12 layers, densifies 315 -> 651 kg m-3, and
produces 5532 kg m-2 of melt, 2437 of refreezing and 4109 of runoff — so accumulation,
split/merge, all three densification regimes, melt, percolation and refreezing are all
covered. Colder or drier settings build layers but never melt; warmer ones ablate to bare ice
and never split.

---

## WP16 — validation harness

### D22. `8.13` in the densification Arrhenius denominators is a typo; corrected
**What:** `DENSIFY_R_GAS` is the universal gas constant, not Chion.jl's `8.13`.
Reported upstream as Chion.jl issue #18.

**Why:** provenance. The low-density branch is **Herron & Langway (1980)**,
*Firn densification: an empirical model*, J. Glaciol. 25(93), stage 1:
`k0 = 11*exp(-10160/(R*T))` with `R = 8.314 J K-1 mol-1`. The activation energy
`10160 J mol-1` matches this module exactly and is H&L's published value for
rho < 550; Chion.jl's leading `0.011` is their `k0 = 11` with the kg/Mg
conversion; `Ea/(R*T)` is dimensionless and equals 4.70 at 260 K. The expression
is only dimensionally coherent with `R`. The mid/high branches carry
`Ea = 60000 J mol-1`, the standard ice-creep value, which expects the same `R`.

PLAN.md §4.1 previously listed this as NOT reconcilable without a modelling
decision, because the substitution was measured at "a factor of ~2.6" at 260 K
and "~1400" on the mid/high branches. **Those figures were wrong.** Recomputed:

| branch | at | factor |
|---|---|---|
| low-density (`-10160`) | 260 K | **1.11** |
| mid/high (`-60000`) | 263 K | **1.86** |

Both branches move in the same direction and by a similar order, which is what a
single shared constant should do. The earlier asymmetry was an artefact of the
arithmetic error, and it was itself the main argument for treating `8.13` as
calibrated rather than mistyped.

**Impact:** the *integrated* effect is small, because densification self-limits
against the `(rho_i - rho)` gap and the `rho_i` cap. Over a 10-year single-column
percolation-zone run (the §4.1 measurement requirement, now discharged):

| | `8.13` | gas constant | change |
|---|---|---|---|
| final bulk density | 651 kg m-3 | 657.0 | +0.9% |
| cumulative refreezing | 2437 | 2456.3 | +0.8% |
| cumulative runoff | 4109 | 4114.8 | +0.14% |
| cumulative melt | 5532 | 5531.4 | -0.01% |
| peak layer count | 12 | 12 | — |

So it is a correction, not a recalibration. Against Chion.jl it is a 1-9%
divergence in the density-driven fields, which is why D24 exists.

### D23. PDD adopts BESSI's ice-facing `smb_ice` convention and a capped reservoir
**What:** chion's PDD no longer reproduces Chion.jl's snowpack budget. `smb_ice`
is ice-facing (`snow_to_ice - ice_melt`); refreezing is capacity-limited by the
snow remaining after melt (`min(f*H_snow, snow_melt)`); refrozen mass becomes
superimposed ice and leaves the melt-able reservoir; and the reservoir is capped
at the new `H_snow_max` parameter (default 5000 kg m-2, matching smbpal/ITM),
with the excess converted to ice. Reported upstream as Chion.jl issue #19.

**Why:** Chion.jl's PDD credits `smb_ice` with `d(snowpack_swe)`, making it a
whole-column mass change despite its own metadata calling it "Net mass forcing
to the ice sheet". That double-counts against a host ice-sheet model. PDD has no
firn representation, so a whole-column `smb_ice` implies a reservoir the model
does not have. BESSI and ITM are both ice-facing; PDD was the odd one out. The
capped one-layer scheme is what smbpal and Chion.jl's own BESSI already use.

**Impact:**
- **The full three-reservoir closure now holds:**
  `snowfall + rainfall == d(snowpack_swe) + d(smb_ice) + d(runoff)`.
  Correction C2 below records that this identity *cannot* hold — that was true
  of Chion.jl's convention and is no longer true of chion's. It is now the WP16
  gate for PDD, measured at 2.3e-15 relative in dp and 3.8e-07 in sp.
- `chion_get_smb`'s PDD special case (D13, `smb_ice - snowpack_swe`) is deleted.
  All three models now share one definition of the ice-facing flux, and the
  sp round-off that case carried is gone with it.
- **PDD can now report a positive ice-facing flux.** D13 recorded that it was
  structurally never positive; the cap's snow-to-ice export is what fixes that.
- `snowpack_swe` is bounded, so `wp` is safe for it. Uncapped it was not.
- Three test assertions were inverted, because they had pinned the upstream
  defects in place — most starkly a "D6 guard" asserting that the snowpack is
  *never exhausted* under sustained melt.
- `pdd_column_apply` / `pdd_column_step` gained optional `snow_melt_out`,
  `ice_melt_out`, `refrozen_out`, `snow_to_ice_out` diagnostics. The test used
  to infer ice melt as `d(snowpack_swe) - d(smb_ice)`, an identity of the old
  convention that became silently wrong when the convention changed.

### D24. `legacy_chion=1` build variant
**What:** `make ... legacy_chion=1` defines `CHION_LEGACY`, which reverts the
deliberate physics corrections (currently D22) to Chion.jl's values. Builds land
in `libchion/{include,bin}[-dp]-legacy`.

**Why:** "is the port faithful?" and "is the reference correct?" are different
questions. Without this, every upstream bug chion fixes turns into a WP16 gate
failure, and the only way to keep the gate green is to stop testing those fields
— so the harness gets weaker exactly as the port gets better. D22 alone would
have ungated 15 of BESSI's 18 fields.

**Impact:** WP16 gates BESSI port fidelity with the dp+legacy build, and reports
the correction's effect separately as chion-dp vs chion-dp-legacy. **Not a
production setting** — it selects physics believed to be wrong.

D23 is deliberately NOT covered by the switch: reproducing Chion.jl's PDD
convention would mean maintaining a second copy of the PDD core, which is
upstream defect 13 (three diverged copies) reintroduced on purpose. PDD is
gated on its own mass closure instead, which is a stronger property than
agreement with a reference that cannot satisfy it.

### D21. `chion_grid.x` stamps output at the end of the step, not the start
**What:** the driver wrote the post-step state under the pre-step time, and its
output test was seeded such that with `dt_out == dt` the after-step-1 record was
never written at all — the first output record held the state after step *two*.
Both are fixed: the stamp is `time + dt_use`, and the test is written on the same
post-step time so `dt_out <= dt` emits a record after every step.

**Why:** labelling a post-step state with the pre-step time shifts the whole
output series one step earlier than the physics, which biases any comparison
against a reference model by exactly one timestep — silently, and in a way that
looks like a small physics disagreement rather than a bookkeeping error.
`chion_column.x` was already correct (`time = time_init + k*dt`), so only the
gridded driver was affected. The restart stamp was corrected to match.

**Impact:** output files from `chion_grid.x` before this change have their time
axis shifted by one step and are missing the first post-step record. Nothing
else consumed them yet.

### D19. `wp` is selectable at compile time (`precision=sp|dp`)
**What:** `wp` is no longer fixed at `sp`. `make ... precision=dp` defines `CHION_DP`, and
`src/chion_defs.F90` (renamed from `.f90`; the capital `F` makes both gfortran and ifort
preprocess it with no `-cpp`/`-fpp` flag and no configme change) selects `wp = dp`. `wp_acc`
stays `dp` in **both** builds — the accumulator argument in D1 is about summation over
10^4–10^5 steps, not about the precision of the state, and it holds regardless of `wp`.

**Why:** WP16 compares chion against Chion.jl, which is `Float64` throughout. With `wp = sp`
the two models share discrete branch points (`mass_max` split, `mass_min` merge,
`surface_has_snow`, the melt gate) but not the round-off that decides which side of them a
given step lands on, so a *correct* port diverges by O(1) in layer structure once any branch
fires a step apart. Validating at `wp = dp` removes that confound: a residual is then
attributable to the port rather than to the precision. Running both builds against the same
Chion.jl reference additionally turns the sp-vs-dp difference into a measured number, which
was previously assumed (PLAN.md §3.1 measured expressions in isolation, never end-to-end).

**Impact:**
- `sp` keeps the unsuffixed `libchion/include` and `libchion/bin`, so the artifact paths in
  `.configme/manifest.toml` and `.runme/info.json` are unchanged and `dp` is purely additive
  (`libchion/include-dp`, `libchion/bin-dp`). The two builds coexist without `make clean`.
- Pre-existing and *not* addressed here: `openmp=0` and `openmp=1` still share the sp
  directories, so switching openmp does still need a clean. Only fesm-utils' include dir is
  swapped by that flag.
- chion declares no generic interfaces, so `wp == wp_acc` in a dp build creates no ambiguous
  overloads internally. It does change ncio resolution: `nc_write` selects its **double**
  variant for `wp` arrays, so a dp build writes double-precision NetCDF variables. Any reader
  must not assume float.
- All 13 acceptance tests pass in all six configurations (sp/dp × `-O2`/`debug=1`/`openmp=1`).
- Tests that assert precision policy must now assert it *per build*. `test_defs` previously
  hard-coded `wp == sp`, which does not verify the policy — it forbids one of two supported
  builds.

### D20. `test_itm`'s smbpal reference is typed `_wp`, not copied verbatim
**What:** the reference implementation inside `tests/test_itm.f90` was a byte-faithful copy of
`smbpal/src/smb_itm.f90`, deliberately retaining smbpal's untyped and dp literals (`0.d0`,
`1.d0`, `1d3`, `0d0`, and a bare `273.15`). Every floating literal is now typed `_wp`.

**Why:** the verbatim form silently tested a *different model* under `precision=dp`. A bare
`273.15` parses as **single** precision — `273.1499938964844` — and is then widened, so
against chion's true-dp `273.15_wp` the reference carried a 6.1035e-06 K offset in the ITM
temperature term. That propagated to ~1.8e-3 mm/d in melt and, on one step, pushed `melt`
across `melt_crit`, flipping the albedo between `alb_snow_dry` and `alb_snow_wet` — a 1.0e-5
jump. Seven of the nine fields failed their tolerance. The predicted offset matched the
measured `tsrf` difference to every printed digit, which is what identified the cause.

Under `wp = sp` the untyped literals *are* the working precision, so the issue is invisible;
it is purely an artefact of the reference, not of chion. This is upstream smbpal defect 5
(mixed sp/dp literals) having a real consequence rather than a stylistic one.

**Decision:** the test now asserts that chion's ITM is **algebraically identical** to smbpal's,
at whatever precision the build uses. What it gives up is the ability to detect divergence
from smbpal-as-compiled. That was judged an acceptable trade because smbpal is not itself a
validated reference — see the smbpal defect list below, in particular defect 1, which biases
SMB by up to −2100 mm w.e./yr in the ablation zone.

**Impact:** the measured residual fell by roughly two orders of magnitude at sp, confirming
that the previous 3.2e-6 was smbpal's literal artefacts and *not* a porting difference:

| build | worst rel(scale) | in ulp of `eps(wp)` |
|---|---|---|
| sp | 1.14e-07 (`refrz`) | 0.96 |
| dp | 5.09e-15 (`refrz`) | 22.9 |

`alb_s`, `melt` and `tsrf` are now bit-identical at sp. The gate was re-derived accordingly,
from 1e-5 to `128*epsilon(wp)` at each field's own scale. Note the sp build shows *fewer* ulp
than dp: sp rounds coarsely enough that most of these differences fall below the last stored
bit and cancel to exactly zero, while dp resolves them instead of discarding them — more ulp
at a far smaller absolute error.

---

## WP-wide build note

### D12. `debug=1` drops `underflow` from the FPE trap set
**What:** `config/Makefile` overrides the compiler fragment's `DFLAGS_DEBUG`, removing
`underflow` from `-ffpe-trap`. yelmo and the shipped configme fragment both trap it.
**Why:** Gradual underflow to zero is correct and expected throughout this model — `exp()` of
a large negative argument appears in the HTESSEL thermal metamorphism term, the PDD
normal-CDF tails, the albedo decay law, and both densification Arrhenius factors. Trapping it
turns correctly-rounded results into SIGILL. Three of the ten acceptance tests (WP4, WP7,
WP9) aborted under `debug=1` while passing under `-O2`, purely from this.
**Impact:** `invalid`, `zero` and `overflow` are still trapped, which is where real bugs
surface. Remove the override to restore yelmo's flag set. Note the agents developing batch 1
tested against `invalid,zero,overflow` only, so this mismatch was invisible to all of them —
worth checking any future WP against the *project's* debug flags, not a hand-rolled set.

---

## Batch 1 acceptance-criteria corrections

Two work-package briefs specified criteria that were wrong. Recorded because the reasoning
matters more than the fix.

### C1. WP5 — "the analytic linear gradient" does not exist
The brief asked for a steady state under constant surface flux with a zero-flux bottom. No
such steady state exists: with no sink, energy accumulates and the column warms indefinitely.
The correct statement is the *quasi*-steady profile, in which the whole column warms at a
uniform rate and the conductive flux across interface `k` is `F*(1 - M_k/M_tot)` — the flux at
depth carries only the energy needed to warm the mass below it. A genuine constant gradient
would need a Dirichlet base, which the scheme deliberately lacks, and faking one with a
massive bottom layer distorts `interface_conductance` because `dz = m/rho`.

### C2. WP9 — the requested mass-balance identity cannot hold *(superseded by D23)*
**Superseded.** This was true of Chion.jl's convention, which chion reproduced
at the time. chion's PDD is now ice-facing (D23) and the full identity
`d(smb_ice) + d(runoff) + d(snowpack_swe) == snowfall + rainfall` holds exactly
and is asserted. The analysis below remains accurate as a description of
Chion.jl.

The brief asked to assert
`d(smb_ice) + d(runoff) == snowfall + rainfall - d(snowpack_swe)`.
Chion.jl credits `smb_ice` with `d(snowpack_swe)`, so `smb_ice` is a whole-column mass change,
not the "net mass forcing to the ice sheet" its own NetCDF metadata claims. **BESSI uses the
opposite convention** (`smb_ice` receives only bottom-export, bare-ice and ice-melt terms), so
the two models in Chion.jl disagree about what their headline output variable means — and
`smb_ice` feeds `ice_sheet_net_forcing_yearly` directly. What does hold, and is asserted, is
`d(smb_ice) + d(runoff) == snowfall + rainfall`, plus a regression guard that the residual
equals `d(snowpack_swe)` exactly.

---

## Open items to report upstream to Chion.jl

Collected across batch 1. Severity: **A** = wrong results, **B** = latent/conditional,
**C** = cosmetic or doc-only.

### The BESSI mass-closure identity (WP8)

Derived by enumerating every mass mutation in `column_step_core!` and splitting runoff into
column-sourced and ice-sourced parts:

```
total_snow_water_mass + runoff + smb_ice - vapor_mass  ==  cumulative accepted precipitation
```

**`mass_base` does not appear** — it is already inside `smb_ice`
(`smb_ice = mass_base + vapor_bare - melt_bare - melt_ice` by construction), so including it
separately double-counts. Worth knowing before writing any conservation check in WP11 or WP16.

Two conditions are required for it to close, both of which are upstream defects rather than
port artefacts: rain must be withheld on steps beginning with `mass(1) <= 0` (defect 11), and
humidity forcing must be off (defect 1). With dry air over a thin pack, defect 1 alone leaves
a 105 kg m-2 residual against 192 kg m-2 of reported sublimation.

Measured relative residual: 9.2e-7 (BESSI densification), 8.8e-7 (HTESSEL), 6.9e-7 to 9.0e-7
across all twelve scheme combinations.

**The residual saturates with run length** — 6.4e-7 at 1 yr, 9.2e-7 at 5 yr, 9.7e-7 at 20 yr,
9.8e-7 at 40 yr. It is a bounded per-step relative bias (~8 ulp of `sp`, from
`mass(1) = m_prev + m_added`), not a random walk, so `sp` is safe here. But it sits just
inside the 1e-6 acceptance threshold with almost no headroom: **that tolerance cannot be
tightened without moving the layer mass arrays to `dp`.**

### A — defects

24. **(A) A CF-compliant numeric time axis is silently discarded.**
    `_read_time_values` (`dataloaders.jl:17`) reads `ds[time_name].var[:]` — the
    **undecoded** NCDatasets accessor — and then tests `eltype(raw) <: DateTime`.
    NCDatasets only returns `DateTime` from the *decoded* accessor `ds[name]`,
    so that test can never be true for a numeric `days since ...` axis, however
    correctly it is written. The function falls through to a synthetic axis of
    **one day per record** (`dataloaders.jl:38`), discarding the file's real
    timing with no warning. Only the `YYYY/MM/DD/HH` branch actually fires.

    Not cosmetic: `dt_days` feeds every rate-to-mass conversion, and the
    PISM/Calov–Greve PDD integral is auto-selected on `27 <= dt_days <= 32`. On
    a 30-day forcing file Chion.jl therefore stepped 30x smaller than chion
    *and* silently used the simple PDD form instead of PISM — every PDD field
    wrong by a factor ~30. Found in WP16; the fix is to read `ds[name][:]`, or
    to convert the numeric axis using its `units` attribute.

25. **(B) NetCDF output is written in Float32 while the model computes in
    Float64.** The output buffers in `io.jl` are `Matrix{Float32}` /
    `Array{Float32,3}`. Roughly seven of the sixteen significant digits carried
    through the calculation are discarded on write, which sets a hard floor of
    ~1.2e-7 relative on any comparison, regression test or restart made through
    these files — including Chion.jl's own. Worth a `Float64` option at least
    for validation output. Found in WP16, where it turned out to be the binding
    constraint on how tightly the port can be verified.

26. **(C) Output files carry no time coordinate variable** — only a bare `t`
    dimension. A Chion.jl output file cannot be interpreted on its own; the
    reader has to already know the forcing that produced it. Found in WP16.

1. **(A) Vapor-mass diagnostics are not mass-closed.**
   `_apply_snow_surface_vapor_mass_flux!` returns the *unclipped* `vapor_mass` while the mass
   it applies is clipped by `max(..., 0)`. When sublimation demand exceeds the surface layer,
   the cumulative `vapor_mass`/`sublimation` diagnostics overstate what was removed.
2. **(A) `free_slot_for_surface_split` reads index 0** when the column has no active layers —
   `_get_layer(mass, _n_active(...), idx)` with no guard. chion raises an explicit error.
3. **(A) `_htessel_thermal_metamorphism` is dead above ~150 kg m-3.**
   `xi = 2.8e-6*exp(-4.2e-2*(T0-T) - 460*max(0, rho-150))` — the `460` multiplies a density
   excess in kg m-3, so the exponential underflows to exactly zero for any density above
   ~150.15. The term only ever acts on the freshest snow. Looks like a missing unit
   conversion (460 per Mg m-3?).
4. **(A) PDD on GPU has never run.** Both GPU call sites in `runtime.jl` pass one argument
   too many; no matching method exists. Separately, `active_indices` is silently ignored on
   every CPU PDD path, while BESSI honours it.
5. **(A) PDD refreezing is uncapped** — a flat fraction of all snow melt refreezes with no
   cold-content or capacity limit. smbpal caps at `min(melt_snow, acc*f_refrz_max)`, where
   capacity scales with *accumulation*. Biases SMB positive, unbounded in combination with 7.
6. **(A) PDD refrozen mass re-enters the melt-able reservoir**, so the snowpack decays only
   as `0.4^n` and is never exhausted.
7. **(A) PDD `snowpack_swe` is unbounded** — no cap, no aging, no densification, so the
   ablation buffer depends on spin-up length.

19. **(A) Diurnal substepping multiplies the albedo aging rate by `n_substeps`.**
    `_update_surface_albedo_arrays!` sits inside `column_step_core!`, and the aging law
    carries no `dt` (trap 5). With `max_substeps = 8` the albedo ages eight times per day.
    Bounded by the `alpha_wet` floor, but it means enabling substepping silently changes the
    albedo scheme, not merely the shortwave resolution. Found in WP8.
20. **(A) The bare-ice path uses `rainfall_rate` in the energy budget but discards its mass.**
    Extends defect 11: rain is a genuine mass leak on *any* bare column, not only on
    massless-surface columns. Found in WP8.

### B — latent

8. **(B) `merge_surface_layer` divides by `subsurface_mass` unguarded.** Safe only because the
   defaults make the divisor positive (600 vs 100); a namelist with
   `mass_min > 2*mass_split` divides by zero.
9. **(B) `enforce_snow_depth_cap` handles zero-density layers inconsistently** between its
   forward and reverse passes — excluded from the depth sum, but the full mass is exported
   without reducing the remaining depth demand, so the cap can over-export.
10. **(B) `continuous_bottom_deplete`'s empty-layer skip drops a layer without consuming
    demand**, so a column of near-empty layers unwinds completely under an arbitrarily small
    depletion request.
11. **(B) Rain on a bare column is silently dropped.** `_apply_accumulation_resolved!` adds
    rain only when `mass[1] > 0` strictly and creates no layer for rain-only forcing, so the
    routine is not mass-closed on its own.
12. **(B) `_state_dict` divides `mass ./ density` with no guard** (diagnostics.jl:32), so a
    zero-density active layer yields `Inf` in `thickness` and `total_thickness`. The kernel
    path is guarded; only the snapshot path is exposed.
13. **(B) Three copies of the PDD core have already diverged** — the `ddf_snow > 0` guard
    exists in only one, and the vectorised methods take `ddf_snow` positionally, bypassing
    constructor validation.
14. **(B) `_apply_melt!`'s fast path returns the full request** while the general path returns
    the amount actually melted. `step.jl:346` relies on the difference to trigger bare-ice
    melt.

21. **(B) Snowfall brightening is non-linear under substepping.** `1 - exp(-dm/3)` means eight
    brightenings of `dm/8` do not equal one of `dm`, so substepping darkens fresh snow
    relative to the daily-mean path. Same root cause as 19. Found in WP8.
22. **(B) `workspace.liquid_water_before_energy` is written for `1:n_before` and read for
    `1:n_after`.** Safe only because `n` cannot grow between the snapshot and the compaction.
    It is a persistent buffer, so a violation would read another column's data rather than
    obvious garbage. chion uses a zero-initialised stack-local snapshot instead.
23. **(B) `_reset_bessi_columns_kernel!` does not reset the four diagnostics**
    (`thickness`, `wet_mass`, `bulk_density`, `liquid_water`) that
    `_initialize_bessi_state_kernel!` does, so a deactivated column carries stale diagnostics
    into output. Preserved, not fixed.

### C — doc and cosmetic

15. **(C) `percolation.md` contradicts `percolation.jl`** on where liquid water in a massless
    layer goes: the doc says the next layer down, the code sends it straight to runoff.
16. **(C) `_column_has_liquid_water` docstring** says "above the empty-layer tolerance" but
    the code uses `EPS_TINY`.
17. **(C) `_surface_liquid_water_content` guards `mass[1]` on `EPS_TINY`** while
    `_update_surface_albedo_arrays!` guards on `EPS_EMPTY_LAYER` one line earlier. Trap 1
    territory, but this particular pair is almost certainly unintentional.
18. **(C) `step.jl:383` recomputes `has_liquid_water` after refreezing and never reads it.**
    Dropped as dead code in chion.
19b. **(C) The `_n_active == 0` guard on the ice-melt shortfall is redundant** — `_apply_melt!`
    can only return `melted < melt_mass` when the column is empty.
20b. **(C) PDD hard-codes `273.15` and `86400.0`** rather than using a constants struct, and
    uses a single `sigma` where smbpal uses three by surface type.

The full PDD analysis, with quantified impacts and recommended fixes, is in
`docs/pdd_defects.md` (12 Chion.jl defects, 1 smbpal defect).

---

## Open items to report upstream to smbpal

1. **(A) `calc_ablation_pdd` counts snow that was never there.** It sets
   `melt = melt_snow + melt_ice` with `melt_snow` the *uncapped potential*, so `melt`,
   `runoff` and `smb` are all wrong by `melt_snow - acc`, always in the ablation direction.
   **This matters for WP12 and WP19**, which use smbpal as the equivalence target — we would
   be validating against a buggy reference. Must be raised before the yelmox cutover.
2. **(B) `calc_albedo_surface` tests `H_ice .eq. 0.0` exactly**, so `H_ice = 1e-6 m` selects
   the ice branch rather than land.
3. **(B) `H_snow_crit` has no zero guard**; `H_snow_crit_desert = 0` divides by zero.
4. **(B)** The `else` branch of `melt_net` computes `refrz - melted_snow` with a comment
   asserting `refrz` is zero there. It is not enforced, and with `H_ice = 0` and snow present
   `refrz` can be nonzero.
5. **(B, was C)** Mixed sp/dp and untyped literals (`0.d0`, `1.d0`, `1d3`, and a bare
   `273.15`) compile only as a gfortran extension. **Upgraded from C to B:** this is not
   cosmetic. If smbpal is ever compiled at `wp = dp`, the bare `273.15` still parses as
   single precision and is then widened, putting a 6.1e-06 K offset into the ITM temperature
   term — enough to shift melt by ~1.8e-3 mm/d and to flip the surface albedo across
   `melt_crit`. smbpal's results are therefore silently precision-dependent beyond ordinary
   round-off. Found while porting chion to a dp build; see D20.
