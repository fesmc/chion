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

### D9. configme cannot install chion until the registry entry lands
**What:** `configme -m macbook -c gfortran` fails with *"package 'chion' is not supported"*.
configme resolves packages only from its own shipped `data/packages/`; unlike machines and
compilers, there is **no user-level or repo-level override tier** (`src/configme/data.py:34`).
**Why it matters:** WP18 (adding `configme/data/packages/chion.toml`) is therefore a
*blocker* for WP2's acceptance test, not a downstream task as originally planned. It also
requires a change to a separate repo plus a `pip install -U` before it takes effect.
**Interim:** `config/Makefile` was verified by assembling it exactly as configme does —
compiler fragment, then machine fragment, then the auto-detected netCDF block, then
`include config/common.mk` — and asserting the placeholder appears exactly once and the
`common.mk` include is present. Default, `debug=1` and `openmp=1` builds all verified.
**Status:** open; needs a decision on committing to the configme repo.

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

## Open items to report upstream to Chion.jl

- `docs/src/processes/percolation.md` contradicts `src/processes/percolation.jl` on where
  liquid water in a massless layer goes: the doc says the next layer down, the code sends it
  straight to runoff.
- `_column_has_liquid_water` docstring says "above the empty-layer tolerance" but the code
  uses `EPS_TINY`, not `EPS_EMPTY_LAYER`.
- `PDDModel` defects — see `docs/PLAN.md` section 3.1; the list is to be completed in WP9.
