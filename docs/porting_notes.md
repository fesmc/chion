# Porting notes: Chion.jl -> chion

Every intentional deviation from `Chion.jl` (branch `main`), and every Julia quirk
deliberately preserved. Required reading before changing any physics module.

Policy for what may and may not be cleaned up: `docs/PLAN.md` section 4.1.
The full list of traps: `docs/PLAN.md` section 5.

Each entry states **what**, **why**, and **impact**.

---

## WP1 — `chion_defs.f90`

### D1. `wp = dp`, not re-exported from fesm-utils
**What:** `chion_defs` imports `sp`/`dp` from fesm-utils `precision` but defines its own
`wp = dp`. fesm-utils (and yelmo) use `wp = sp`.
**Why:** Chion.jl is `Float64` throughout, and `EPS_TINY = 1e-12` is below single-precision
resolution — in `sp` it would be indistinguishable from zero, silently disabling several
guards. Decision recorded in `docs/PLAN.md` section 3, item 2.
**Impact:** Host code linking against both yelmo and chion must convert at the interface.
`wp_chion` is exported for that purpose, mirroring yelmo's `wp_yelmo`.

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
