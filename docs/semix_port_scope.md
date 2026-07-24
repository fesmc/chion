# Porting SEMIX's surface scheme into chion — staged scope

Goal: make chion span the whole snowpack-model design space — from minimal-input
bulk melt models up to CLIMBER-X's SEMIX surface energy balance with
spectral, dust-aware albedo — as **selectable, orthogonal options** on top of the
existing firn/layer physics, which stays untouched. The near-term target is
**offline physics equivalence** validated on Greenland; the eventual swap of
chion in for `semi` inside CLIMBER-X is a separate, later effort (see
[Deferred: CLIMBER-X integration](#deferred-climber-x-integration)).

SEMIX itself (CLIMBER-X's `src/smb/`, `semi_m`/`smb_ebal_mod`; Willeit, Calov,
Ganopolski) is characterized in the memory note and in
[docs/steady_state_snowpack.md](steady_state_snowpack.md) (model comparison
section). This document is the *how*.

## Design: orthogonal flags, not a model list

Layering and surface scheme are independent axes. `model = bessi` becomes the
**energy-balance family**, configured by flags; `pdd`/`itm` remain the separate
bulk melt-parameterization family (no energy balance), unchanged.

| axis | flag | values | status |
|---|---|---|---|
| column structure | `Ntot` *(exists)* | `1` (single layer) … `N` (firn column) | ✅ works today |
| surface energy balance | `seb_scheme` *(new)* | `bessi` \| `semix` | ✅ rungs 2–3 done |
| albedo | `albedo_scheme` *(extend enum)* | `constant` \| `dynamic` \| `prescribed` \| `semix` | rung 1 (dust) |
| net shortwave | `has_q_sw_net` *(exists)* + internal spectral | prescribed **or** chion-owns | both supported |
| background ice albedo | `has_prescribed_ice_albedo` *(new, optional)* | chion's own **or** passed | rung 1 |

Every SEMIX-ward permutation is then a flag combination:

| column × SEB | `bessi` SEB | `semix` SEB |
|---|---|---|
| **N-layer** (`Ntot=N`) | today ✅ | ✅ |
| **1-layer** (`Ntot=1`) | today ✅ (2.1× faster, R² 0.86) | ✅ |

### The 1-layer column is already here (`Ntot=1`)

No new flag. `Ntot=1` runs chion's dedicated single-layer closed-form energy
solve (`snow_energy.f90:315-341`), not a degenerate multilayer run. On GRL-16KM
it is the fastest config (22.5 s, 2.1× vs n=15) and the best MAR match (bias
−2.3, R² 0.86) — see the layer-count comparison in
[steady_state_snowpack.md](steady_state_snowpack.md). The single layer is an
evolving-density *firn* layer (in-layer refreezing/densification), not SEMIX's
fixed-density bulk reservoir; that difference is physical, not structural.
`Ntot=1, seb_scheme=semix` is therefore the near-1:1 SEMIX analog and costs
nothing beyond the SEB port. A strict fixed-density bulk variant is a possible
later option, not required for flexibility or a close comparison.

### Why `seb_scheme` is a sub-flag, not a new `model`

The whole SEB is a single call, `snow_energy_flux(...)` at `snow_bessi.f90:667`,
and the firn conduction solver (`solve_tridiagonal_thomas`, assembly at
`snow_energy.f90:352-393`) is fully reusable. Everything surface-specific is
steps 1–2 of that routine (`snow_energy.f90:242-313`); everything
conduction-specific is steps 3b–6 (`:343-464`), joined only through `rhs(1)` /
`diag(1)`. So a SEMIX surface scheme swaps the *flux formulas* feeding `rhs(1)`
/ `diag(1)` and leaves the firn untouched — mirroring how chion already carries
`albedo_scheme` and `low_density_densification` as in-BESSI sub-schemes.

### Coupling decision: α (confirmed)

chion has **no separate skin node** — the top firn-layer temperature *is* the
surface, and the SEB is linearized into row 1 of the firn conduction matrix.
SEMIX's `ebal` uses a **massless skin node** solved analytically
(`t_skin = num/denom`, `smb_ebal.f90:132-136`), coupled to its subsurface by a
ground heat flux. Both decompose every flux into the same constant + linear-in-Tₛ
form, so:

**α — match SEMIX's flux formulas inside chion's row-1 linearization.** SEMIX's
per-flux `num`/`denom` map directly onto chion's `q_const`/`q_lin`
(`snow_energy.f90:306-307`). No new node, firn solver reused, minimal surgery.
The top layer keeps a small thermal inertia vs SEMIX's massless skin —
negligible at daily steps, and chion's diurnal substepping already resolves
sub-daily. (Rejected: β, porting SEMIX's massless-skin + ground-flux boundary
wholesale — deeper surgery, changes the conduction top BC, partly duplicates the
existing top layer.)

## Extension mechanics (reference for every rung)

**New optional forcing field** — ~7–8 templated lines across two files
(`chion_defs.F90`, `chion_model.f90`), all behind a `has_*` flag so minimal-input
mode keeps working:
1. value + flag in `chion_step_forcing_class` (`chion_defs.F90:237-252`)
2. `(ncol)` value + flag arrays in `chion_forcing_class` (`:274-290`)
3. allocate both in `chion_forcing_alloc` (`:499-515`)
4. defaults (value + `.FALSE.`) in `chion_forcing_alloc` (`:528-544`)
5. deallocate in `chion_forcing_dealloc` (`:561+`)
6. per-column scatter in `chion_pack_step_forcing` (`chion_model.f90:123-138`)

**New scheme enum** (mirror `albedo_scheme`): enum constants in `chion_defs.F90`,
a string→flag helper like `chion_albedo_scheme_flag` (`:735-763`), loaded in
`chion_const_load` (`chion_api.f90:761+`). Branch at the one relevant call site.

**Shared constants to add** (SEMIX pulls these from CLIMBER-X `constants`):
`Ls`, `q_sat_i`/`dqsat_dT_i` (saturation over ice), `rho_a(T,p)`, `cap_a`,
`karman`, `g`, `z_sfl`, roughness `z0m_*`, `zm_to_zh`. Land in `chion_defs.F90`
constants / a small `snow_surface_par` analog.

---

## Rung 1 — Spectral + dust albedo *(build first)*

The priority: dust darkening for glacial inception/deglaciation. Independent of
the SEB rungs — it can land first and be used with the existing `bessi` SEB.

**Port** (`smb_surface_par.f90`): `surface_albedo` (4-band driver, `:45-134`),
`snow_grain_size` (`:143-172`), `dust_in_snow` (`:179-218`), and one of
`snow_albedo_ww` (Warren & Wiscombe, `:225-281`) or `snow_albedo_dang` (Dang
2015, `:288-389`). Start with WW; add Dang behind an `isnow_albedo` sub-option.

**New albedo scheme:** `CHION_ALBEDO_SEMIX` enum value + branch in
`snow_albedo.f90`. Note the existing branch tests `== CHION_ALBEDO_CONSTANT` and
falls through to dynamic (trap 9) — the new value needs an explicit branch.

**New prognostic per-column state** (add to `bessi_state_class`,
`snow_bessi.f90:116-148`; today only `t_srf`/`albedo` are instantaneous scalars):
- `snow_grain(:)` — grain size / aging proxy (needs `t_skin`, snowfall rate)
- `dust_con(:)` — dust concentration in snow (needs `dust_dep`, snowfall, and
  `w_snow_max` seasonal-max SWE for the melt-amplification term)
- `w_snow_max(:)` — seasonal max column SWE (drives dust melt amplification)
- `dt_snowfree(:)` — snow-free timer
- albedo state broadened scalar → **4 bands** {vis,nir}×{dir,dif}

**New forcing inputs** (each `has_*`-gated): `dust_dep`, `coszm`, `cloud`
(dir/dif weighting), `z_sur_std` (subgrid orography σ). Optional
`prescribed_ice_albedo` + `has_prescribed_ice_albedo` for the background (else
chion's own ice albedo — the default).

**Net shortwave becomes spectral.** The two SW-absorption sites
(`snow_energy.f90:262`, `snow_surface_fluxes.f90:431`) sum band contributions
instead of one `(1-albedo)` scalar. chion owns the albedo→net-SW loop (needs
spectral SW↓ inputs — CLIMBER-X already carries `swd_sur_{vis,nir}_{dir,dif}`).
The prescribed path (`has_q_sw_net`) remains as the "pass swnet directly"
alternative, so both can be compared; at daily steps the feedback is weak enough
that prescribing net SW is a fair approximation.

**Validate:** GRL-16KM, `Ntot=1` and `Ntot=15`, with vs without dust forcing;
sensitivity of margin/percolation SMB to `dust_con`. (No Greenland dust obs to
match against — validate behaviour and, later, against SEMIX output.)

**Effort: L.** Commits: (1) enum + scalar→band albedo state + net-SW spectral
sum, dust off; (2) `snow_grain_size` + aging; (3) `dust_in_snow` + `w_snow_max`
state + `dust_dep`/`coszm`/`cloud` inputs; (4) WW bands; (5) optional prescribed
ice albedo; (6) Dang variant.

## Rung 2 — Aerodynamic turbulent fluxes

**Port** `resistance` (`smb_surface_par.f90:396-438`): snow-weighted roughness,
neutral `Ch`, bulk-Richardson stability, `r_a = 1/(Ch·wind)`. Replace chion's
bulk `D_sh` sensible + latent with `f_sh = ρa·cap_a/r_a`, `f_lh = Ls/r_a·ρa`,
plus SEMIX's dew inhibition (`l_dew`: zero the latent flux when `q2m > qsat`).

**Inputs:** `wind` (exists), humidity (`relative_humidity`/`has_*` exists → derive
`q2m`), `air_pressure` (exists); roughness params `z0m_snow`/`z0m_ice`/`zm_to_zh`
as new constants. `h_snow` from the column.

**Gate under `seb_scheme=semix`** (first piece of that branch). `bessi` SEB keeps
`D_sh`.

**Validate:** GRL-16KM, `seb_scheme=semix` vs `bessi`, both `Ntot=1/15`.

**Effort: M.**

### Rung 2 result (done)

`seb_scheme` swaps the turbulent exchange at **all three** flux sites, not just
the energy solve: the linearized surface row (`snow_energy.f90`), the exact
bare-ice fluxes and the post-solve vapour mass (`snow_surface_fluxes.f90`).
Wiring only the first would leave bare-ice ablation — much of the margin melt —
on `D_sh` while the snow column used `r_a`.

`f_sh = ρa·cp_air/r_a` is **not** a near-equivalent of `D_sh = 10`:

| regime | `r_a` [s m⁻¹] | `f_sh` [W m⁻² K⁻¹] |
|---|---|---|
| stable (surface 5 K colder than air) | 168 | 8.0 |
| unstable (surface 5 K warmer) | 67 | 19.9 |

(deep snow, 5 m s⁻¹, 263 K air, CLIMBER-X defaults). SEMIX leaves the stable
branch at the neutral coefficient and only enhances the unstable one.

Over Greenland the **stable** side dominates — a melting surface sits at T0
under warmer summer air — so the scheme damps sensible heating and reduces
melt. GRL-16KM, 50 yr, surface SMB vs MAR:

| config | bias | RMSE | R² | melt [mm/yr] |
|---|---|---|---|---|
| `bessi`, `Ntot=1` | −2.3 | 197 | 0.86 | 204 |
| `semix`, `Ntot=1` | +7.0 | 204 | 0.85 | 197 |
| `bessi`, `Ntot=15` | +18.3 | 214 | 0.84 | 203 |
| `semix`, `Ntot=15` | +28.9 | 227 | 0.82 | 195 |

The cost concentrates at the margin (`z<800`: bias +57 → +126, R² 0.62 → 0.56
at `Ntot=1`); the interior is untouched. Two caveats on reading this as a
verdict on SEMIX: the latent flux is **identically zero** in these runs
(`has_relative_humidity` is false in the GRL driver), and `D_sh = 10` is a
tuned BESSI value while `z0m_snow`/`z_sfl` are CLIMBER-X's. This is a
faithful-port checkpoint, not a calibration.

**Composability:** `seb_scheme` is one switch among several — running a full
SEMIX configuration means setting `albedo_scheme`, `seb_scheme` and their
sub-options together, which is the intended design (orthogonal axes, not a
model list).

Two sub-options land with it: `semix_qsat` (`"semix"` = CLIMBER-X's `q_sat_i`,
`"bessi"` = chion's ice vapour pressure through the same 0.622/p; they agree to
0.1%), and `l_neutral`/`l_dew` carried over from `smb_par`.

**Humidity forcing (`rh_default`).** No domain loader carries a humidity field,
so `has_relative_humidity` was false everywhere and the turbulent latent flux
was identically zero in every benchmark to date — under *both* schemes. The
grid driver now has an `rh_default` knob (mirroring `wind_default`), off at
zero. GRL-16KM, `Ntot=1`, all ice:

| config | bias | RMSE | R² | melt | subl |
|---|---|---|---|---|---|
| `bessi`, rh off | −2.3 | 197 | 0.86 | 204 | 0.0 |
| `bessi`, rh=0.7 | +0.5 | 219 | 0.83 | 167 | 47.0 |
| `semix`, rh off | +7.0 | 204 | 0.85 | 197 | 0.0 |
| `semix`, rh=0.7 | −9.7 | 222 | 0.82 | 168 | 54.2 |

A uniform 0.7 is not a good humidity field — R² falls in both schemes — but the
latent path is large (≈50 mm/yr sublimation, melt down by a fifth) and was
previously untested end-to-end. Caveat: the two schemes read the same number
differently, BESSI relative to saturation over **water** and SEMIX over **ice**,
so varying `rh_default` across `seb_scheme` is not a controlled comparison of
the turbulent exchange alone.

**One asymmetry is deliberate.** SEMIX builds `f_lh` with the latent heat of
sublimation at every temperature, so under `seb_scheme=semix` the post-solve
vapour *mass* conversion uses `Lv+Lm` unconditionally, while the *reservoir*
choice (solid `mass(1)` vs liquid `mass_w(1)`) still turns on T0. Under `bessi`
both still turn on T0, as before. Bare ice already used `Lv+Lm`
unconditionally, so this makes the snow and bare-ice budgets agree.

## Rung 3 — Full SEMIX flux set → `seb_scheme=semix`

Complete the SEMIX surface balance under α: emissivity-based longwave (snow vs
ice `emissivity_*`), the exact `num`/`denom` decomposition of net SW, LW,
sensible, latent, ground flux from `ebal` (`smb_ebal.f90:100-145`), folded into
chion's `q_const`/`q_lin`. Wire the `update_tskin` second-correction equivalent
if needed (chion's re-solve at the melting point, `snow_energy.f90:400-454`,
already plays this role — confirm the mapping).

After this rung, `seb_scheme=semix` is complete and composes with any `Ntot` and
either albedo scheme.

**Validate:** term-by-term flux comparison against a SEMIX single-column run if
feasible; GRL-16KM skill vs the `bessi` SEB.

**Effort: M.**

### Rung 3 result (done)

Rung 3 turned out to be **one substantive change**, not four. Taking the
`ebal` flux set term by term:

| `ebal` term | status under coupling α |
|---|---|
| `num_sw` | already identical — chion's `sw_abs` is the same quantity |
| `num_sh`/`denom_sh` | landed in rung 2 |
| `num_lh`/`denom_lh` | landed in rung 2 |
| `num_g`/`denom_g` | **no analog** — chion's row-1 conduction coupling is the ground flux |
| `num_lw`/`denom_lw` | **the actual rung-3 work** |
| `update_tskin` | **nothing to port** — see below |

**Longwave.** SEMIX absorbs the downwelling flux with the surface emissivity,
`emiss·lwdown`; BESSI takes it at face value and applies `eps_snow` only to the
emitted term — an absorptivity of 1 against an emissivity of 0.98, which is not
Kirchhoff-consistent. That asymmetry is the whole of rung 3, and it is worth
about **5 W m⁻² less energy into the surface**, in one direction, all year.

SEMIX also carries separate snow and ice emissivities. chion reuses `eps_snow`
and adds `eps_ice`, both defaulting to chion's 0.98; CLIMBER-X uses 0.99 for
both. The *value* barely matters — raising it increases absorption and emission
together — but the asymmetry does:

| config (`Ntot=1`, all ice) | bias | RMSE | R² | melt | margin bias |
|---|---|---|---|---|---|
| `bessi` | −2.3 | 197 | 0.86 | 204 | +57 |
| `semix`, eps = 0.98 | +20.0 | 211 | 0.84 | 187 | +189 |
| `semix`, eps = 0.99 | +21.5 | 212 | 0.84 | 186 | +196 |

**`update_tskin` has no chion analog, and that is a structural consequence of
coupling α, not an omission.** SEMIX solves its massless skin node *before* the
subsurface step, so once `smb_temp` has updated `t_prof` the skin temperature is
stale and must be re-diagnosed against the new ground flux with `flx_melt`
removed. chion has no separate skin node: `t_srf` is `temperature(1)`, which
comes out of the *same* implicit solve as the conduction, simultaneously — there
is nothing to go stale. The melting-point re-solve
(`snow_energy.f90:400-454`) already plays the role of `update_tskin`'s
`- flx_melt` term, dropping the surface-flux feedback from row 1 and pinning it
at T0. Confirmed, not ported.

### Cumulative rungs 2+3

GRL-16KM, 50 yr, surface SMB vs MAR:

| config | bias | RMSE | R² | melt |
|---|---|---|---|---|
| `bessi`, `Ntot=1` | −2.3 | 197 | 0.86 | 204 |
| `semix`, `Ntot=1` | +20.0 | 211 | 0.84 | 187 |
| `bessi`, `Ntot=15` | +18.3 | 214 | 0.84 | 203 |
| `semix`, `Ntot=15` | +42.9 | 243 | 0.79 | 186 |

Both rungs push the same way — less energy into the surface, less melt, more
positive SMB bias — and both concentrate at the margin (`z<800`, `Ntot=1`: bias
+57 → +189, R² 0.62 → 0.51). The interior is untouched at every rung.

This is a *faithful-port* checkpoint and should not be read as SEMIX performing
worse than BESSI. The two schemes are not being compared on equal terms:
`D_sh = 10` and `eps_snow = 0.98` are tuned BESSI values, the roughness and
surface-layer parameters are CLIMBER-X's untuned defaults, the latent flux is
zero without `rh_default`, and SEMIX in CLIMBER-X runs with its own spectral
albedo rather than chion's `dynamic` one. A like-for-like comparison needs
`albedo_scheme=semix` and a real humidity field.

## Rung 4 — SEMIX diurnal / statistical melt *(optional)*

SEMIX resolves sub-daily melt statistically from `tstd` (daily T std dev) and
`swnet_min` via the Krapp et al. 2016 cycle (`smb_ebal.f90:153-243`). chion
already has diurnal *shortwave substepping* (a different, arguably better
mechanism). Add the tstd scheme as an alternative only if SEMIX-equivalent melt
statistics are wanted. Inputs: `tstd`, `swnet_min`.

**Effort: M.** Low priority.

---

## Composed: the full SEMIX configuration

The rung-by-rung numbers above each vary one axis. This is the 2×2×2 with both
schemes and humidity together, GRL-16KM, `Ntot=1`, 50 yr, surface SMB vs MAR.
`rh` means `rh_default = 0.7`; `dry` means no humidity forcing.

| albedo | SEB | rh | bias | RMSE | R² | melt | subl | wall [s] |
|---|---|---|---|---|---|---|---|---|
| dynamic | bessi | dry | −2.3 | 197 | **0.86** | 204 | 0.0 | 26.5 |
| dynamic | bessi | rh | +0.5 | 219 | 0.83 | 167 | 47.0 | 35.5 |
| dynamic | semix | dry | +20.0 | 211 | 0.84 | 187 | 0.0 | 38.1 |
| dynamic | semix | rh | +7.5 | 236 | 0.80 | 160 | 48.2 | 42.2 |
| semix | bessi | dry | −9.4 | 199 | **0.86** | 212 | 0.0 | 33.5 |
| semix | bessi | rh | −4.6 | 224 | 0.82 | 167 | 53.6 | 41.9 |
| semix | semix | dry | +17.1 | 213 | 0.84 | 191 | 0.0 | 44.7 |
| semix | semix | rh | +0.9 | 241 | 0.79 | 161 | 56.6 | 48.4 |

(The `semix`/`bessi`/`dry` row reproduces the standalone rung-1 reference,
−9.4 / R² 0.86, exactly.)

### The effects are not additive, and the reason is a measurement artifact

Summing the single-axis deltas predicts a full-configuration bias of +15.7;
the actual value is +0.9. The interaction is almost all in the humidity column,
where the sign of the humidity effect flips with the SEB scheme (+2.8 and +4.8
under `bessi`, −12.5 and −16.2 under `semix`).

That flip is **not physics**. The two schemes read `rh_default` against
different saturation references — BESSI over water, SEMIX over ice — so at the
same nominal 0.7 SEMIX is given systematically drier air:

| T | es_water/es_ice | rh over ice matching rh_water = 0.7 |
|---|---|---|
| 273.15 K | 1.000 | 0.700 |
| 268.15 K | 1.059 | 0.741 |
| 263.15 K | 1.121 | 0.785 |
| 253.15 K | 1.259 | 0.881 |

Re-running the full configuration at `rh_default = 0.785` (the ~263 K
equivalent) removes the gap almost entirely:

| config | all ice | margin<800 | 800–1500 | 1500–2200 | interior>2200 |
|---|---|---|---|---|---|
| `dynamic`/`bessi`, rh=0.7 | +1 (0.83) | +189 (0.54) | +72 (0.68) | −10 (0.93) | −50 (0.96) |
| `semix`/`semix`, rh=0.7 | +1 (0.79) | +249 (0.45) | +104 (0.61) | −18 (0.92) | −67 (0.92) |
| `semix`/`semix`, rh=0.785 | −2 (0.83) | +179 (0.54) | +75 (0.67) | −9 (0.93) | −57 (0.94) |

**With comparable humidity forcing the full SEMIX configuration is
indistinguishable from BESSI** — R² 0.83 both, RMSE 221 vs 219, and matching
zone for zone. No single scalar can correct this properly (the ratio is 1.0 at
T0 and 1.26 at −20 °C), so 0.785 is a bound rather than a fix; the honest
reading is that the rh=0.7/0.785 spread, R² 0.79–0.83, **brackets a forcing
ambiguity as large as the scheme differences it was being used to measure.**

### What the matrix actually says

- **Dry, the SEMIX SEB costs a little skill**: R² 0.86 → 0.84, entirely at the
  margin. That result is clean — no humidity is involved.
- **The SEMIX albedo alone costs nothing and helps the margin**: R² 0.86
  either way, margin bias +57 → +28, the only component that improves it.
- **Uniform humidity is the dominant error source, in both schemes**: R² 0.86 →
  0.83 under `bessi`, and it is what drags the composed configuration down.
  A real humidity field matters more than either scheme choice.
- **Near-zero bias is not skill.** The full configuration at rh=0.7 has a global
  bias of +0.9 — better than the baseline's −2.3 — while being the *worst* R² of
  the set, with margin +249 against interior −67. The mean cancels; the pattern
  does not.
- **Cost**: the full configuration is 1.8× the baseline (48.4 s vs 26.5 s).

At `Ntot=15` the full configuration gives bias +19.5 / R² 0.72 (`rh=0.7`,
uncorrected, so pessimistic on the same grounds). The layer-count axis continues
to behave as the rung-2/3 tables show — `Ntot=1` remains the better match.

## Deferred: CLIMBER-X integration

Out of near-term scope (offline physics first). When chion is fully capable, the
swap is an adapter, not a physics change:

- A driver-layer adapter (like `libs/domains/`) maps CLIMBER-X's `smb_in_class`
  (2-D, post-downscaling) → chion column list → `chion_update` → back to
  `s_out`/`ts_out`; replace `call semi` at `smb_model.f90:520`.
- **Stays host-side in CLIMBER-X, feeds chion as inputs:** all downscaling
  (`downscaling.f90`), bias correction (`smb_bias_corr.f90`), ice-albedo/`f_ice`
  (`ice.f90`), topography (`topo.f90`), grid remapping, `fake_atm_hires`. This is
  exactly chion's "host supplies per-column forcing" contract.
- **Caveat:** SEMIX computes albedo *before* it builds net SW (albedo feeds the
  spectral→net-SW assembly, `semi.f90:210-236`). With chion owning the dust
  albedo, either chion receives spectral SW↓ and closes the loop internally
  (preferred), or the host runs an albedo→swnet pre-step. Both supported via the
  net-shortwave axis above.

## Validation harness

Every rung uses the existing standalone Greenland path
(`scripts/run_layer_comparison.sh` / `diagnostics/compare_models.jl`, SMB vs MAR,
global + per-elevation-zone). Run each rung at **both `Ntot=1` and `Ntot=15`** so
the surface-scheme and albedo effects are separated from the layering effect.
Paleo/dust behaviour (rung 1) has no Greenland obs target and is validated on
sensitivity and, later, against SEMIX output directly.
