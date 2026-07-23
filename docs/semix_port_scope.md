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
| surface energy balance | `seb_scheme` *(new)* | `bessi` \| `semix` | rungs 2–3 |
| albedo | `albedo_scheme` *(extend enum)* | `constant` \| `dynamic` \| `prescribed` \| `semix` | rung 1 (dust) |
| net shortwave | `has_q_sw_net` *(exists)* + internal spectral | prescribed **or** chion-owns | both supported |
| background ice albedo | `has_prescribed_ice_albedo` *(new, optional)* | chion's own **or** passed | rung 1 |

Every SEMIX-ward permutation is then a flag combination:

| column × SEB | `bessi` SEB | `semix` SEB |
|---|---|---|
| **N-layer** (`Ntot=N`) | today ✅ | rungs 2–3 |
| **1-layer** (`Ntot=1`) | today ✅ (2.1× faster, R² 0.86) | rungs 2–3, free at `Ntot=1` |

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

## Rung 4 — SEMIX diurnal / statistical melt *(optional)*

SEMIX resolves sub-daily melt statistically from `tstd` (daily T std dev) and
`swnet_min` via the Krapp et al. 2016 cycle (`smb_ebal.f90:153-243`). chion
already has diurnal *shortwave substepping* (a different, arguably better
mechanism). Add the tstd scheme as an alternative only if SEMIX-equivalent melt
statistics are wanted. Inputs: `tstd`, `swnet_min`.

**Effort: M.** Low priority.

---

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
