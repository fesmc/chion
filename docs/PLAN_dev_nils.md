# Plan: sync chion with Chion.jl `dev_nils`

Status: decisions recorded (§4), nothing implemented. Written 2026-10-09.

## 0. Baseline

| | commit | date |
|---|---|---|
| chion port base (PLAN.md "v0.2.0 @ main"; `validation/Manifest.toml` develops `~/models/Chion.jl` by path) | `a9ec154` | 2026-07-09 |
| `origin/main` (= dev_nils up to PR #26) | `3ac2525` (merge of `12407a3`) | 2026-08-10 |
| `origin/dev_nils` (current) | `27113b6` | 2026-09-01 |
| `origin/nils-devel` | `4d0ff00` — **stale**, ancestor of dev_nils, 0 commits of its own | 2026-03-23 |

- `a9ec154..origin/dev_nils`: 37 non-merge commits. `origin/main` has nothing that dev_nils lacks (merge commits only); main lacks dev_nils' last 6 (`5b4fe49`..`27113b6`).
- The local checkout `~/models/Chion.jl` is still on `a9ec154` (39 behind origin/main) with untracked `src/SUMupData.jl`, `test/Project.toml`. Do not move it; use a worktree (WP0).

## 1. Upstream changes, classified

15 commits change results or the API; the rest are GPU/perf, scripts, gitignore or docs
(`1637b9c cc0bf6b a11a7dc eb091f6 b8f8ccf 03a5805 2900e56 e740f67 208dcf7 811f074 397a86f 75af9bf 5fb25e9 fc56827 8f0d5c8 e50e42d 75a73d7 cd0f785 e5d38b9 13b7adc`),
plus `aba07a0` (type-stability refactor, bit-neutral) and `6fca5d7` (PDD honours active columns; chion always has).

### 1a. Upstream converged on chion — no Fortran physics change

| commit | change | chion |
|---|---|---|
| `408e91c` | gas constant 8.13 → **8.314** | D22 already uses 8.31446. `legacy_chion` still reverts to 8.13 → must become 8.314 (WP1). |
| `ce6a68d` | PDD rewritten to chion's D23 budget: capped reservoir `H_snow_max`, refrozen → ice, ice-facing `smb_ice`, `pdd_method` simple/pism, `erfc` | Same algorithm as `snow_pdd.f90:pdd_column_apply`. PDD can now be **gated** against Chion.jl (WP3). Defects D1–D8, D10, D12 of `pdd_defects.md` closed upstream; D11 (one sigma) still open on both. |
| `29eb867`, `f2c46d0` | `ITMModel` added, ported from chion | Identical to `snow_itm.f90:itm_step` **before D27** (`04ff915`): Julia's Tsrf uses the daily `melt_net` rate. ITM can be compared against Chion.jl (WP4). |
| `12407a3` | restores BESSI bare-ice early return | Fortran always had it (`snow_bessi.f90:620-649`). |
| `6d077c5` | dynamic albedo aging × `dt_days` | Fixes our upstream defect 19. Bit-identical at daily dt without substeps. Fortran lacks it (WP5). |
| `385d452` → `5b4fe49` | percolation `max_lwc` 0.1 → 0.05, then back to 0.1 (inside an unrelated commit) | Net 0.1 = Fortran. Keep 0.1; ask Nils (N1). |

### 1b. Port (changes BESSI results or adds options)

| commit | change | Fortran target | default BESSI results | host contract |
|---|---|---|---|---|
| `8fff530` | bare-ice step: `runoff += rain + melt` (was melt only); closes defects 11/20 | `snow_bessi.f90:640`, `snow_accumulation.f90:19-23` note | runoff on bare-ice steps with rain | no |
| `6d06af6` | new albedo `:aging`: snowfall>0 → `alpha_dry`, age=0; else `α = α_wet + (α_prev−α_wet)·exp(−dt/τ)`, τ = 2 d if T_top ≥ T0 else 20 d; diagnostic `snow_age_days` (output) | `snow_albedo.f90` (new branch in `albedo_update` + refresh), `chion_defs.F90` flags/consts, nml, `chion_io.f90` | yes (default from WP14) | no |
| constants/models | defaults: `BESSIModel` albedo `:dynamic` → `:aging`; `alpha_wet` 0.70 → 0.60; `turbulent_flux_scheme` → `:semix` | `chion_defs.F90:533-557`, `input/chion_defaults.nml:119-143` | yes (WP14; turbulence excepted) | no |
| `27113b6` | `refreezing_correction` (scales cold content, default 1) + `min(T0,…)` clamp on full-refreeze temperature | `snow_refreezing.f90:89,113-117`, `bessi_par` | none at default (ulp near T0) | no |
| `81034fa` + `d0146e1` | interface conductance arithmetic → **harmonic** `G = 2KiKj/(Kj dzi + Ki dzj)`, with `β = −dt/ci` (the 2 was dropped in d0146e1; `81034fa` alone doubled diffusion, so its "doesn't change much" is void) | `snow_energy.f90:113-136, 418` | yes (sharp K contrasts) | no |
| `49990e6` | conductivity → **Calonne et al. (2019)**, T- and ρ-dependent; `Ki` becomes dead | `snow_energy.f90:97-111, 423, 430`; `Ki` in `chion_defs.F90:213,509`, `chion_api.f90:1048`, nml:98 | yes (+14 % at 500, +36 % at 700 kg m⁻³) | no |
| `d0146e1` | BESSI latent heat: L_v (not L_v+L_m) when T_s ≥ T0; vapour mass from the humidity gradient, independent of L (prescribed `q_lh` and semix still use Q/L(phase)) | `snow_vapor.f90:166-244`, `snow_surface_fluxes.f90:477-483` | yes at melting surfaces (~12 %) | no |
| `d0146e1` | `turbulent_flux_scheme` split from `seb_scheme` (seb → longwave/`eps_ice` only) | `snow_energy.f90:268`, `chion_defs.F90` | no (flag only) | no |
| `d0146e1` | diurnal T amplitude `A = clamp(A0 + γ·max(z_s − z_ref,0), 0, A_max)`; defaults neutral | `snow_bessi.f90:117,212,916`; pack `surface_height` into `chion_step_forcing_class` (`chion_model.f90:117` deliberately omits it) | no | internal only: host already fills `forc%surface_height` |
| `ecd4992`, `29eb867`, `6d06af6` | output: `t` coordinate (days since 1970, proleptic_gregorian), dims now (t,y,x) like chion; `snow_age_days`; ITM vars `alb_s`, `smb_cum`, `smbi`, `melt_net`, `*_cum` | `chion_io.f90` (comment :50-56 now stale), `input/chion-variables-itm.md` | no | no |

### 1c. Do NOT port — report upstream instead

1. **Julia's SEMIX turbulence (`5b4fe49`/`d0146e1`) is no longer CLIMBER-X SEMIX**: stability `1/(1+10Ri)` / `√(1−16Ri)`, wind floor 0.1, fixed z0m 0.001/0.01, z0h = z0m/10, z_sfl 10 m, phase-dependent L, no `l_dew`, BESSI q_sat, exchange factors. Neutral C_h ≈ +27 % vs ours. Porting it into `snow_seb_semix.f90` breaks `semix_port_scope.md`. Also internally inconsistent on bare ice (energy with L_v at T0, mass with L_v+L_m). Not ported for now; N3.
2. **Julia SEMIX albedo bugs** (ours matches CLIMBER-X): Dang direct-vis dust `−0.05` vs `−0.049`; WW 1–10 ppm band drops the `0.01` intercept and uses log 1 vs log 1.001; Dang `dalb_snow_vis/nir` never applied; `w_snow_max` reset on `day_of_year ≤ dt` (never fires with mid-day stamps; wrong for Antarctica); no host bare-ice albedo.
3. **ITM D27** (Tsrf scaled to annual `melt_net`) not upstream; Julia ITM firn warming ~360× too small.
4. `f2c46d0`: with prescribed albedo, a bare column keeps a stale albedo (Fortran applies it, `snow_bessi.f90:616-618`). Liquid-water guards removed (differs only for ≤1e-12 kg water) — keep Fortran's guards. Diurnal wrapper with `n_substeps==1` returns 0 SW at polar night.
5. `8fff530` double-count window: for `0 < mass(1) ≤ EPS_EMPTY_LAYER` rain is already in `mass_w(1)` and is added to runoff again.
6. `d0146e1`: `diurnal_temperature_cycle=true` with missing `surface_height` (NaN) gives NaN amplitude even with γ=0.
7. Output: ITM `melt/runoff/refreezing` written as rates but labelled "Cumulative"; `smb` name clashes with chion's ice-facing `smb`.
8. Stale upstream tests/docs: `test_case_api.jl:108,110` (`alpha_wet==0.70`, τ_melt 5), `albedo.md:85,105`.
9. Gravity 9.81 in densification (D25) still upstream.

## 2. Work packages

One commit each. "Gate" = `validation/validate.jl` dp+`legacy_chion=1` vs Chion.jl at 4 ulp Float32; "tests" = `make tests` all pass; "bit-id" = `chion_column.x` output identical to the previous commit (cmp of NetCDF data via the harness' `compare.jl`).

### Stage A — reference = `origin/main` (`12407a3`, before Calonne/harmonic/SEMIX/latent changes)

Tracking rule (decided): the reference is the `dev_nils` tip. Each re-point (WP0, WP8, any later sync) first runs
`git -C ~/models/Chion.jl log --no-merges origin/dev_nils..origin/main` and reports anything non-empty to the user. Empty as of 2026-10-09. Stage A uses `12407a3` only as an intermediate commit *on* dev_nils, so Stage A can be checked against Chion.jl before Stage B's three changes, which upstream gives no switch to turn off.

**WP0 — Re-point the harness (validation/ only).**
`git -C ~/models/Chion.jl worktree add <path> 12407a3`; `Pkg.develop` that path in `validation/`; re-pin `Manifest.toml` (Chion.jl now needs SpecialFunctions). Update `runners.jl` to the new API (typed `BESSIModel{…}`, `SnowpackForcing{calendar,fields}`, `t` coordinate now written → simplify record alignment). Pass explicit `albedo=:dynamic, alpha_wet=0.70, turbulent_flux_scheme=:bessi` so old chion physics is selected.
Verify: harness runs; record the field-by-field result. Expected failures: density-driven fields (8.13 vs 8.314, → WP1) and `runoff` on `bare_recover` (→ WP2). Everything else at ≤ 0.5 ulp as before.

**WP1 — `legacy_chion` gas constant → 8.314.** Production keeps our 8.31446261815324 (decided). Upstream's literal 8.314 differs by ~3e-4 relative in the Arrhenius term, far above the gate, so `CHION_LEGACY` (`chion_defs.F90:160-163`) reverts to 8.314 instead of 8.13. Once Nils switches to 8.31446 (N2), drop R from `CHION_LEGACY` and mark D22 "fixed upstream". Update D22/D24 text.
Verify: gate green except bare-ice runoff; production (non-legacy) build bit-id.

**WP2 — Bare-ice rain to runoff (`8fff530`).** `snow_bessi.f90:640`; remove the "rain dropped" note in `snow_accumulation.f90`; update the WP8 closure identity in `porting_notes.md` (no longer needs rain withheld). Decided: route rain exactly once — no double count when `0 < mass(1) ≤ EPS_EMPTY_LAYER` (1c.5). New D-entry; report upstream. The gate can differ from Chion.jl only in that window; assert the harness columns never enter it, or exclude it explicitly.
Verify: gate fully green at `12407a3`; BESSI closure test in `tests/test_bessi.f90` updated to include rain on bare columns and passes ≤1e-6.

**WP3 — PDD: default `simple`, gated against Chion.jl.** Decided: `pdd_method = "simple"` in `input/chion_defaults.nml:252` (and rewrite the "DEFAULT DIVERGES" note), matching Chion.jl. Move PDD from "reported" to "gated" in `validate.jl`, running both methods. Update `pdd_defects.md` "upstream" column and the README.
Verify: PDD fields within gate at dp for both methods; closure test unchanged; `par/` and test namelists that relied on the old default set `pism` explicitly where intended.

**WP4 — ITM compared against Chion.jl.** Add an ITM runner (Chion.jl needs `surface_height`, `ice_thickness`, `annual_pdd`, `latitude_deg` in the forcing). Put D27 under `legacy_chion` (decided; consistent with D24's rule) so `tsrf` can be gated too. Keep the smbpal comparison.
Verify: all ITM fields within gate at dp+legacy; production ITM bit-id (yelmox default unaffected).

**WP5 — Dynamic albedo aging × dt (`6d077c5`).** `snow_albedo.f90:203-205`; drop quirk 5 / PLAN §4.1 note and the nml defect-19 note (decided).
Verify: bit-id at daily dt, no substeps; gate green; new unit test: two half-day steps == one daily step (to round-off) in `tests/test_surface.f90`.

**WP6 — `:aging` albedo scheme (`6d06af6`).** New `CHION_ALBEDO_AGING`; params `aging_cold_timescale_days`, `aging_melting_timescale_days`; validation `0 ≤ α_wet ≤ α_dry ≤ 1`; `snow_age_days` state, reset on bare column and final fixup, written to output and restart (decided). Default unchanged here (switched in WP14).
Verify: default bit-id; gate with `albedo=:aging` added as a 5th BESSI configuration in the harness; unit test of the exponential relaxation.

**WP7 — `refreezing_correction` + clamp (`27113b6`) — DEFERRED (decided: wait).** Neutral at Chion.jl's default (1.0) up to an ulp near T0, so the Stage B gate is unaffected. Revisit if upstream changes the default.

### Stage B — reference = `origin/dev_nils` (`27113b6`)

**WP8 — Re-point harness to `27113b6`.** Same as WP0. Expected failures: all temperature-driven fields (Calonne, harmonic, latent heat).

**WP9 — Harmonic interface conductance + `β = −dt/ci`.** Together, one commit. `snow_energy.f90:113-136, 418`.
Verify: unit test — uniform column bit-id to old form up to round-off; analytic two-layer steady flux test in `tests/test_energy.f90`. Gate still red until WP10–11 (report diff reduction).

**WP10 — Calonne (2019) conductivity.** `snow_energy.f90:97-111`, layer-1 at T^n, others at their own T, `rho_i` from constants. Decided: replace outright, no switch (upstream has none; a switch would carry a dead `Ki` path that nothing validates). Remove `Ki` from `chion_defs`, `chion_api.f90:1048`, nml. Same for WP9: harmonic replaces arithmetic.
Verify: unit test against hand-computed values at (ρ,T) = (300,253), (500,253), (917,263); gate expected green only after WP11 (report remaining diff). Report 10-yr column change.

**WP11 — Phase-dependent latent heat; gradient-based vapour mass (BESSI turbulence).** `snow_vapor.f90`, `snow_surface_fluxes.f90:477-483`; bare-ice keeps L_v+L_m (as upstream).
Verify: gate fully green at `27113b6` (with explicit `albedo=:dynamic`, `alpha_wet=0.70`, `turbulent_flux_scheme=:bessi`); BESSI closure tests pass (re-check defect-1 vapour closure — note the harness runs humidity-off, so add a humidity-on unit test).

**WP12 — `turbulent_flux_scheme` flag split from `seb_scheme`.** Flags + parser; `seb_scheme` keeps longwave/`eps_ice`; `turbulent_flux_scheme = "bessi" | "semix"` where `semix` = our CLIMBER-X `snow_seb_semix.f90` (not Julia's, 1c.1).
Verify: all four combinations run; `seb_scheme=semix` (both) bit-id to current `semix`; `(bessi,bessi)` bit-id to current default.

**WP13 — Elevation-dependent diurnal T amplitude.** Add `surface_height` to `chion_step_forcing_class` and pack it in `chion_model.f90:117`; params γ [K km⁻¹], z_ref, A_max in `bessi_par`/nml. Guard the missing-height case properly (no NaN, 1c.6). No new host field: yelmox already fills `forc%surface_height`.
Verify: default bit-id; unit test of the clamp; diurnal stays out of the gate (harness note).

**WP14 — Adopt upstream defaults (decided).** `albedo = "aging"`, `alpha_wet = 0.60`. **Exception:** `turbulent_flux_scheme` stays `"bessi"` — upstream's default `:semix` is Julia's modified scheme, which we are not porting yet (§1c.1, N3). Harness drops the albedo overrides but keeps `turbulent_flux_scheme=:bessi`. Also re-check test/par namelists that assumed the old defaults.
Verify: gate green with Chion.jl defaults; record default-run deltas (10-yr column, ANT/GRL chion_grid) in CHANGELOG.

**WP15 — Output/naming alignment.** Decided: ITM output takes Julia's names (`alb_s`, `smb`, `smbi`, `melt_net`, `smb_cum`, `melt_cum`, `runoff_cum`, `refreezing_cum`, …) and units, in `chion_io.f90` and `input/chion-variables-itm.md`. ITM's `smb` then means total SMB [mm w.e. d-1], unlike BESSI/PDD's ice-facing `smb` — document it. Write correct long_names (Julia labels the rates "Cumulative", 1c.7) and report that upstream. Fix stale dim-order comment `chion_io.f90:50-56`. Keep `dust_dep`/`alb_ice_host` Fortran names (host contract) and document the mapping.
Verify: `tests/test_io.f90`; yelmox build unaffected (names it reads unchanged).

**WP16 — Docs.** `PLAN.md` header (new reference commit), `porting_notes.md` (new D-entries for WP2/5/9–11; D22/D24/D27 updates; upstream-defect status: 19, 20, 11, 24?, 26 closed), `pdd_defects.md`, CHANGELOG, `validation/README.md`. Draft the upstream issue list from §1c for Nils.

### Stage C — chion-originated, coordinated with Nils

**WP17 — Thin-snow albedo in BESSI (PLAN ONLY; design open, Q-T1..T7).**
*Problem* (yelmox GRL-8KM, MAR forcing, `model=bessi`): the host spreads monthly precipitation over every day, and its smooth snow fraction gives trace snowfall even at +2 to +5 °C. BESSI treats any surface layer (`mass(1) > TOL_EMPTY_LAYER`) as full snow cover. `albedo_update` step 2 clamps the stored albedo into `[alpha_wet, alpha_dry]`, so a bare column (`alpha_ice`) jumps to ≥ 0.70 on the first trace of snow (`snow_albedo.f90:151-225`, bare test in `snow_bessi.f90:629-633, 822`). `alpha_ice` is almost never seen: a sweep of 0.2/0.3/0.4 gives identical melt (478 Gt/yr). With a near-step snow fraction (`sf_a=5`), melt is 818 vs 670 Gt/yr. Under ssp585 BESSI gives −0.149 m SLE by 2300 vs ITM's −0.344, and the ablation area does not grow. The upstream `:aging` scheme makes this worse: any snowfall > 0 resets the albedo to `alpha_dry`.
*Precedent:* ITM blends `alb_bg + min(H_snow/H_crit,1)·(alb_snow − alb_bg)` (`snow_itm.f90:calc_albedo_surface`). CLIMBER-X SEMIX blends `f_snow = tanh(h_snow/(c_fsnow·z0m_ice))·f_snow_orog` against `alb_bg` (`climber-x/src/smb/smb_surface_par.f90:108-129`). chion's SEMIX port left that blend out on purpose (`semix_port_scope.md` §"3. snow-cover-fraction blend"). Chion.jl `dev_nils` has nothing equivalent: bare is a hard switch at `EPS_EMPTY_LAYER` in all schemes (`albedo.jl`). Raise with Nils (N5).
*Sketch* (to be fixed by the answers below):
`alpha_eff = f·alpha_snow + (1−f)·alb_ice_use`, with `f = f(column SWE or depth)` and `f(0)=0`, `f→1`. `alpha_snow` stays the prognostic snow albedo (aging/refresh memory, unblended). `alpha_eff` feeds the SEB and is written as `albedo`. The bare-ice branch is unchanged: f→0 makes the switch continuous.
*Files:* `snow_albedo.f90` (blend function, applied after each scheme), `snow_albedo_semix.f90` (if SEMIX uses its own f_snow), `snow_bessi.f90` (store `albedo_snow` and `albedo`; SEB uses `albedo`), `chion_defs.F90` + nml (parameter), `chion_io.f90` (restart field `albedo_snow`; output).
*Verification:*
(i) Unit tests: f(0)=0, monotone, f(deep)=1 to round-off; alpha_eff continuous across the bare transition; blend off gives the old step bit-identical.
(ii) Chion.jl gate unchanged: the blend is off under `legacy_chion` (or the parameter sets it off in the harness) until Chion.jl adopts it.
(iii) yelmox GRL-8KM at the default host `sf_a`: an `alpha_ice` 0.2/0.3/0.4 sweep must give clearly different melt, with a spread comparable to the `sf_a=5` case (818 vs 670 Gt/yr). Check runoff against MAR (645 Gt/yr reference) and the ablation-area growth under ssp585.
*Ordering:* after WP6 (needs `:aging`) and WP2 (bare-ice runoff); independent of Stage B. Do not implement before T1–T7 are answered.

Order rationale: harness first so every physics WP has a gate; Stage A changes are small and mostly bit-id; Stage B's three temperature-physics changes cannot be gated individually against one upstream commit (upstream has no switches), so WP9–11 go green together at WP11.

## 3. yelmox impact

- ITM (current default): no result change from this plan. WP4 adds verification only; WP15 may rename ITM output variables (yelmox does not read chion output files).
- BESSI host (next step): WP13 needs no new host field. yelmox never sets `solar_longitude_deg` (stays 0) — must be fixed before BESSI with diurnal SW or SEMIX `coszm`. Upstream drives SEMIX with 10 m wind; our semix uses z_sfl = 100 m.

## 4. Decisions (user, 2026-10-09)

1. Reference: `dev_nils` tip; report anything on `main` that dev_nils lacks (§2 tracking rule).
2. Adopt upstream defaults (WP14), except `turbulent_flux_scheme` (see 5).
3. Albedo aging × dt: adopt (WP5).
4. Calonne + harmonic: replace outright, no switch (WP9, WP10).
5. Julia's modified SEMIX turbulence: not ported for now (N3).
6. `max_lwc`: keep **0.1** (current value on both sides at dev_nils tip); ask Nils (N1).
7. `refreezing_correction`: wait (WP7 deferred).
8. PDD default `pdd_method = "simple"` (WP3).
9. Gas constant: keep ours (8.31446); `legacy_chion` reproduces upstream's 8.314 until Nils changes it (WP1, N2).
10. Bare-ice rain: fix the double count, route rain once (WP2).
11. ITM D27 under `legacy_chion`: yes (WP4).
12. `snow_age_days` in restart: yes (WP6).
13. ITM output names: Julia's (WP15).

## 5. For Nils (to send in parallel)

- **N1** `max_lwc`: 0.05 (`385d452`) or 0.1 (restored inside `5b4fe49`)? chion keeps 0.1 meanwhile.
- **N2** Gas constant: please use 8.31446261815324 (already `DEFAULT_UNIVERSAL_GAS_CONSTANT` in `constants.jl`) instead of the literal 8.314 in `densification.jl:10,151,174`.
- **N3** `:semix` turbulence (now BESSI's default) is no longer CLIMBER-X SEMIX (§1c.1). Intended as a new scheme?
- **N5** Thin-snow albedo (WP17): BESSI has no partial snow cover, so `alpha_ice` is almost never seen under trace snowfall, and `:aging` resets to `alpha_dry` on any snowfall. Would Chion.jl adopt the same blend (T1–T7), so the two stay comparable?
- **N4** Bugs in §1c: SEMIX albedo (1c.2), ITM D27 (1c.3), bare-column prescribed albedo (1c.4), rain double count (1c.5), NaN diurnal amplitude (1c.6), ITM output long_names (1c.7), stale tests/docs (1c.8), gravity 9.81 (1c.9).

## 6. Open questions — WP17 thin-snow albedo

- **T1 Blend function:** ITM-style linear `min(1, SWE/SWE_crit)`, or SEMIX-style `tanh(h/(c·z0))`? Recommended: one function for all schemes; SEMIX gets its own CLIMBER-X form only if exact SEMIX fidelity matters.
- **T2 Variable:** column SWE `sum(mass+mass_w)` [kg m-2] (density-independent, matches ITM's mm w.e.), surface-layer mass only, or snow depth [m] (SEMIX)? Whole-column SWE treats thin firn over ice as snow; surface-layer mass reacts only to fresh snow.
- **T3 Parameter name/default:** e.g. `swe_crit_albedo` [kg m-2], default 10 (ITM's `H_snow_crit_desert`; ITM goes up to 100 with PDDs)? Default on or off? On changes default BESSI results vs Chion.jl, so `legacy_chion` must turn it off.
- **T4 Schemes:** dynamic and aging yes; constant yes (also hard-switches today); SEMIX: chion's blend or CLIMBER-X's f_snow (incl. `f_snow_orog`, which needs `z_sur_std`); prescribed: no — a prescribed albedo is already the total surface albedo.
- **T5 State:** keep the prognostic snow albedo (`albedo_snow`, aging/refresh memory) separate from the effective `albedo`, which needs a new restart field? Or blend in place, which is simpler but makes aging start from a blended value?
- **T6 `:aging` reset:** separately require a minimum snowfall for the `alpha_dry` reset (e.g. the existing `1−exp(−dm/3)` refresh scale), or rely on the blend alone?
- **T7 Background albedo:** `alb_ice_use` (`alpha_ice` or the host's `alb_ice_host`) everywhere, including tundra/land columns with `H_ice=0`? ITM uses a land/forest background there.
