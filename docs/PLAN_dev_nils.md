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
6. `d0146e1`: `diurnal_temperature_cycle=true` with missing `surface_height` (NaN) gives NaN amplitude even with γ=0. **Fixed upstream in `03bb445`; chion WP13 matches (non-finite height = no excess).**
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

**WP6b — Continuous snowfall rejuvenation for `:aging` (chion deviation; T6 decided, T10–T12 open).** Upstream resets to `alpha_dry` (age 0) on any snowfall > 0, so a trace of snow fully rejuvenates the surface. Replace that with a partial refresh, in proportion to the step's fresh snowfall `S` [kg m-2] over a reference amount `S_ref`:
- refresh fraction `f_new = min(1, S/S_ref)`;
- aging progress `E = −ln((α − α_wet)/(α_dry − α_wet))` (dimensionless; `E = Σ dt/τ`, so it is exact for `:aging` even though τ switches between 20 d and 2 d) is reduced by that fraction: `E ← (1 − f_new)·E`;
- equivalently, in albedo space and with no new state, `α ← α_wet + (α_dry − α_wet)·((α − α_wet)/(α_dry − α_wet))^(1 − f_new)`;
- the diagnostic `snow_age_days ← (1 − f_new)·snow_age_days`.

This gives `S = 0`: no change; `S ≥ S_ref`: the upstream reset; a trace: almost none. The refresh acts on the prognostic snow albedo (`albedo_snow`, WP17), before the step's aging relaxation, i.e. in the existing `albedo_refresh_from_snowfall` slot (`snow_accumulation.f90`). New parameter `aging_snowfall_ref` [kg m-2] in `&bessi`/`chion_const`; proposed default 10 kg m-2 (≈ 3 cm of fresh snow, the e-folding depth of Oerlemans & Knap 1998) — T10. `legacy_chion` reverts to upstream's binary reset, so the gate keeps covering `:aging`. New D-entry; raise with Nils (N5).
Verify: unit tests — `S = 0` leaves α unchanged; `S ≥ S_ref` gives `alpha_dry` and age 0; successive refreshes compose as `E·(1−f₁)(1−f₂)`; monotone in `S`. Gate with `albedo=:aging` under legacy unchanged. In yelmox, the albedo no longer recovers to `alpha_dry` on trace-snowfall days in the ablation zone.

**WP7 — `refreezing_correction` + clamp (`27113b6`) — DEFERRED (decided: wait).** Neutral at Chion.jl's default (1.0) up to an ulp near T0, so the Stage B gate is unaffected. Revisit if upstream changes the default.

### Stage B — reference = `origin/dev_nils` (`27113b6`)

**WP8 — Re-point harness to `27113b6`.** Same as WP0. Expected failures: all temperature-driven fields (Calonne, harmonic, latent heat).

**WP9 — Harmonic interface conductance + `β = −dt/ci`.** Together, one commit. `snow_energy.f90:113-136, 418`.
Verify: unit test — uniform column bit-id to old form up to round-off; analytic two-layer steady flux test in `tests/test_energy.f90`. Gate still red until WP10–11 (report diff reduction).

**WP10 — Calonne (2019) conductivity.** `snow_energy.f90:97-111`, layer-1 at T^n, others at their own T, `rho_i` from constants. Decided: replace outright, no switch (upstream has none; a switch would carry a dead `Ki` path that nothing validates). Remove `Ki` from `chion_defs`, `chion_api.f90:1048`, nml. Same for WP9: harmonic replaces arithmetic.
Verify: unit test against hand-computed values at (ρ,T) = (300,253), (500,253), (917,263); gate expected green only after WP11 (report remaining diff). Report 10-yr column change.

**WP11 — Phase-dependent latent heat; gradient-based vapour mass (BESSI turbulence).** `snow_vapor.f90`, `snow_surface_fluxes.f90:477-483`; bare-ice keeps L_v+L_m (as upstream).
Verify: gate fully green at `27113b6` (with explicit `albedo=:dynamic`, `alpha_wet=0.70`, `turbulent_flux_scheme=:bessi`); BESSI closure tests pass (re-check defect-1 vapour closure — note the harness runs humidity-off, so add a humidity-on unit test).

**WP12 — `turbulent_flux_scheme` flag split from `seb_scheme`. SUPERSEDED by C6/C7.** Flags + parser; `seb_scheme` keeps longwave/`eps_ice`; `turbulent_flux_scheme = "bessi" | "semix"` where `semix` = our CLIMBER-X `snow_seb_semix.f90` (not Julia's, 1c.1).
Verify: all four combinations run; `seb_scheme=semix` (both) bit-id to current `semix`; `(bessi,bessi)` bit-id to current default.

**WP13 — Elevation-dependent diurnal T amplitude.** Add `surface_height` to `chion_step_forcing_class` and pack it in `chion_model.f90:117`; params γ [K km⁻¹], z_ref, A_max in `bessi_par`/nml. Guard the missing-height case properly (no NaN, 1c.6). No new host field: yelmox already fills `forc%surface_height`.
Verify: default bit-id; unit test of the clamp; diurnal stays out of the gate (harness note).

**WP14 — Adopt upstream defaults. VOID: upstream 03bb445 reverted to dynamic albedo / `alpha_wet` 0.70; replaced by C11.** `albedo = "aging"`, `alpha_wet = 0.60`. **Exception:** `turbulent_flux_scheme` stays `"bessi"` — upstream's default `:semix` is Julia's modified scheme, which we are not porting yet (§1c.1, N3). Harness drops the albedo overrides but keeps `turbulent_flux_scheme=:bessi`. Also re-check test/par namelists that assumed the old defaults.
Verify: gate green with Chion.jl defaults; record default-run deltas (10-yr column, ANT/GRL chion_grid) in CHANGELOG.

**WP15 — Output/naming alignment.** Decided: ITM output takes Julia's names (`alb_s`, `smb`, `smbi`, `melt_net`, `smb_cum`, `melt_cum`, `runoff_cum`, `refreezing_cum`, …) and units, in `chion_io.f90` and `input/chion-variables-itm.md`. ITM's `smb` then means total SMB [mm w.e. d-1], unlike BESSI/PDD's ice-facing `smb` — document it. Write correct long_names (Julia labels the rates "Cumulative", 1c.7) and report that upstream. Fix stale dim-order comment `chion_io.f90:50-56`. Keep `dust_dep`/`alb_ice_host` Fortran names (host contract) and document the mapping.
Verify: `tests/test_io.f90`; yelmox build unaffected (names it reads unchanged).

**WP16 — Docs.** `PLAN.md` header (new reference commit), `porting_notes.md` (new D-entries for WP2/5/9–11; D22/D24/D27 updates; upstream-defect status: 19, 20, 11, 24?, 26 closed), `pdd_defects.md`, CHANGELOG, `validation/README.md`. Draft the upstream issue list from §1c for Nils.

### Stage C — reference = Chion.jl `main` `9ec6cc7` (= `03bb445`, "calibrated GrIS surface setup")

Added 2026-10-09 after Nils merged `03bb445` to main (main == dev_nils). Full review:
side-session-notes `chionjl-03bb445-review.md`. Stage B (WP8–WP11, WP13) finishes against `27113b6` first; WP12 is folded into C6/C7, WP14 into C11, WP17 into C12. Levante worktree: `/work/ba1442/robinson/Chion.jl-main`.

- **C0** Harness against `9ec6cc7`, with pins reproducing chion's physics where options still exist (graybody LW, bessi SEB/turbulence, `ice_substrate_layers=0`, near-surface thicknesses off, diurnal off, albedo values). Expected red: the old surface formulation is gone upstream.
- **C1** Non-switchable small fixes: fresh snow on bare ice takes air T in every new layer; vapour-mass clipping closed (our defect 1); wet-albedo relaxation `(1−r)^dt`; depth cap constant 22.5 m; NaN guards; drop namelist aliases (canonical `warren_wiscombe`).
- **C2** Robin surface boundary: Tsrf = interface temperature via half-cell conductance `2K1/dz1`, no surface heat capacity, melt energy `Q(T0) − Gs(T0−T1)`, vapour at Tsrf. C1+C2 gate green together. **Done (C0–C2): gate vs `9ec6cc7` green, worst 0.47 ulp.**
- **C3** Ice substrate (`ice_substrate_layers`, default 5; 0.05 m doubling to 1.55 m; insulated base; one matrix with snow); bare ice carries cold content. chion: reset `ice_temperature` on column reset (Q7, deviation, N8); none on land columns `H_ice = 0` (Q6); old restarts init `min(t_srf, T0)` (Q8). Gated as its own configuration. **Done:** default 0 (C11 switches), D34; substrate configuration green (0.47 ulp) with a dry `bare_ice` column for cold-content coverage; default column bit-identical. yelmox must fill `forc%H_ice` for BESSI (only done for ITM).
- **C4** Fine near-surface layers (0.02/0.05/0.10/0.30 m, cap-down/fill-up remesh twice per step, 100 kg m-2 surface merge disabled), ported faithfully, gated. Remesh uses D31's exact mean (Q16).
- **C4b** chion deviation: split/merge resumes on layer 5 (upstream accumulates everything below the fine layers in one layer up to the 22.5 m cap); reverted under `legacy_chion`; N7.
- **C5** Cloud-proxy longwave `ε = 0.624 + 0.0032(Ta−T0) + 0.613n`, `n = 1 − SWD/(TOA·(0.85+0.075 z/km))`, internal TOA from lat/solar longitude/day; optional host `toa_shortwave` (Q5, N9); needs `surface_height` in step forcing (WP13).
- **C6** `seb_scheme` = longwave only (`bessi` graybody | `semix`), turbulence separate (supersedes WP12).
- **C7** Port Julia's bulk turbulence as `turbulent_flux_scheme = "semix"` (sensible factor, stable coefficient); chion's CLIMBER-X SEMIX surface scheme renamed `"climberx"` (Q2, reverses §4.5 / 1c.1 because the calibrated upstream defaults depend on it). chion fixes: bare-ice SEMIX latent uses Lv+Lm (Q12) and `R_dry` = 287.058 (Q13), both reverted under `legacy_chion`.
- **C8** First gated diurnal configuration (8 substeps, 1 K cycle).
- **C9** Prescribed melt forcing — deferred (Q9). `surface_smb` monthly output — skipped (Q10).
- **C10** Output additions folded into WP15.
- **C11** Adopt the full calibrated upstream default set (`alpha_ice` 0.40, `alpha_wet` 0.70, dynamic albedo, `seb=semix`, Julia SEMIX turbulence 2.5/40, 8 substeps with 1 K cycle, substrate 5, fine layers with C4b); harness then gates with Chion.jl defaults (+ legacy for chion deviations).
- **C12** WP17 rework: snow cover from column SWE, not `mass(1)` (Q15 revises T2: fine layers pin `mass(1)` ~6 kg m-2); WP6b refresh shape `1 − exp(−S/S_ref)` (Q14 revises T11: linear is not substep-invariant, 66 % of a full refresh over 8 substeps).
- **C13** Docs, CHANGELOG, list for Nils (N7–N10).

Then: WP15 (outputs), WP18 (build/performance), WP19 (MAR/RACMO evaluation incl. `aging_snowfall_ref`), WP16 (docs), PR.

Decisions Q1–Q16 (review recommendations, adopted autonomously 2026-10-09 under the user's instruction to take main as the reference and continue): all as recommended — see the review file.

### Stage D — chion-originated, coordinated with Nils (reworked as C12)

**WP17 — Thin-snow albedo in BESSI (PLAN ONLY; T1–T7 decided, T8–T9 open).**
*Problem* (yelmox GRL-8KM, MAR forcing, `model=bessi`): the host spreads monthly precipitation over every day, and its smooth snow fraction gives trace snowfall even at +2 to +5 °C. BESSI treats any surface layer (`mass(1) > TOL_EMPTY_LAYER`) as full snow cover. `albedo_update` step 2 clamps the stored albedo into `[alpha_wet, alpha_dry]`, so a bare column (`alpha_ice`) jumps to ≥ 0.70 on the first trace of snow (`snow_albedo.f90:151-225`, bare test in `snow_bessi.f90:629-633, 822`). `alpha_ice` is almost never seen: a sweep of 0.2/0.3/0.4 gives identical melt (478 Gt/yr). With a near-step snow fraction (`sf_a=5`), melt is 818 vs 670 Gt/yr. Under ssp585 BESSI gives −0.149 m SLE by 2300 vs ITM's −0.344, and the ablation area does not grow. The upstream `:aging` scheme makes this worse: any snowfall > 0 resets the albedo to `alpha_dry`.
*Precedent:* ITM blends `alb_bg + min(H_snow/H_crit,1)·(alb_snow − alb_bg)` (`snow_itm.f90:calc_albedo_surface`). CLIMBER-X SEMIX blends `f_snow = tanh(h_snow/(c_fsnow·z0m_ice))·f_snow_orog` against `alb_bg` (`climber-x/src/smb/smb_surface_par.f90:108-129`). chion's SEMIX port left that blend out on purpose (`semix_port_scope.md` §"3. snow-cover-fraction blend"). Chion.jl `dev_nils` has nothing equivalent: bare is a hard switch at `EPS_EMPTY_LAYER` in all schemes (`albedo.jl`). Raise with Nils (N5).
*Design (decided):*
- **Blend** (T1, T2): `f = min(1, mass(1)/swe_crit_albedo)` on the surface layer's solid mass; `alpha_eff = f·alpha_snow + (1 − f)·alpha_bg`. With `n > 1`, merging keeps `mass(1) ≥ mass_min` (100 kg m-2), so firn columns get `f = 1` by construction. The blend acts where it should: seasonal snow over ice, and fresh snow on a bare column (`n = 1`, small `mass(1)`).
- **Parameter** (T3): `swe_crit_albedo` [kg m-2], default 10, **on by default**. `legacy_chion` sets it off (f ≡ 1), so the Chion.jl gate is unchanged.
- **Schemes** (T4): dynamic, aging and constant use the blend. SEMIX uses CLIMBER-X's own `f_snow = tanh(h_snow/(c_fsnow·z0m_ice))·f_snow_orog`, with `f_snow_orog = h_snow/(h_snow + c_fsnow_orog·z_sur_std)` when `z_sur_std` is supplied (`climber-x/src/smb/smb_surface_par.f90:108-129`); this closes `semix_port_scope.md` item 3. Prescribed albedo is untouched.
- **State** (T5): `albedo_snow` is prognostic (aging, refresh and WP6b memory, unblended) and gets a restart field. `albedo` is the effective value: it feeds the SEB and is written as output. Restarts without `albedo_snow` initialise it from `albedo`.
- **Background** (T7): `alpha_bg = alb_ice_use` (`alpha_ice` or the host's `alb_ice_host`) where `H_ice > 0`; a land albedo where `H_ice = 0` (T8). Needs `H_ice` in `chion_step_forcing_class`. The host already fills `forc%H_ice` for ITM, so this is packing only, not a new host field.
- **Bare-ice branch**: unchanged; `f → 0` as `mass(1) → TOL_EMPTY_LAYER` makes the switch continuous. The bare branch uses `alpha_bg` too (T9).
- **`:aging` snowfall** (T6): continuous rejuvenation, see WP6b.

*Files:* `snow_albedo.f90` (blend, background), `snow_albedo_semix.f90` (SEMIX f_snow and orography), `snow_bessi.f90` (`albedo_snow`/`albedo`, bare branch), `chion_defs.F90` + nml (`swe_crit_albedo`, `alpha_land`, SEMIX `c_fsnow`, `c_fsnow_orog`; `H_ice` in step forcing), `chion_model.f90` (pack `H_ice`), `chion_io.f90` (restart `albedo_snow`; output).
*Verification:*
(i) Unit tests: f(0) = 0, linear, f = 1 above `swe_crit_albedo`; `alpha_eff` continuous across the bare transition; with f ≡ 1 (legacy) the old step is bit-identical; SEMIX f_snow matches CLIMBER-X's formula; restart round-trip of `albedo_snow`.
(ii) Chion.jl gate unchanged under `legacy_chion` (blend off; WP6b reverted).
(iii) yelmox GRL-8KM at the default host `sf_a`: an `alpha_ice` 0.2/0.3/0.4 sweep gives clearly different melt (target: a spread comparable to the `sf_a=5` case, 818 vs 670 Gt/yr). Check runoff against MAR (645 Gt/yr) and ablation-area growth under ssp585.
*Ordering:* after WP6/WP6b and WP2; independent of Stage B. Implement after T8–T9.

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
- **N7** Fine layers: everything below the 4 fine layers accumulates in one layer up to the 22.5 m cap. **N8** substrate temperature not reset on column reset. **N9** internal fixed-orbit TOA vs host insolation (paleo). **N10** substep-invariant snowfall refresh.
- **N5** Thin-snow albedo (WP17): BESSI has no partial snow cover, so `alpha_ice` is almost never seen under trace snowfall, and `:aging` resets to `alpha_dry` on any snowfall. Would Chion.jl adopt the same blend (T1–T7), so the two stay comparable?
- **N4** Bugs in §1c: SEMIX albedo (1c.2), ITM D27 (1c.3), bare-column prescribed albedo (1c.4), rain double count (1c.5), NaN diurnal amplitude (1c.6), ITM output long_names (1c.7), stale tests/docs (1c.8), gravity 9.81 (1c.9).

## 6. WP17 / WP6b decisions and open questions

Decided (user, via the yelmox session, 2026-10-09): T1 linear `min(1, SWE/SWE_crit)`; T2 surface-layer mass; T3 `swe_crit_albedo` = 10 kg m-2, on by default; T4 dynamic/aging/constant blend, SEMIX uses CLIMBER-X f_snow with orography, prescribed untouched; T5 separate `albedo_snow` with a restart field; T6 continuous snowfall rejuvenation (WP6b); T7 land background where `H_ice = 0`.

Decided 2026-10-09 (user accepted the proposals): T8 constant `alpha_land` = 0.2; T9 skip ice ablation on bare land columns (`H_ice = 0`), inside WP17; T10 `aging_snowfall_ref` = 10 kg m-2, separate parameter; T11 linear; T12 dynamic refresh unchanged (Chion.jl parity).
Levante test failures under ifx `-Ofast` (test_layers miscompiled, sp 1-ulp in test_energy/test_itm): option (a) — tests and validation builds use `-O2 -fp-model precise`; production stays `-Ofast`.

Former open questions (kept for the record):
- **T8 Land background:** a constant `alpha_land` (e.g. 0.2), or ITM's PDD-dependent land/forest blend `alb_land·(1000 − PDD)/1000 + alb_forest·PDD/1000`? The latter needs annual PDDs in the BESSI step forcing (yelmox already computes them for ITM).
- **T9 Bare land columns:** with `H_ice = 0` the bare branch still computes ice melt and credits `smb_ice`. Skip ice ablation on land columns, i.e. no ice to melt, the surface energy goes nowhere? Separate fix, or part of WP17?
- **T10 `aging_snowfall_ref` default:** 10 kg m-2 (≈ 3 cm of fresh snow; Oerlemans & Knap 1998 e-folding depth)? Same value as `swe_crit_albedo`, or tied to it (one parameter)?
- **T11 Refresh shape:** linear `min(1, S/S_ref)` (proposed; it saturates at a finite snowfall) or `1 − exp(−S/S_ref)`, like the dynamic scheme's existing refresh (`ALBEDO_SNOWFALL_EFOLD_MASS` = 3 kg m-2)?
- **T12 Dynamic scheme:** it already refreshes continuously (`1 − exp(−dm/3)`, linear in albedo). Leave it, or switch it to the same age-space form so both schemes share one refresh?
