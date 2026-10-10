# Chion.jl: issues and suggestions from the chion port

For Nils. Collected while syncing chion (the Fortran port) with Chion.jl `main` `9ec6cc7`
(= `03bb445`). File and function references are at `9ec6cc7`. "chion" points to the entry
in chion's `docs/porting_notes.md` (D-numbers) or `docs/PLAN_dev_nils.md`, where the
analysis and measurements are. chion carries each physics fix below as a deviation that its
`legacy_chion` build reverts, so the validation against Chion.jl stays green either way;
adopting a fix upstream lets chion drop the deviation.

Status of earlier reports: fixed upstream since the port base are the vapour diagnostics
(defect 1), the NaN diurnal amplitude, the bare-ice rain leak (defects 11/20, but see 1),
aging × `dt` (defect 19), the CF time axis (24), the output time coordinate (26) and the PDD
budget (#12-#14, #19). Open from the PDD list: #17's single `temperature_sigma` and no
`melt`/`refreezing` in `PDDState`.

## A. Wrong results

**1. Rain on a near-empty column is counted twice.**
When `0 < mass(1) <= EPS_EMPTY_LAYER`, accumulation has already put the step's rain into
`mass_w(1)`; the bare-ice branch then adds `rainfall_mass` to `runoff` again (since
`8fff530`). The substrate bare path does the same (`liquid_mass = rain + melt`).
*Where:* `src/step.jl:column_step_core!` (bare branch, `runoff += rainfall_mass + melt_mass`,
l. 465); `src/step.jl:_step_bare_ice_substrate!` (l. 354).
*Fix:* route rain once, where it falls: to `runoff` in accumulation when there is no layer
to hold it (`n = 0` or `mass(1) <= 0`), into `mass_w(1)` otherwise; add no rain in either
bare path. *chion:* D29.

**2. ITM firn warming is ~360× too small.**
`Tsrf = min(T0, Tair + firn_fac*max(melt_net_rate, 0))` uses the daily rate [mm w.e. d-1],
but `firn_fac` is smbpal's calibration against the annual net melt [mm w.e. yr-1, 360-day
year].
*Where:* `src/processes/itm.jl:_itm_apply_column!` (l. 104).
*Fix:* `firn_fac*max(melt_net_rate, 0)*360` (smbpal's year), and document `firn_fac` as
`[K (mm w.e. yr-1)-1]`. *chion:* D27.

**3. Fine near-surface layers leave the layer below them unbounded.**
With `near_surface_layer_max_thicknesses_m` set, the surface split only fires when
`mass(1) > mass_max` right after snowfall, and the surface merge is off. The remesh pushes
every increment down the four fine layers into layer 5, which nothing splits or merges: it
grows until the 22.5 m depth cap exports its base. A firn column becomes 4 thin layers over
one cell of up to ~22 m, so densification, cold content, percolation and refreezing below
0.47 m happen in one cell. The docs ("Without these limits the surface layer can hold up to
`mass_max`") suggest this is not intended. 10-yr percolation-zone column, chion's fix vs
Chion.jl's behaviour: mean layer count 4.9 -> 10.0, runoff +7.7 %, refreezing -6.0 %, liquid
water -40 %.
*Where:* `src/processes/layer_structure.jl:_remesh_near_surface_layers!` (l. 553),
`_split_surface_layer!`, `_merge_surface_layer!`.
*Fix:* after the remesh, apply the surface split/merge by mass to the first unlimited layer
(`k0` = number of limited layers + 1): split while `mass(k0) > mass_max`, merge while
`mass(k0) < mass_min` (generalise the two functions and `_free_slot_for_surface_split!` to
a layer index). *chion:* D32.

**4. Column reset leaves the ice substrate temperature stale.**
`_reset_bessi_columns_kernel!` resets the snow state but not `ice_temperature`, so a
re-activated column starts on the substrate of its previous life.
*Where:* `src/runtime.jl:_reset_bessi_columns_kernel!` (l. 135).
*Fix:* set `ice_temperature` to `temperature_init`, as `_initialize_bessi_state_kernel!`
does. *chion:* D34.

**5. Unsplit diurnal days lose their shortwave in the polar night.**
With `diurnal_shortwave_substeps` on, a day the criterion does not split (`n_substeps = 1`:
cold, no sun, dt outside 0.75-1.25 d) still goes through the interval average over
`[-π, π]`, which rebuilds the shortwave from the solar geometry. It returns the forcing's
value where the sun rises but zero where the geometry has no daylight (polar night, missing
latitude), whatever shortwave the forcing carries. Affects synthetic or smoothed forcing and
calendar mismatches at the polar-night edge.
*Where:* `src/step.jl:_column_step_diurnal!` (l. 260).
*Fix:* step an unsplit day with its forcing unchanged. *chion:* D39.

**6. SEMIX snow albedo differs from CLIMBER-X.**
(a) Dang direct visible dust term uses `-0.05 p²`; CLIMBER-X `-0.049`. (b) WW 1-10 ppm band
drops the `0.01` intercept and uses `log(d)` from 1 instead of 1.001 (CLIMBER-X table
`(1.001, 10, 100, 1000)` -> `(0.01, 0.05, 0.15, 0.30)` aged). (c) Dang's
`dalb_snow_vis/nir` offsets are never applied. (d) The seasonal SWE maximum resets on
`day_of_year <= dt_days`: never fires with mid-day stamps, and 1 January is mid-summer in
Antarctica. (e) The old bare-ice path ignores the host ice albedo (the substrate path uses
`_bare_ice_albedo`).
*Where:* `src/processes/albedo.jl:_semix_ww_bands` (l. 27), `_semix_dang_bands` (l. 40),
`_update_semix_surface_albedo!` (l. 80); `src/step.jl:column_step_core!` bare branch.
*Fix:* follow `climber-x/src/smb/smb_surface_par.f90`; reset the seasonal maximum on a
hemisphere-aware date. *chion:* `src/physics/snow_albedo_semix.f90` matches CLIMBER-X.

**7. Bare column with prescribed albedo keeps a stale albedo** on the old bare-ice path
(no substrate): `use_prescribed_albedo || set alpha_ice` leaves the previous value.
*Where:* `src/step.jl:column_step_core!` (l. 458). *Fix:* apply the prescribed albedo there
too (the substrate path does via `_bare_ice_albedo`). *chion:* `snow_bessi.f90` applies it.

**8. HTESSEL thermal metamorphism is dead above ~150 kg m-3.**
`exp(-4.2e-2 (T0-T) - 460 max(0, ρ-150))`: the 460 multiplies a density excess in kg m-3,
so the term underflows to zero for ρ > 150.15. Probably a missing unit conversion (460 per
Mg m-3?).
*Where:* `src/processes/densification.jl:_htessel_thermal_metamorphism` (l. 19-24).
*Fix:* check the source (Vionnet et al. / HTESSEL); `ρ` in g cm-3 or `0.46` per kg m-3.
*chion:* porting_notes open item 3 (ported as is).

## B. Physics choices to consider

**9. Thin-snow albedo (snow-cover fraction).**
BESSI treats any surface layer as full snow cover: a trace of snow on bare ice lifts the
albedo from `alpha_ice` to at least `alpha_wet`. With a host or forcing that snows a trace on
most melt days (monthly means interpolated to days: snowfall on 98 % of ablation-zone melt
days, 83 % below 1 mm d-1), `alpha_ice` is almost never seen: GRL-16KM melt 335/334/333
Gt/yr for `alpha_ice` 0.3/0.4/0.5. With fine layers `mass(1)` is ~6 kg m-2 on every column,
so a fraction must use the column SWE.
*Where:* `src/processes/albedo.jl`, `src/step.jl` (bare switch at `EPS_EMPTY_LAYER`).
*Suggestion:* `albedo = f*alpha_snow + (1-f)*alpha_bg`, `f = min(1, SWE/swe_crit)` (10 kg m-2;
CLIMBER-X's `tanh` form for the SEMIX albedo), with the snow albedo a separate state. With
it, `alpha_ice` acts again (melt 562/500/440) and `alpha_ice = 0.50` matches MAR (SMB 340
vs 348 Gt/yr, R² 0.86); `0.40` was calibrated without a blend. *chion:* D40.

**10. Land columns.** Chion.jl has no ice thickness, so it puts the ice substrate under
tundra and melts bare ice under every column, crediting `smb_ice`.
*Suggestion:* an optional `ice_thickness` (it exists for ITM): where it is 0, no substrate,
no ice ablation, a land background albedo (0.2), snow-free surface at the air temperature.
*chion:* D34, D41.

**11. `:aging` resets on any snowfall.** A trace of snow fully rejuvenates the surface; with
daily-interpolated monthly precipitation the ablation-zone albedo never leaves `alpha_dry`.
*Where:* `src/processes/albedo.jl:_update_aging_surface_albedo_arrays!` (l. 205).
*Suggestion:* scale the aging progress `E = -ln((α-α_wet)/(α_dry-α_wet))` and the age by
`exp(-S/S_ref)` (`S_ref` 10 kg m-2). The exponential composes exactly across diurnal
substeps (a linear `min(1, S/S_ref)` gives 66 % of a full refresh over 8 substeps).
*chion:* D30. Related: the dynamic scheme's `1 - exp(-dm/3)` brightening is not
substep-invariant (`_refresh_dynamic_albedo_from_snowfall!`, l. 172; porting_notes defect 21),
now default-relevant with 8 substeps.

**12. Bare-ice latent heat in the semix turbulence.** The energy solve uses
`_surface_vapor_latent_heat(T0) = Lv` on bare ice at T0, while the mass is converted with
`Lv + Lm`; each kg sublimated from ice costs the energy budget `Lv` but the mass budget
`Lv + Lm`. Same for BESSI turbulence over the substrate.
*Where:* `src/processes/energy_flux.jl:_surface_vapor_latent_heat` (l. 67) and the semix
latent flux (l. 159); `src/step.jl:_step_bare_ice_substrate!` (l. 346).
*Fix:* `Lv + Lm` for any solid surface (bare ice, substrate), at any temperature.
*chion:* D35 (semix only).

**13. Cloud-proxy TOA is a fixed modern orbit.** `_daily_toa_shortwave` computes TOA from
S0 1361, obliquity 23.44° and a 365-day eccentricity term. A host whose shortwave is a
transmissivity times its own (paleo, 360-day) insolation gets a biased `SWD/TOA` ratio.
*Where:* `src/processes/surface_fluxes.jl:_daily_toa_shortwave` (l. 11), `_cloud_proxy_emissivity`.
*Suggestion:* an optional `toa_shortwave` forcing field that replaces it. *chion:* D33.

**14. `max_lwc` is a hard-coded keyword default (0.1).** It went 0.1 -> 0.05 (`385d452`) and
back to 0.1 inside `5b4fe49`. Which is intended? chion keeps 0.1.
*Where:* `src/processes/percolation.jl:_go_percolation!` (l. 21, 35).
*Suggestion:* a named parameter in `SnowpackPhysicalConstants`.

## C. Constants

**15. Gas constant.** Densification uses the literal `8.314`;
`DEFAULT_UNIVERSAL_GAS_CONSTANT = 8.31446261815324` exists (`src/constants.jl:23`).
*Where:* `src/processes/densification.jl` l. 10, 151, 174. *Fix:* use the constant (3e-4 in
the Arrhenius term). *chion:* D22.

**16. Gravity.** Densification uses `9.81` (l. 225); everything else 9.80665.
*Fix:* standard gravity. *chion:* D25.

**17. Dry-air gas constant.** `_semix_air_density` divides by the literal `287.05`.
*Where:* `src/processes/energy_flux.jl` l. 141. *Fix:* a named `R_dry` (287.058, as
CLIMBER-X). *chion:* D38.

## D. Numerics and robustness

**18. Mass-weighted mean of equal values is not exact.** `(m1 x1 + m2 x2)/(m1 + m2)` can
merge two layers at T0 to `T0 - 1 ulp`; the `:aging` timescale test (`T >= T0`) then depends
on round-off (one harness column diverged by 926 ulp). The remesh now calls it twice per
substep.
*Where:* `src/processes/layer_structure.jl:_mass_weighted_mean` (l. 17).
*Fix:* `x1 + m2/(m1+m2)*(x2 - x1)`. *chion:* D31.

**19. Bottom deplete sets `Tsrf = T(1)`.** Every other path treats `Tsrf` as the Robin
interface temperature; after the depth cap the next step linearizes at the top cell's
centre.
*Where:* `src/processes/layer_structure.jl:_continuous_bottom_deplete!` (l. 361).
*Fix:* leave `Tsrf` (or `T0` when the column empties). *chion:* porting_notes defect 27
(ported as is).

**20. `_merge_surface_layer!` divides by `subsurface_mass` unguarded** (l. 128): a namelist
with `mass_min > 2*mass_split` divides by zero. *chion:* defect 8.

**21. Fine-layer fill-up loop re-tests the deficit against `EPS_TINY`** after a transfer
that already filled the receiver; harmless in Float64, a round-off loop in Float32.
*Where:* `_fill_near_surface_layer_thicknesses!` (l. 513-519). *Fix:* stop after a transfer
that leaves the donor non-empty. *chion:* D36.

**22. Column reset leaves the four diagnostics stale** (`thickness`, `wet_mass`,
`bulk_density`, `liquid_water`), unlike initialisation (`src/runtime.jl` l. 135).
*chion:* defect 23.

## E. I/O and API

**23. ITM rates labelled cumulative.** `ITM_OUTPUT_VARS` `melt`, `runoff`, `refreezing` are
the step's rates [mmWE day-1] but take BESSI's metadata ("Cumulative melt", mmWE).
*Where:* `src/io.jl` `NETCDF_METADATA` (l. 17-19). *Fix:* ITM entries ("ITM melt rate",
"mmWE day-1"). *chion:* D42. Note also that ITM's `smb` is the total SMB rate while BESSI's
`smb_ice` is ice-facing; worth stating in the long names.

**24. `H_snow` cannot be selected by name.** `normalize_netcdf_variables` lowercases every
token (`:h_snow`), which is not in `ALL_OUTPUT_VARS`, so `"H_snow"` errors; only `"all"`
writes it (`N` and `Tsrf` have special cases).
*Where:* `src/io.jl:normalize_netcdf_variables` (l. 187). *Fix:* match case-insensitively
against the field names.

**25. `load_forcing_file` reads no ITM inputs.** `ice_thickness` and `annual_pdd` cannot be
given in a forcing file; chion's harness merges them in by hand.
*Where:* `src/dataloaders.jl:load_forcing_file` (l. 234). *Fix:* optional static fields, as
`surface_height`.

**26. Output is Float32.** The buffers are `Matrix{Float32}`/`Array{Float32,3}`, which sets
a ~1.2e-7 relative floor on any comparison or restart through the files.
*Where:* `src/io.jl` (l. 76-79). *Suggestion:* a Float64 option. *chion:* defect 25.

## F. Tests and docs

**27.** `docs/src/processes/percolation.md` ("Pore-Collapse Routing") says liquid water in
a layer with zero solid mass is routed onward like a collapsed pore; the code sends it
straight to runoff (`src/processes/percolation.jl:_go_percolation!`, l. 44-47). chion
defect 15.
