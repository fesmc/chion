# PDD defect list — `Chion.jl/src/processes/pdd.jl`

WP9 deliverable. `docs/PLAN.md` §3.2 states that Chion.jl's `PDDModel` is known not to be
fully working, so this file is the primary intellectual output of WP9: what is wrong, what
`smbpal/src/smb_pdd.f90` does instead, what the physical consequence is, and what should be
fixed upstream.

Sources compared, line numbers as of `Chion.jl@main`:

| role | file |
|---|---|
| structure (ported) | `Chion.jl/src/processes/pdd.jl`, `src/models.jl:80-116`, `src/state.jl:180-195`, `src/runtime.jl:250-305` |
| physics (authoritative where they disagree) | `smbpal/src/smb_pdd.f90` — `calc_ablation_pdd`, `calc_temp_effective`; call sites `smbpal/src/smbpal.f90:290-310`, `:495-505` |
| port | `chion/src/physics/snow_pdd.f90`, `chion/tests/test_pdd.f90` |

Severity: **A** = wrong physics or wrong published quantity; **B** = fragile / latent bug;
**C** = hygiene.

Summary:

| id | severity | one line |
|---|---|---|
| D1 | A | refreezing has no cold-content and no capacity limit |
| D2 | A | `smb_ice` is credited with the snowpack, so it is not "net mass forcing to the ice sheet" |
| D3 | A | `snowpack_swe` is uncapped, un-aged, un-densified and never becomes ice |
| D4 | A | PDD flavour is selected implicitly by timestep length |
| D5 | B | two entry points apply different physics |
| D6 | A | refrozen mass re-enters the melt-able snow reservoir; the reservoir is never exhausted |
| D7 | B | the six-line core exists in three copies that have already diverged |
| D8 | B | the active-column mask is ignored, and the GPU call sites are `MethodError`s |
| D9 | C | no `melt` / `refreezing` diagnostics, unlike `BESSIState` |
| D10 | C | `273.15` and `86400.0` hard-coded rather than taken from the constants struct |
| D11 | C | one global `temperature_sigma`; smbpal uses three, by surface type |
| D12 | C | Abramowitz–Stegun `_normal_cdf` polynomial instead of `erfc` |

One defect was also found **in smbpal**, running the other way — see S1 at the end.

---

## D1 — `refreezing_fraction` has no cold-content and no capacity limit

**Where.** `pdd.jl:52`, `:232`, `:274`:

```julia
refrozen = refreezing_fraction * snow_melt
```

**(a) Is it wrong, or defensible?** Wrong as written, and not merely a simplification.
A flat fraction of *all* snow melt refreezes, unconditionally: regardless of the pack
temperature, regardless of how much pore space exists, and regardless of how much snow fell.
There is no state in `PDDState` that could limit it — no temperature, no density, no
cold content — so the model cannot express "this pack is isothermal at 0 °C and can refreeze
nothing".

A degree-day model does have to parameterise refreezing crudely, so a constant fraction is a
legitimate *form*. What is not defensible is that the fraction multiplies the melt of an
arbitrarily old, arbitrarily large reservoir (D3), with nothing bounding the total.

**(b) What smbpal does.** `smb_pdd.f90:37,47`:

```fortran
refrz_max = acc*f_refrz_max          ! capacity, proportional to ACCUMULATION
refrz     = min(melt_snow, refrz_max)
```

The cap is a fraction of the accumulation, not of the melt. That is the physically meaningful
statement: the refreezing capacity of a snowpack scales with the mass of cold snow that was
deposited (its pore space and its cold content), not with how much of it happened to melt.
`f_refrz_max` is passed as `par%itm%Pmaxfrac` (`smbpal.f90:500`). smbpal takes `acc` to be the
*annual* snowfall, consistent with its annual-diagnostic framing.

**(c) Physical consequence.** Sign: **SMB biased positive, runoff biased low**, and the bias
grows with the size of the snow reservoir.

- *Temperate / isothermal pack* (maritime margin, late summer, percolation zone at melting
  point): true refreezing tends to zero once the pack is isothermal and saturated. Chion.jl
  keeps refreezing 60 % of melt indefinitely. Error is the full `0.6 × snow_melt`; for a
  Greenland percolation-zone summer with 400 mm w.e. of snow melt that is ~240 mm w.e. of
  spurious superimposed ice per year.
- *Multi-year pack* (compounded with D3): `snow_melt` is drawn from a reservoir that includes
  decades of old snow, so `0.6 × snow_melt` can exceed the year's snowfall. smbpal's cap makes
  that impossible by construction.
- The two schemes can also disagree the *other* way in a single cold high-accumulation year
  (smbpal caps at `0.6 × acc` = 300 mm where Chion.jl gives `0.6 × melt` = 240 mm), so this is
  not a uniform offset that could be tuned away with a single factor.

**(d) Recommended upstream fix.** Adopt smbpal's form:

```julia
refrozen = min(snow_melt, refreezing_fraction * snowfall)
```

i.e. cap the refrozen mass by a capacity proportional to the snow *deposited*, not to the snow
*melted*. Optionally add a second cap from an explicit cold-content proxy if `PDDState` ever
gains a pack temperature. Keep `refreezing_fraction` as the tunable, and rename it
`refreezing_capacity_fraction` so the meaning is not misread.

**Port status.** chion reproduces Chion.jl's uncapped form so that WP16 can compare the two
models field for field. It is not the recommended physics; do not carry it into WP19.

---

## D2 — refrozen mass in `snowpack_swe` **and** `smb_ice` **and** `runoff`

This is PLAN §3.2's "check this is not double-counting". The answer is subtler than
double-counting, and worse in practice.

**Where.** `pdd.jl:54-56`:

```julia
snowpack_swe[idx] = available_snow - snow_melt + refrozen
smb_ice[idx]     += snowfall - snow_melt + refrozen - ice_melt
runoff[idx]      += rainfall + snow_melt - refrozen + ice_melt
```

**(a) Is it wrong, or defensible?** Mass is **not** double-counted in the naive sense. Adding
the three increments gives exactly

```
d(smb_ice) + d(runoff) = snowfall + rainfall
```

term for term, with `snow_melt`, `refrozen` and `ice_melt` all cancelling. That closure is
exact and is asserted to 1e-10 in `tests/test_pdd.f90`.

The defect is what `smb_ice` then *means*. Note that

```
d(snowpack_swe) = snowfall - snow_melt + refrozen
d(smb_ice)      = snowfall - snow_melt + refrozen - ice_melt
                = d(snowpack_swe) - ice_melt
```

so `smb_ice` is credited with the entire change in the snow reservoir. `smb_ice` is therefore
the mass change of the **whole column** (snow reservoir plus ice), not the flux delivered to
the ice. The three-reservoir identity that a coupled ice-sheet host would expect,

```
d(smb_ice) + d(runoff) = snowfall + rainfall - d(snowpack_swe)
```

**does not hold**; its residual is exactly `d(snowpack_swe)`. `tests/test_pdd.f90` asserts that
residual *equals* `d(snowpack_swe)`, which makes the check a regression guard: it will break
the moment this is fixed upstream, rather than silently passing.

Two independent pieces of evidence that this is a defect and not a convention:

1. `Chion.jl/src/io.jl:9` documents the variable as
   `long_name = "Net mass forcing to the ice sheet"`. It is not that.
2. **BESSI uses the opposite convention.** In `BESSIModel`, `smb_ice` is incremented *only*
   when mass leaves the bottom of the snowpack or the column is bare —
   `layer_structure.jl:218` (`exported_excess` from a bottom merge), `:337`/`:348`
   (`_continuous_bottom_deplete!`), `step.jl:223` (bare-ice `net_mass_change`), `step.jl:348`
   (`- ice_melt`). Snowfall retained in the pack is *never* credited to `smb_ice` in BESSI.
   So the same output variable means two different things depending on which model produced it,
   and the two are not comparable — which is precisely what WP16 needs to do.

**(b) What smbpal does.** smbpal has no persistent snow reservoir at all. It computes
`smb = sf - runoff` (`smbpal.f90:502`) over an annual window, so the question does not arise;
by construction its SMB is a flux, not a storage change.

**(c) Physical consequence.** `Chion.jl/src/simulation.jl:432` exposes
`ice_sheet_net_forcing_yearly` as the `smb_ice` delta, so this value goes straight into a
coupled ice sheet.

- Over a closed seasonal cycle (pack builds and fully melts) the annual delta is nearly right,
  because `d(snowpack_swe)` over the year is ~0. The error is intra-annual: the ice sheet is
  told it gained the whole winter snowpack in winter and lost it again in summer.
- Any *trend* in the pack leaks directly into the ice-sheet forcing. In a column that never
  melts — which, by D3, means the pack grows forever — `smb_ice` grows without bound and the
  ice sheet is forced with a mass gain that is physically still sitting on the surface as snow.
  For a 300 mm w.e. yr⁻¹ accumulation site that is 30 m w.e. of spurious ice-sheet forcing over
  100 years.
- Monthly output (`monthly_output.jl:77`) differences `smb_ice`, so monthly SMB fields inherit
  the same contamination.

**(d) Recommended upstream fix.** Separate the two budgets. `smb_ice` should be credited only
with mass that actually leaves the snowpack for the ice — in a bulk PDD model that means the
refrozen superimposed ice and the ice melt:

```julia
snowpack_swe[idx] = available_snow - snow_melt          # refrozen leaves the pack
smb_ice[idx]     += refrozen - ice_melt
runoff[idx]      += rainfall + snow_melt - refrozen + ice_melt
```

which restores `d(smb_ice) + d(runoff) + d(snowpack_swe) = snowfall + rainfall` and, as a
bonus, fixes D6 (the refrozen mass no longer re-enters the melt-able reservoir). Note this
also makes `smb_ice` mean the same thing as in BESSI. If the current whole-column quantity is
wanted as well, publish it separately as `column_mass_change`.

**Port status.** chion reproduces Chion.jl. The test asserts *both* the closure that does hold
and the exact size of the violation of the one that does not.

---

## D3 — `snowpack_swe` is uncapped, un-aged and un-densified

**Where.** Nothing anywhere in `pdd.jl` bounds `snowpack_swe`, converts it to ice, or ages it.
Its only sink is melt.

**(a) Is it wrong, or defensible?** Wrong. A bulk PDD scheme needs *some* mechanism to retire
old snow, because the reservoir is what determines whether the next warm season melts snow (at
`ddf_snow = 3`) or ice (at `ddf_ice = 8`). Without one, the accumulation zone accumulates a
reservoir with no upper bound and no memory decay.

**(b) What smbpal does.** smbpal's `calc_ablation_pdd` is stateless: `acc` is the current
year's accumulation and there is no carry-over at all. That is the opposite extreme — arguably
too little memory — but it cannot diverge.

**(c) Physical consequence.** Three distinct effects, all in the same direction (ablation
under-predicted):

1. **Unbounded growth.** At 300 mm w.e. yr⁻¹ with no melt, `snowpack_swe` reaches
   3 × 10⁴ kg m⁻² in 100 years and 3 × 10⁵ kg m⁻² in 1000. There is no equilibrium.
2. **Spurious ablation buffer.** A warming applied after a long cold spin-up must first melt
   the entire accumulated reservoir at `ddf_snow` before any ice melt begins. A 100-year
   spin-up buys ~100 years of suppressed ice ablation. This makes the transient response to a
   warming scenario dependent on the spin-up length, which is a serious modelling artefact.
3. **Feeds D1 and D2.** The unbounded reservoir is what makes uncapped refreezing (D1) and
   snowpack-contaminated `smb_ice` (D2) unbounded rather than merely biased.

There is also a **precision** consequence for the chion port: `docs/PLAN.md` §2.2 specifies
`snowpack_swe` as `wp = sp`, which is correct for a bounded reservoir. At 10⁵ kg m⁻², `sp`
resolution is ~0.008 kg m⁻², so daily snowfall increments below that are lost outright — the
exact failure mode `wp_acc` exists to prevent. The right fix is D3 itself; if D3 is not fixed
upstream, `snowpack_swe` must be promoted to `wp_acc` in chion. Flagged at its declaration in
`snow_pdd.f90`.

**(d) Recommended upstream fix.** Any one of:

- a cap `snowpack_swe = min(snowpack_swe, snowpack_max)` with the excess exported to `smb_ice`
  as firn-to-ice conversion (cheapest, and mirrors BESSI's depth cap in spirit); or
- an explicit conversion timescale, exporting a fixed fraction per year to ice; or
- smbpal's stateless annual framing, if the model is only ever driven annually.

The cap is preferred: it is one line, it bounds the state, and it gives `smb_ice` a physically
meaningful source term once D2 is fixed.

**Port status.** Reproduced. Flagged at the `snowpack_swe` declaration.

---

## D4 — PDD flavour selected implicitly by timestep length

**Where.** `pdd.jl:13-14`:

```julia
_is_monthly_pdd_step(forcing, time_index) = 27.0 <= _step_dt(forcing.dt_days, time_index) <= 32.0
```

Dispatched at `pdd.jl:95` and `runtime.jl:254`, `:281`.

**(a) Is it wrong, or defensible?** Wrong, and PLAN §3.2 already calls it "implicit and
fragile". A change of *physics* — from a truncation `max(T-T₀,0)` to a statistical expectation
integral — is triggered by the *numerics* of the timestep. Consequences:

- A 26-day step and a 28-day step use different melt physics.
- A 365-day (annual) step falls outside the window and uses the truncation, which is exactly
  the regime where the truncation is worst.
- A 5-day or 10-day step uses the truncation, so a model refined from monthly to 10-daily
  forcing silently loses all sub-freezing melt.
- The window is not even a correct test for "monthly": `dt_days` may be a vector
  (`forcing.jl:313`), so different months in the same run can take different branches. A
  365.25/12 = 30.44-day step qualifies, a 26-day February-ish step does not.

**(b) What smbpal does.** smbpal **always** uses the expectation integral, at every timestep
length. `calc_temp_effective` is called for monthly forcing (`smbpal.f90:300`, accumulated as
`teff*30.0`) and for 10-daily forcing (`smbpal.f90:383`, `teff*10.0`). It has no truncated
branch at all. The scheme is only used with sigma > 0, and the integral reduces smoothly to
the truncation as sigma → 0, so there is nothing to switch on.

**(c) Physical consequence.** The truncation `max(T̄ - T₀, 0)` systematically **under-predicts
melt** wherever the mean temperature is near or below freezing, which is most of an ice sheet
for most of the year. With sigma = 5 K:

| T̄ (°C) | simple, K d⁻¹ | PISM, K d⁻¹ | melt missed at `ddf_snow=3` |
|---|---|---|---|
| −10 | 0 | 0.0425 | 0.13 kg m⁻² d⁻¹ |
| −5 | 0 | 0.4166 | 1.25 kg m⁻² d⁻¹ |
| −2 | 0 | 1.152 | 3.46 kg m⁻² d⁻¹ |
| 0 | 0 | 1.995 | 5.98 kg m⁻² d⁻¹ |
| +5 | 5 | 5.417 | 1.25 kg m⁻² d⁻¹ |
| +40 | 40 | 40 | 0 |

(the −5 °C and 0 °C rows are asserted in `tests/test_pdd.f90`). Integrated over a melt season
in the percolation zone this is a large, one-signed error: a daily-forced Chion.jl PDD run
gives systematically less ablation than the same forcing through smbpal, and the difference is
concentrated exactly at the equilibrium-line altitude, where it moves the ELA.

**(d) Recommended upstream fix.** Two changes:

1. Replace `_is_monthly_pdd_step` with an explicit `pdd_method::Symbol` field on `PDDModel`
   (`:simple` | `:pism`), validated in the constructor. Never branch physics on `dt_days`.
2. Make `:pism` the default, matching smbpal. `:simple` is then available as
   `temperature_sigma → 0` in all but name, and should probably be deleted outright.

**Port status.** chion implements (1): `pdd_par_class%pdd_method`, flags `CHION_PDD_SIMPLE` /
`CHION_PDD_PISM`, string mapping via `pdd_method_flag`. This is an approved §4.1 cleanup and
is a **deliberate divergence** from Chion.jl's dispatch. chion keeps `simple` as the default
so that the daily path reproduces Julia bit-for-bit modulo D12; **WP13 should set `pism` as the
namelist default for production**, matching smbpal, and WP16 must run both.

---

## D5 — two entry points apply different physics

**Where.** Four `pdd_step!` methods and two `pdd_monthly_step!` methods, in one file:

| line | signature | monthly detection? | implementation |
|---|---|---|---|
| `pdd.jl:92` | `(model, state, forcing, time_index)` | **yes** | vectorised (monthly) or KA kernel (daily) |
| `pdd.jl:198` | `(swe, smb, runoff, pdd_sum, forcing, time_index, ddf…, scratch)` | **no** | CPU vectorised, always simple |
| `pdd.jl:283` | `(swe, smb, runoff, pdd_sum, forcing, ddf…, scratch)` | **no** | loops all time indices, always simple |
| `pdd.jl:301` | `(swe, smb, runoff, pdd_sum, forcing, time_index, ddf…)` | **no** | KA kernel, always simple |
| `pdd.jl:241` | `pdd_monthly_step!(… , model, scratch)` | n/a | CPU vectorised, always PISM |
| `pdd.jl:334` | `pdd_monthly_step!(… , model)` | n/a | KA kernel, always PISM |
| `pdd.jl:366` | `pdd_step!(model, state, forcing)` | **yes**, via `:92` | loops all time indices |

**(a) Is it wrong, or defensible?** Wrong. `pdd_step!(model, state, forcing)` (`:366`) and
`pdd_step!(swe, smb, runoff, pdd_sum, forcing, ddf…, scratch)` (`:283`) are both "run the whole
forcing", and they give **different answers** on monthly forcing: the first switches to the
expectation integral, the second does not. Which one a caller gets depends on whether they
happen to have a `PDDModel` in hand or the four raw vectors. There is no documentation of the
difference.

Worse, `pdd_step!(model, state, forcing, time_index)` at `:92` is itself inconsistent between
its two branches: the monthly branch calls the CPU-vectorised `pdd_monthly_step!` (`:241`),
while the daily branch falls through to the 9-argument `pdd_step!` at `:301`, which launches a
KernelAbstractions kernel *even on CPU* — a different code path from the vectorised
`pdd_step!` at `:198` that `runtime.jl` actually uses for the same case.

**(b) What smbpal does.** One entry point, `calc_ablation_pdd`, `elemental`, no dispatch.

**(c) Physical consequence.** Reproducibility, not physics per se: the same model, same state
and same forcing produce different mass balance depending on which exported method the caller
reaches. This is the kind of defect that makes a validation harness (WP16) untrustworthy,
because the reference answer depends on the harness's call style.

**(d) Recommended upstream fix.** Collapse to one kernel and one dispatcher. Once D4 gives
`PDDModel` an explicit `pdd_method`, the PDD flavour is a model property, the monthly-vs-daily
split disappears, and `pdd_monthly_step!` can be deleted entirely. Keep exactly two public
methods: `pdd_step!(model, state, forcing, time_index)` and
`pdd_step!(model, state, forcing)`; make everything else internal.

**Port status.** chion has one core (`pdd_column_apply`), one per-column wrapper
(`pdd_column_step`) and one domain driver (`pdd_step`). Removing duplicate entry points that
exist only for GPU dispatch is explicitly allowed by PLAN §4.1.

---

## D6 — refrozen mass re-enters the melt-able snow reservoir

**Where.** `pdd.jl:54`: `snowpack_swe = available_snow - snow_melt + refrozen`.

**(a) Is it wrong, or defensible?** Wrong. Refrozen melt water is **superimposed ice**. It is
denser, it is at the melting point, and it must not subsequently be ablated with the *snow*
degree-day factor. Here it is added straight back as ordinary snow, indistinguishable from the
rest of the reservoir, and can be melted again next step at `ddf_snow` and refrozen again at
`refreezing_fraction`.

**(b) What smbpal does.** `refrz` is removed from runoff and never re-enters any reservoir
(`smb_pdd.f90:47-50`); there is no reservoir to re-enter.

**(c) Physical consequence.** Two effects.

1. **The snow reservoir can never be exhausted.** Whatever the PDD, at most a fraction
   `(1 - refreezing_fraction)` of the pack is removed per step: with the default 0.6, the
   reservoir decays geometrically as 0.4ⁿ. `tests/test_pdd.f90` asserts this directly — a step
   with `T = T₀ + 1000 K` leaves exactly `0.6 × swe` behind, and 21 such steps still leave a
   non-zero pack. The onset of ice melt at a site is therefore delayed by roughly
   `ln(swe/ε)/ln(1/0.4)` steps, ~10–15 days at the start of every melt season, on top of D3's
   much larger delay.
2. **Ablation at the wrong degree-day factor.** Superimposed ice is melted at `ddf_snow = 3`
   rather than `ddf_ice = 8`, understating ablation by a factor of ~2.7 for that mass.

Both push SMB positive, compounding D1 and D3.

**(d) Recommended upstream fix.** The fix for D2 also fixes this: move `refrozen` out of
`snowpack_swe` and into `smb_ice`. If a separate superimposed-ice reservoir is wanted, add a
fifth state field with its own degree-day factor, but that is a redesign; the D2 fix is
sufficient and correct.

---

## D7 — the six-line core exists in three copies, already divergent

**Where.** `_pdd_apply_column_pdd!` (`pdd.jl:32-58`), the vectorised `pdd_step!`
(`pdd.jl:222-237`), and `pdd_monthly_step!` (`pdd.jl:263-279`).

**(a) Is it wrong?** It is hygiene, but the copies have **already diverged**. The kernel
version guards the division:

```julia
remaining_pdd = ddf_snow > zero(ddf_snow) ? max(pdd - snow_melt / ddf_snow, zero(pdd)) : zero(pdd)
```

while both vectorised versions do not:

```julia
@. remaining_pdd = max(pdd - snow_melt / ddf_snow, 0.0)
```

`PDDModel`'s constructor enforces `ddf_snow > 0` (`models.jl:104`), so the unguarded form is
currently unreachable — but the vectorised methods take `ddf_snow` as a **bare positional
argument** (`pdd.jl:206`), bypassing the constructor entirely, so a caller can pass
`ddf_snow = 0` and get `Inf`/`NaN` in one path and `0` in another.

**(b) smbpal.** One `elemental` routine, no duplication.

**(c) Consequence.** Latent `NaN` on an unvalidated call path; guaranteed drift between the
copies as the model is edited.

**(d) Fix.** Delete the two vectorised copies and broadcast the single kernel, or delete the
kernel and keep one vectorised form. Do not expose `ddf_snow`/`ddf_ice`/`refreezing_fraction`
as positional arguments at all — take the `PDDModel`, so validation cannot be bypassed.

**Port status.** chion has one copy. It keeps the `ddf_snow > 0` guard (harmless, and it
documents the inconsistency) on top of `pdd_par_validate`.

---

## D8 — the active-column mask is ignored; the GPU call sites do not compile

**Where.** `runtime.jl:250-305`.

**(a) Two distinct bugs, both real.**

1. **`active_indices` is silently ignored on CPU.** `step_model!(::PDDModel, …)` reaches the
   vectorised `pdd_step!` / `pdd_monthly_step!`, neither of which has an `active_indices`
   parameter; both loop `eachindex(snowpack_swe)`. So `set_active_mask!` has no effect on a
   PDD run — deactivated columns keep accumulating snow, PDDs, runoff and `smb_ice`. BESSI
   honours the mask (`_step_range!(…, model_runtime.active_indices, …)`), so the two models
   behave differently under the same API call.
2. **The GPU branches are `MethodError`s.** `runtime.jl:292-303` calls `pdd_step!` with **ten**
   positional arguments, the tenth being `model_runtime.active_indices`. The only ten-argument
   method is `pdd.jl:198`, which is constrained to `Vector{Float64}` and so cannot match GPU
   arrays; the `AbstractVector` method at `pdd.jl:301` takes nine. Likewise `runtime.jl:281-292`
   calls `pdd_monthly_step!` with eight arguments where `pdd.jl:334` takes seven. Neither GPU
   path can ever have been executed.

**(b) smbpal.** N/A — smbpal is `elemental` over the whole grid and uses `where` masks.

**(c) Consequence.** (1) is a correctness bug for any host that switches columns off — e.g.
yelmox retreating the ice margin. Deactivated columns continue to report growing `smb_ice`,
which is then fed back to the ice sheet. (2) means PDD-on-GPU is untested and non-functional;
it also means the `_pdd_step_kernel!` / `_pdd_monthly_step_kernel!` definitions are dead code
on that path.

**(d) Fix.** Give the PDD steppers an `active_indices` argument and loop over it, exactly as
BESSI does; then fix the GPU call arity. Add a test that asserts a deactivated column's state
is unchanged after `step!` — this would have caught both.

**Port status.** chion honours the mask: `pdd_step` loops `active_idx` only, and
`tests/test_pdd.f90` asserts a deactivated column does not advance. **Deliberate divergence
from Chion.jl**, justified under §4.1 "fix outright bugs, provided they are reported back".

---

## D9 — no `melt` / `refreezing` diagnostics

`PDDState` (`state.jl:180-185`) carries only `snowpack_swe`, `smb_ice`, `runoff`, `pdd_sum`.
`BESSIState` additionally carries `melt` and `refreezing`, and `PDD_OUTPUT_VARS`
(`io.jl:31`) accordingly cannot report them. Total melt and refreezing are the two quantities
a PDD scheme most obviously produces, and smbpal returns both (`calc_ablation_pdd` outputs
`melt` and `refrz`). Recommend adding `melt` and `refreezing` accumulators — they are two
lines in the core and make PDD and BESSI output directly comparable, which WP16 needs.

*Port status:* not added; adding state fields is a design change, not a cleanup. chion's test
recovers `ice_melt` as `d(snowpack_swe) - d(smb_ice)`, which works only because of D2 and
would stop working once D2 is fixed — an argument for D9.

---

## D10 — hard-coded `273.15` and `86400.0`

`pdd.jl:9`, `:11`, `:29`, `:222`, `:225`, `:263` bypass `SnowpackPhysicalConstants`, which
carries both values (`T0`, `seconds_per_day`). PLAN §5 item 3 lists these. Consequence: a run
that changes `T0` or `seconds_per_day` in the constants struct silently gets inconsistent
physics between BESSI and PDD.

*Port status:* chion takes both from `chion_const_class` (`pdd_degree_days`, `pdd_step_mass`).

---

## D11 — one global `temperature_sigma`

`PDDModel` has a single scalar `temperature_sigma = 5.0`. smbpal has three
(`sigma_snow`, `sigma_land`, `sigma_melt`) and selects between them per grid cell by surface
type and temperature (`smbpal.f90:296-298`):

```fortran
smb%now%sigma = smb%par%sigma_snow
where (z_srf > 0 .and. H_ice == 0)           smb%now%sigma = smb%par%sigma_land
where (H_ice  > 0 .and. t2m >= 273.15)       smb%now%sigma = smb%par%sigma_melt
```

The physical argument is that sub-grid and sub-monthly temperature variability is smaller over
a melting, isothermal, high-thermal-inertia snow surface than over bare land. Since `sigma`
enters the expectation integral directly, and `E[T⁺] = sigma/√(2π)` at `T̄ = 0`, this is a
first-order control on melt at the ELA: halving sigma from 5 K to 2.5 K halves the melt at
`T̄ = 0 °C`.

*Recommendation:* make `temperature_sigma` a per-column forcing field rather than a scalar
model parameter, so the host can supply a surface-type-dependent value. Not urgent, but
required before chion-PDD can reproduce smbpal-PDD quantitatively.

*Port status:* scalar, matching Chion.jl. Noted for WP13/WP16.

---

## D12 — Abramowitz–Stegun `_normal_cdf` polynomial

`pdd.jl:16-21` implements Φ(x) with the Abramowitz–Stegun 26.2.17 rational polynomial, whose
absolute error is ~7.5 × 10⁻⁸, and which saturates to exactly 0 or 1 in the far tails.
`erfc` is available in both languages and is exact to round-off. smbpal uses `erfc`
(`smb_pdd.f90:94`) and additionally promotes to `real(8)` around it, with the comment that this
avoids underflow at very low temperatures.

*Port status:* **deliberate deviation, explicitly permitted by PLAN §4.1.** chion uses
`0.5_wp_acc*erfc(-z/sqrt(2))`, evaluated in `wp_acc`. `tests/test_pdd.f90` reimplements the
A–S polynomial solely to measure the difference: max |Φ_erfc − Φ_AS| = **7.45 × 10⁻⁸** over
z ∈ [−6, 6], i.e. an order of magnitude below `sp` round-off and irrelevant to any mass
budget. This is the one place where chion and Chion.jl will not agree bit-for-bit on the PISM
branch; WP16 should set the PDD tolerance accordingly (1e-6 relative is ample).

*Recommendation upstream:* replace `_normal_cdf` with `0.5 * erfc(-x / sqrt(2))`. It is
shorter, faster on any modern libm, and correct in the tails.

---

## S1 — a defect found in **smbpal**, not Chion.jl

Running the comparison the other way turned up one place where **Chion.jl is right and smbpal
is wrong**, which is worth recording because PLAN §3.2 says "where they disagree, smbpal wins"
and this is the exception.

`smb_pdd.f90:42-44`:

```fortran
melt_snow = pdds*mm_snow                                  ! POTENTIAL, uncapped
melt_ice  = max(0.0, (melt_snow-acc)*mm_ice/mm_snow)
melt      = melt_snow + melt_ice
runoff    = melt - refrz
```

`melt_snow` is the *potential* snow melt from the available PDDs and is never capped at the
available accumulation `acc`. When the PDDs exceed what the snow can absorb, `melt` — and
hence `runoff`, and hence `smb = sf - runoff` — counts snow that was not there.

Worked example, `mm_ice = mm_snow`, `acc = A`, PDDs such that `melt_snow = 2A`:

| | smbpal | correct |
|---|---|---|
| snow melt | 2A | A (only A exists) |
| ice melt | A | A |
| total melt | 3A | 2A |
| runoff (f = 0) | 3A | 2A |
| smb = sf − runoff | −2A | −A |

The overstatement is exactly `melt_snow − acc`, i.e. it grows without bound as the site gets
warmer, and it is always in the ablation direction. Note `melt_ice` itself is correct — it is
the same "degree days left over after the snow is gone" construction Chion.jl uses; only the
`melt` and `runoff` totals are wrong.

Chion.jl gets this right via `snow_melt = min(available_snow, ddf_snow*pdd)` (`pdd.jl:49`), and
chion follows Chion.jl here. **Recommended fix in smbpal:**

```fortran
melt_snow_potential = pdds*mm_snow
melt_snow           = min(melt_snow_potential, acc)
melt_ice            = max(0.0, (melt_snow_potential-acc)*mm_ice/mm_snow)
melt                = melt_snow + melt_ice
refrz               = min(melt_snow, acc*f_refrz_max)
runoff              = melt - refrz
```

This should be raised against smbpal before WP19 uses smbpal as the equivalence target, since
it affects the ablation zone directly. Note smbpal's default `abl_method` in yelmox is `itm`,
not `pdd`, so the production impact is limited to PDD configurations.

---

## Deviations of the chion port from Chion.jl

To be folded into `docs/porting_notes.md` under a WP9 heading (this WP was scoped to create
only `src/physics/snow_pdd.f90`, `tests/test_pdd.f90` and this file).

| id | deviation | authority | impact |
|---|---|---|---|
| P1 | `pdd_method = "simple" \| "pism"` replaces the implicit `27 <= dt_days <= 32` trigger | PLAN §4.1 "fix outright bugs"; explicitly requested in the WP9 brief | Same physics for daily and 27–32-day forcing with the default `simple`; different (and better) for every other step length. Set `pism` in WP13 for production. See D4. |
| P2 | `0.5*erfc(-z/sqrt(2))` replaces the Abramowitz–Stegun polynomial | PLAN §4.1, first bullet, names this case | ≤ 7.45e-8 absolute in Φ; below `sp` round-off. See D12. |
| P3 | `T0` and `seconds_per_day` come from `chion_const_class` | PLAN §4.1 "promote magic numbers"; requested in the brief | None at default constants. See D10. |
| P4 | one kernel + one wrapper + one driver, instead of six overlapping entry points | PLAN §4.1 "remove the duplicate/inconsistent entry points that exist only for GPU dispatch" | Removes the D5/D7 ambiguity. |
| P5 | `pdd_step` honours `active_idx` | PLAN §4.1 "fix outright bugs, provided they are reported back" — reported as D8 | Deactivated columns no longer advance. Differs from Chion.jl, which advances them. WP16 must compare with all columns active. |
| P6 | an unrecognized `pdd_method` stops with a message rather than falling through | consistency with `chion_densify_scheme_flag`, `docs/porting_notes.md` D5 | Behaviour differs only for input that was already invalid. |
| P7 | `pdd`, `snowfall`, `rainfall` and all intermediate melt terms are `wp_acc` locals | PLAN §3.1 | Makes the closure identity hold to 1e-10 rather than to `sp`. `snowpack_swe` remains `wp` per PLAN §2.2 — see the precision note under D3. |

**Preserved deliberately (do NOT "fix" in chion):** the uncapped `refreezing_fraction` (D1),
`smb_ice` crediting the snowpack (D2), the unbounded reservoir (D3), and refrozen mass
re-entering the snow reservoir (D6). All four are Chion.jl physics and are reproduced so that
WP16 can compare the two models. They are *not* the recommended physics and must be resolved —
upstream, or by an explicit documented divergence — before WP19 uses chion-PDD for anything.
