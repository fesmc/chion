# validation — WP16

Compares chion against its reference models, field by field, and fails on a
threshold.

```sh
julia --project=validation validation/validate.jl          # full
julia --project=validation validation/validate.jl --quick  # 90 d / 12 months
```

First run only:

```sh
julia --project=validation -e 'using Pkg; Pkg.develop(path=joinpath(homedir(),"models","Chion.jl")); Pkg.instantiate()'
```

The reference is whichever Chion.jl the active environment develops; select
another without editing code by pointing `CHION_VALIDATION_PROJECT` at an
environment that develops it (e.g. one per reference commit, each made as
above with `Pkg.develop(path=<Chion.jl worktree>)`). The log's header prints the
reference's path and commit. Options a newer reference adds with defaults that
differ from chion's physics are pinned to chion's choice when present
(`BESSI_SCHEME_PINS` in `runners.jl`: `turbulent_flux_scheme = seb_scheme =
:bessi`, `refreezing_correction = 1` since dev_nils `27113b6`; `longwave_scheme =
:graybody`, `ice_substrate_layers = 0`, no fine near-surface layers since main
`9ec6cc7`).

Requires three chion builds, all with `fpsafe=1`:

```sh
make all fpsafe=1                                  # libchion/bin-fpsafe
make all fpsafe=1 precision=dp                     # libchion/bin-dp-fpsafe
make all fpsafe=1 precision=dp legacy_chion=1      # libchion/bin-dp-legacy-fpsafe
```

`fpsafe=1` compiles with value-safe optimization (`-O2 -fp-model precise` for
ifx/ifort, `-O2` for gfortran) instead of the machine fragment's flags, in its
own directories, so the production `-Ofast` objects are untouched. The gates
resolve fractions of a Float32 ulp and `test_itm.x` checks sp to one ulp; ifx
`-Ofast` (Levante) reorders and contracts arithmetic enough to fail those
checks, which would measure the compiler rather than the port. The acceptance
tests use the same builds (`libchion/bin*-fpsafe/test_*.x`); production
(`libchion.a` for yelmox, `chion_grid.x` runs) keeps the machine's flags.

## What it runs

| target | reference | authority |
|---|---|---|
| BESSI | Chion.jl | authoritative — tight tolerances; four configurations: `albedo = :dynamic`, `:aging` (`6d06af6`, timescales set explicitly, plus `snow_age_days`), `:dynamic` with humidity on (uniform `rh = 0.7`, sea-level pressure; exercises the latent flux and vapour mass of `d0146e1`), `:dynamic` with the thermal ice substrate (`ice_substrate_layers = 5`, `03bb445`; uniform `HI = 1000 m` for chion, D34; plus a dry `bare_ice` column), `:dynamic` with fine near-surface layers (`(0.02, 0.05, 0.10, 0.30)` m, `03bb445`), alone and with the substrate, and `:dynamic` with the cloud-proxy longwave (`longwave_scheme = :cloud_proxy`, `03bb445`; chion's internal TOA), `seb_scheme = :semix` with BESSI turbulence, Chion.jl's surface scheme `seb_scheme = turbulent_flux_scheme = :semix` (humidity on; alone and with the substrate), and diurnal substeps at Chion.jl's calibrated `03bb445` set (8 substeps, 1 K cycle) |
| PDD | Chion.jl, and its own mass closure | authoritative since Chion.jl adopted chion's budget (D23, `ce6a68d`); both `pdd_method`s gated |
| ITM | Chion.jl, and smbpal | Chion.jl's `ITMModel` (ported from chion, `29eb867`): gated at dp+legacy, all 8 written fields (D27 reverted). smbpal, the production reference: runs `test_itm.x`, not a reimplementation |

One forcing file drives both models per target, so a difference is attributable
to the models rather than to two generators drifting apart.

## Why there is a separate Julia environment

Chion.jl's checked-in `Manifest.toml` is stale — `CUDA` is a direct dependency
in its `Project.toml` but absent from the manifest, so `Pkg.instantiate()` fails
inside that repository. This environment `Pkg.develop`s Chion.jl by path and
pins its own manifest, which leaves the reference model's repository untouched
and means a future change to its manifest cannot silently alter what WP16
validates against.

## The three comparisons

| comparison | question | gated? |
|---|---|---|
| chion dp+legacy vs Chion.jl | is the **port** faithful? | **yes** |
| chion sp vs chion dp | what does `wp = sp` cost? | reported |
| chion dp vs chion dp+legacy | what does the gas-constant correction do? | reported |

`legacy_chion=1` reverts chion's deliberate physics corrections to Chion.jl's
values (D24). It exists because "is the port faithful?" and "is the reference
correct?" are different questions: without it, every upstream bug chion fixes
becomes a gate failure, and the only way to stay green is to stop testing those
fields — the harness would weaken exactly as the port improved. The gas-constant
fix alone (D22) would have ungated 15 of BESSI's 18 fields.

PDD needs no switch: Chion.jl `ce6a68d` adopted chion's budget (D23), so the
plain dp build is gated against Chion.jl, once per `pdd_method` (`simple`,
`pism`), with every `&pdd` parameter, including `H_snow_max`, set explicitly on
both sides. chion's output is also gated on the full closure

    snowfall + rainfall == d(snowpack_swe) + d(smb_ice) + d(runoff)

at dp and sp, since it does not depend on the reference.

## Why chion is built at more than one precision

`wp` is a compile-time switch (`precision=sp|dp`, see porting_notes D19).

The **dp build is gated**. chion ships at `wp = sp`, but sp and Chion.jl's
`Float64` do not differ by a bounded per-field offset: the two models share
discrete branch points — the `mass_max` split, the `mass_min` merge,
`surface_has_snow`, the melt gate — and not the round-off that decides which
side of one a given step falls on. Once any branch fires a step apart, the layer
structure diverges by O(1) even though the physics is identical. Gating sp would
therefore fail on a correct port, and PLAN.md's instruction is that a tolerance
needing to be loosened means re-read the code — which would send you chasing a
non-bug. Validating at dp removes the confound.

The **sp build is reported, not gated**. That makes the sp-vs-dp difference a
measured end-to-end number (PLAN.md §3.1 measured candidate expressions in
isolation but never a full run), and the `onset` column shows the record at
which sp first parts company with the reference.

## Tolerance derivation

**The reference files are `float32`.** Chion.jl computes in `Float64` but writes
`Float32` (`src/io.jl`: the output buffers are `Matrix{Float32}` /
`Array{Float32,3}`). Every reference value is quantized before the harness sees
it, so ±0.5 ulp of Float32 separates the two models even if their arithmetic is
bit-identical. **No comparison through these files can resolve below that**,
whatever precision either model runs at. Tolerances are therefore expressed in
ulp of Float32 — the file, not the model, sets the resolution.

At `wp = dp` chion's own arithmetic contributes ~1e-16, four orders below that
floor, so it drops out. What remains is the ±0.5 ulp write quantization plus its
accumulation through fields built from other fields. The gate is **4 ulp of
Float32 (4.8e-07)**, which bounds that with margin while staying far below what
a real divergence produces — a shifted layer split or a flipped melt branch
moves a field by O(1) relative, six-plus orders above this.

The same gate covers the layer-resolved fields. That is a deliberately strong
claim: it asserts the layer *structure* is identical, not merely similar.

### Measured (365 daily steps, 4 columns)

Worst field, port-fidelity gate (dp+legacy vs Chion.jl): **0.47 ulp**
(`temperature`, `density`); the `:aging` configuration also 0.47 ulp (`Tsrf`),
`snow_age_days` exact. Before D31 one `:aging` column diverged: a merge of two
layers at `T0` gave `T0 - 1 ulp` in chion, flipping the `T >= T0` timescale test. `N` is exactly 0 — the layer counts agree at every
step of every column. ITM agrees with smbpal to 1.1e-07 relative at sp and
5.1e-15 at dp. PDD vs Chion.jl (dp, 60 monthly steps): worst 0.37 ulp
(`simple`) and 0.41 ulp (`pism`). ITM vs Chion.jl (dp+legacy, 1095 daily steps):
worst 0.44 ulp (`Tsrf`). PDD mass closure: 2.3e-15 at dp, 3.8e-07 at
sp.

**Re-pointed to dev_nils `27113b6` (WP8), before WP9–WP11:** PDD (both methods)
and ITM unchanged and green (≤0.44 ulp). BESSI fails in every
temperature-driven field, both configurations, as expected from upstream's
harmonic conductance, Calonne conductivity and phase-dependent latent heat:
`temperature`, `Tsrf`, `density` and `bulk_density` from record 1; `mass_w`
2.6e-02 relative, `liquid_water`, `refreezing` ~7e-03, `mass`, `mass_base`,
`smb_ice`, `thickness`, `wet_mass` 1.5–2.7e-03, `runoff`/`melt` 6e-06–3e-05,
`albedo` 1.4e-04 (`:dynamic` only; `:aging` stays 0.30 ulp). `N`,
`snow_age_days`, `sublimation` and `latent_heat_flux_sum` agree exactly — the
last two are identically zero, since humidity is off.

**After WP9–WP11 (dev_nils `27113b6`):** all gated fields pass. WP9 (harmonic
conductance) alone roughly halved the temperature-driven differences
(`temperature` 2.4e-03 → 1.5e-03, `mass_w` 2.6e-02 → 1.5e-02); WP10 (Calonne)
closed the rest: worst 0.47 ulp (`temperature`, `Tsrf`), both configurations.
The humidity-on configuration (WP11) agrees field for field, worst 0.47 ulp
(`temperature`), including `sublimation` (0.34 ulp) and
`latent_heat_flux_sum` (0.44 ulp); every column has a non-zero latent flux and
the `melting` column sublimates. Defect 1 (the unclipped vapour diagnostic)
did not split the two models at `27113b6`: chion reproduced it, so it was gated
rather than reported. Both sides fix it since Chion.jl `03bb445` / chion C1
(the diagnostic is the change applied); `tests/test_bessi.f90` test 1c asserts
closure with the surface layer exhausted every step.

**Re-pointed to main `9ec6cc7` (= `03bb445`, C0), before C1/C2:** run with the
pins `seb_scheme = turbulent_flux_scheme = :bessi`, `longwave_scheme =
:graybody`, `ice_substrate_layers = 0`, `near_surface_layer_max_thicknesses_m
= Inf`, diurnal off, `alpha_ice = 0.3`, `alpha_wet = 0.70`. PDD (both methods)
and ITM green. BESSI red in every configuration (16 + 16 + 18 fields), as
expected from the unswitchable Robin surface boundary: `Tsrf`, `temperature`,
`density`, `bulk_density`, `albedo` from record 1 (`Tsrf` 2.2 K, `temperature`
21 K, `albedo` 0.40 absolute), `thickness` record 2, `mass`/`wet_mass` record
5, melt-season fields from record 127 (`melt` 3.4e-03, `runoff` 6.8e-03,
`refreezing` 0.26, `mass_w` 0.68 relative), `N` from record 140–153; with
humidity on also `sublimation` 2.1e-02 and `latent_heat_flux_sum` 1.2e-02 from
record 1. The same runner against `27113b6` still passes every gate (the
`03bb445` pins are feature-detected and absent there).

**After C1–C2 (main `9ec6cc7`, same pins):** all gated fields pass. Worst
0.47 ulp (`temperature`, `Tsrf`) in each BESSI configuration (dynamic, aging,
humidity on); PDD ≤0.41 ulp, ITM ≤0.44 ulp. Not covered by the pins, so not
gated yet: cloud-proxy longwave, ice substrate, fine near-surface layers,
Julia SEMIX turbulence, diurnal substeps (Stage C3–C8).

**After C3 (ice substrate):** the substrate configuration (5 layers, both
sides) passes every field, worst 0.47 ulp (`Tsrf`); the other configurations
are unchanged (all ≤0.47 ulp). `bare_recover` goes bare only in the melt
season, over ice already brought to T0 by the melting snow (119 bare steps, all
melting), so cold bare ice is covered by the dry `bare_ice` column: 220 bare
steps below T0 without melt (Tsrf down to 247 K), melt only after re-warming.

**After C4 (fine near-surface layers):** two configurations, fine layers
(0.02, 0.05, 0.10, 0.30 m) on the default forcing and fine layers + ice
substrate on the substrate forcing, pass every field, worst 0.48 ulp
(`bulk_density`); the remesh uses chion's exact mixing mean (D31) and shows no
branch flips, so D31 stays out of `legacy_chion`. Coverage: layers 1-4 at
their targets on all 1334 column-steps with a layer below; cap-down into layer
5 on 336 `cold_dry` steps; fill-up from layer 5 on 34 `melting` steps; N never
exceeds 5 and layer 5 reaches 11646 kg m⁻² (`ntot_capacity`, depth cap), the
upstream behaviour C4b changes (reverted under legacy).

**After C4b (split/merge below the fine layers, D32):** unchanged, every
configuration green at the same worst ulp; the fine-layer configurations still
see Chion.jl's unsplit layer 5 under `legacy_chion`. The chion behaviour is
covered by `tests/test_layers` and `tests/test_bessi` (test 15) in the
non-legacy builds.

**After C5–C8 (main `9ec6cc7`):** every configuration green, worst 0.48 ulp. New
configurations: cloud-proxy longwave (worst 0.48 ulp, `refreezing`; 241 days on the
shortwave proxy, 124 on the night cloudiness, Tsrf up to 4.6 K from graybody);
`seb_scheme = semix` with BESSI turbulence (0.47 ulp; 0.6 K); `seb = turbulent_flux_scheme =
semix` with humidity on, alone (0.47 ulp; 2.7 K from BESSI turbulence) and over the
substrate (0.47 ulp; 6.2 K), under `legacy_chion` (Julia's 287.05 and bare-ice `Lv`,
D35/D38); diurnal substeps at the calibrated set (0.47 ulp; 330 column-days split into 8
substeps, 240 unsplit polar-night days with shortwave on the D39 legacy path; 6.3 K from
the daily run).

Reported, not gated:

- **`wp = sp` costs** ~4e-06 relative worst case, first divergence typically
  within a few records.
- **The gas-constant correction** moves the density-driven fields by 2-3%
  against Chion.jl's `8.13`. Integrated over a 10-year column run it is <1% on
  every cumulative quantity, because densification self-limits against the
  `(rho_i - rho)` gap and the `rho_i` cap.

## Coverage is asserted, not assumed

An agreement result is only as strong as the coverage behind it — two models
agree perfectly on a column where nothing happens. The shipped column example
had exactly that problem before the D18 retuning (peak layer count 1, no
refreezing) while still looking healthy. So the harness reports what each
column actually reached and asserts it:

| column | exercises | measured |
|---|---|---|
| `cold_dry` | split/merge, densification, no melt | 6 layers, 300→379 kg m⁻³, melt 0 |
| `melting` | energy solve, melt, percolation, refreezing | 9 layers, melt 547, refreeze 232, runoff 678 |
| `bare_recover` | ablation to bare ice, early-return bare path, recovery | N→0 then back to 2, melt 10139 |
| `ntot_capacity` | `Ntot` capacity, bottom merge, depth cap | N=15, base export 8350 kg m⁻² |
| ice substrate `bare_ice` | dry bare ice over the substrate: below T0 without melt, melts only once re-warmed | asserted |
| fine near-surface layers | layers 1-4 at target thickness; cap-down (`cold_dry`) and fill-up (`melting`) fire; layer 5 unsplit (N ≤ 5) | asserted |
| BESSI `:aging` | snow ages; relaxes below `alpha_dry`; ages on a melting surface (`tau_melt`) | asserted |
| PDD monthly | `simple` and `pism` (Calov–Greve) integrals | 60 steps at dt=30 d each |
| ITM `ice_ablation` | ice background; melts out, then bare-ice melt; rain refreezing | asserted |
| ITM `ice_accum_cap` | `H_snow_max` cap, excess to ice | asserted |
| ITM `land_seasonal`, `ocean` | land (tundra PDDs) and ocean backgrounds | asserted |

This caught two real weaknesses: three of the four columns originally never
split a layer, and `mass_base` was identically zero on both sides — so its
"0.00 ulp" pass was vacuous until `ntot_capacity` was tuned to reach the cap.

## Things the harness has to work around

- **Axis order** (D14). Chion.jl declares `("t","x","y")`; chion writes
  `(time,yc,xc)`. `read_canonical` permutes by dimension *name*, so neither
  writer's order is hard-coded.
- **`MV` vs NaN** (D15). Both map to `missing`, and the two files must agree
  about *which* cells are missing before any value is compared.
- **Record alignment.** chion writes an initial pre-step record; Chion.jl does
  not. chion record k+1 ↔ Chion.jl record k, asserted on counts.
- **Chion.jl writes no time coordinate variable** — only a bare `t` dimension,
  so its output cannot be interpreted without knowing the forcing that produced
  it. Alignment is checked structurally as a result.
- **Chion.jl's `load_forcing_file` reads no ITM ice thickness or annual PDDs.**
  The ITM forcing carries them as `HI` / `PDDA`; `run_julia_itm` merges them
  into the loaded forcing, and `chion_grid.x` reads them via `name_hice` /
  `name_pdds`.
- **ITM output names differ.** chion writes the cumulative accumulators under
  BESSI's names (`melt`, `runoff`, `refreezing`, `smb_total`, `albedo`),
  Chion.jl as `*_cum`, `smb_cum` and `alb_s`; `ITM_PAIRS` maps them.
- **Chion.jl cannot read a CF time axis.** The forcing file carries the time
  axis twice — CF numeric for `chion_grid.x`, and YYYY/MM/DD/HH for Chion.jl.
  See the note in `forcing.jl` and the upstream defect list.

## Not covered

Diurnal substepping is off in every configuration but its own (Chion.jl's
calibrated set: 8 substeps, 1 K cycle capped at 1 K, none below -8 C). Snowfall
brightening stays non-linear under substepping (upstream defect 21), the same on
both sides; a day that is not split keeps its forcing in chion (D39, legacy
reverts). The elevation-dependent temperature amplitude (Chion.jl `d0146e1`) is
covered by `test_wp7` only (gradient 0 in the calibrated set).

`day_of_year` and `solar_longitude_deg` (diurnal substeps, cloud-proxy TOA) agree
by construction: Chion.jl derives both from the calendar axis, `chion_grid.x` the
day as `modulo(time,365)+1` and the longitude with Chion.jl's calendar-day formula
(`calendar_solar_longitude_deg`), identical on the one-year axis from 1 January
2000. A multi-year forcing would part at the leap day. Humidity forcing is uniform (`rh_default`; `chion_grid.x` has no
humidity reader) and absent from the first two configurations.
