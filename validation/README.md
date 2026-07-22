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

Requires both chion builds:

```sh
make all                 # libchion/bin
make all precision=dp    # libchion/bin-dp
```

## What it runs

| target | reference | authority |
|---|---|---|
| BESSI | Chion.jl | authoritative — tight tolerances |
| PDD | Chion.jl | structure only; Chion.jl's PDD is not fully working (PLAN.md §3.2) |
| ITM | smbpal | no Chion.jl ITM exists; runs `test_itm.x`, not a reimplementation |

One forcing file drives both models per target, so a difference is attributable
to the models rather than to two generators drifting apart.

## Why there is a separate Julia environment

Chion.jl's checked-in `Manifest.toml` is stale — `CUDA` is a direct dependency
in its `Project.toml` but absent from the manifest, so `Pkg.instantiate()` fails
inside that repository. This environment `Pkg.develop`s Chion.jl by path and
pins its own manifest, which leaves the reference model's repository untouched
and means a future change to its manifest cannot silently alter what WP16
validates against.

## Why chion is built twice

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

Worst field, dp build: **0.47 ulp** (`temperature`, `density`). `N` is exactly
0 — the layer counts agree at every step of every column. PDD worst: **0.73
ulp** (`pdd_sum`). ITM agrees with smbpal to 1.1e-07 relative at sp and 5.1e-15
at dp.

sp build, reported: worst ~4e-06 relative (`liquid_water`, `density`), first
divergence typically within a few records. That is the cost of `wp = sp`,
measured end to end.

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
| PDD monthly | PISM/Calov–Greve expectation integral | 60 steps at dt=30 d |

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
- **Chion.jl cannot read a CF time axis.** The forcing file carries the time
  axis twice — CF numeric for `chion_grid.x`, and YYYY/MM/DD/HH for Chion.jl.
  See the note in `forcing.jl` and the upstream defect list.

## Not covered

Diurnal substepping is off in both models, deliberately: it changes the albedo
scheme rather than only the shortwave resolution (upstream defects 19 and 21),
and it is the only consumer of `day_of_year` / `solar_longitude_deg`, which the
two drivers derive differently. Comparing it would compare two known-divergent
schemes. Humidity forcing is likewise absent — with it on, Chion.jl's vapour
diagnostic is unclosed (upstream defect 1).
