#!/usr/bin/env bash
# Greenland GRL-16KM model/layer comparison (speed + SMB skill).
#
# Runs, serially and identically forced, 50-yr spin-ups for:
#   BESSI Ntot = 15 (default), 7, 5, 3   -- layer-count degradation study
#   PDD, ITM                             -- alternative model comparison
#
# All share the GRL-16KM default par file (par/chion_grl16.nml): same forcing,
# same diurnal-substep BESSI physics, only `model` and `Ntot` change. Serial
# (OMP_NUM_THREADS=1) so wall time and per-column-step are directly comparable.
#
# Output: output/cmp_<name>/{chion_grl16.nc,run.log}, sharing one maps/ cache.
set -euo pipefail

chion=/Users/alrobi001/models/chion
base_par=$chion/par/chion_grl16.nml
bin=$chion/libchion/bin/chion_grid.x
out=$chion/output
mapsrc=$chion/output/grl16_bessi/maps          # pre-built ERA5 weight cache

# Shared maps cache so no run pays the conservative-regrid cost.
mkdir -p "$out/cmp_maps"
cp -n "$mapsrc"/*.nc "$out/cmp_maps"/ 2>/dev/null || true

# name  model  Ntot
configs=(
  "bessi_n15 bessi 15"
  "bessi_n7  bessi 7"
  "bessi_n5  bessi 5"
  "bessi_n3  bessi 3"
  "pdd       pdd   15"
  "itm       itm   15"
)

export OMP_NUM_THREADS=1

for cfg in "${configs[@]}"; do
  read -r name model ntot <<< "$cfg"
  rundir=$out/cmp_$name
  echo "=== $name (model=$model Ntot=$ntot) -> $rundir"
  rm -rf "$rundir"
  mkdir -p "$rundir"
  ln -sfn "$chion/input" "$rundir/input"
  ln -sfn "$out/cmp_maps" "$rundir/maps"

  # Patch model and Ntot into a per-run par file.
  sed -e "s/^\( *model *= *\)\"[a-z]*\"/\1\"$model\"/" \
      -e "s/^\( *Ntot *= *\)[0-9]*/\1$ntot/" \
      "$base_par" > "$rundir/chion_grl16.nml"

  ( cd "$rundir" && "$bin" "$rundir/chion_grl16.nml" ) > "$rundir/run.log" 2>&1
  # Echo the summary block for progress visibility.
  tail -n 12 "$rundir/run.log" | sed 's/^/    /'
done

echo "=== all runs complete"
