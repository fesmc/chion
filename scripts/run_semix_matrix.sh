#!/usr/bin/env bash
# Full SEMIX configuration matrix, GRL-16KM.
#   albedo_scheme x seb_scheme x rh_default, at Ntot=1, plus the full config at Ntot=15.
set -euo pipefail
chion=/Users/alrobi001/models/chion
bin=$chion/libchion/bin/chion_grid.x

#      name              albedo   seb     rh    Ntot
configs=(
  "m_dyn_bessi_dry     dynamic  bessi   0.0   1"
  "m_dyn_bessi_rh      dynamic  bessi   0.7   1"
  "m_dyn_semix_dry     dynamic  semix   0.0   1"
  "m_dyn_semix_rh      dynamic  semix   0.7   1"
  "m_alb_bessi_dry     semix    bessi   0.0   1"
  "m_alb_bessi_rh      semix    bessi   0.7   1"
  "m_alb_semix_dry     semix    semix   0.0   1"
  "m_alb_semix_rh      semix    semix   0.7   1"
  "m_full_n15          semix    semix   0.7  15"
)

for cfg in "${configs[@]}"; do
  read -r name alb seb rh nt <<< "$cfg"
  d=$chion/output/$name
  rm -rf "$d"; mkdir -p "$d/input"
  for f in "$chion"/input/*; do ln -sfn "$f" "$d/input/$(basename "$f")"; done
  rm -f "$d/input/chion_phys_const.nml"
  sed -e "s/^\( *seb_scheme *= *\)\"bessi\"/\1\"$seb\"/" \
      -e "s/^\( *albedo_scheme *= *\)\"dynamic\"/\1\"$alb\"/" \
      "$chion/input/chion_phys_const.nml" > "$d/input/chion_phys_const.nml"
  ln -sfn "$chion/output/cmp_maps" "$d/maps"
  sed -e "s/^\( *Ntot *= *\)[0-9]*/\1$nt/" \
      -e "s/^\( *rh_default *= *\)0\.0/\1$rh/" \
      "$chion/par/chion_grl16.nml" > "$d/chion_grl16.nml"

  echo "=== $name : albedo=$alb seb=$seb rh=$rh Ntot=$nt"
  ( cd "$d" && OMP_NUM_THREADS=1 "$bin" "$d/chion_grl16.nml" > run.log 2>&1 )
  grep -m1 'albedo_scheme =' "$d/run.log" || true
  grep -m1 'seb_scheme =' "$d/run.log" || true
  grep -m1 'rh_default =' "$d/run.log" || true
  grep -m1 'wall time' "$d/run.log" || true
done
echo "=== matrix complete"
