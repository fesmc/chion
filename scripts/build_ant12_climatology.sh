#!/usr/bin/env bash
# ============================================================================
# Build 12-month RACMO2.4/ANT-12 climatologies (1981-2010) for the chion
# Antarctica domain, conservatively regridded onto each target ANT grid.
#
# Source: CORDEX ANT-12 monthly time series downloaded to ~/data/racmo/<var>/
#         (tas, pr, rsds, clt; three decade files each, 1981-2010).
#
# Stage 1 -- native rotated-pole climatology (one cdo pipeline per variable):
#   - mergetime the three decade files       -> continuous 360-month series
#   - accumulation -> mean flux:  divdpm (/days-per-month) then divc,86400 (/s)
#     applied to pr [kg m-2 / month] and rsds [J m-2 / month]
#   - clt [%] -> fraction:        divc,100
#   - ymonmean                                -> 12-month climatology
#   - chname to the standardized names (t2m, pr, swd, tcc)
#
# Stage 2 -- conservative regrid onto each ANT grid (cdo remapcon):
#   The target grids are the standard ANT-XXKM south polar stereographic grids
#   (central meridian 0, standard parallel 71 S), described in METERS via an
#   explicit griddes -- the BedMachine files carry xc/yc in "kilometers", which
#   cdo misreads and mis-scales, so we specify the grid here instead. remapcon
#   is flux-conserving, which is why the regrid is done conservatively rather
#   than by sampling. (coords cannot yet remap from a rotated-pole source; see
#   the fesm-utils issue linked in libs/domains/chion_domain.f90.)
#
# Output (read directly by load_antarctica):
#   ~/data/racmo/clim/ANT-12_RACMO24P_monclim_1981-2010.nc     (native, rotated)
#   ~/data/racmo/clim/<GRID>_RACMO24P_monclim_1981-2010.nc     (per ANT grid)
#   with t2m [K], pr [kg m-2 s-1], swd [W m-2], tcc [1] on (month, yc, xc).
# ============================================================================
set -euo pipefail

root="$HOME/data/racmo"
out_dir="$root/clim"
tmp_dir="$root/clim/tmp"
mkdir -p "$out_dir" "$tmp_dir"

decades="198101-199012 199101-200012 200101-201012"
stem="ANT-12_ERA5_evaluation_r1i1p1f1_UU-IMAU_RACMO24P-NN_v1-r1_mon"

files_for() {  # var -> space-separated list of its three decade files
    local v="$1" d
    for d in $decades; do echo "$root/$v/${v}_${stem}_${d}.nc"; done
}

# ---- Stage 1: native rotated-pole climatology --------------------------------
echo "== tas -> t2m [K] =="
cdo -O -chname,tas,t2m -ymonmean -mergetime $(files_for tas) "$tmp_dir/t2m.nc"

echo "== pr [kg m-2/mon] -> pr [kg m-2 s-1] =="
cdo -O -setattribute,pr@units="kg m-2 s-1" -ymonmean -divc,86400 -divdpm -mergetime $(files_for pr) "$tmp_dir/pr.nc"

echo "== rsds [J m-2/mon] -> swd [W m-2] =="
cdo -O -setattribute,swd@units="W m-2" -chname,rsds,swd -ymonmean -divc,86400 -divdpm -mergetime $(files_for rsds) "$tmp_dir/swd.nc"

echo "== clt [%] -> tcc [fraction] =="
# units must be a TEXT attribute -- a bare "1" is stored numerically by cdo and
# ncio then fails reading it (NetCDF: convert between text & numbers).
cdo -O -setattribute,tcc@units="fraction" -chname,clt,tcc -ymonmean -divc,100 -mergetime $(files_for clt) "$tmp_dir/tcc.nc"

native="$out_dir/ANT-12_RACMO24P_monclim_1981-2010.nc"
echo "== merge native -> $native =="
cdo -O merge "$tmp_dir/t2m.nc" "$tmp_dir/pr.nc" "$tmp_dir/swd.nc" "$tmp_dir/tcc.nc" "$native"
rm -rf "$tmp_dir"

# ---- Stage 2: conservative regrid onto each ANT grid -------------------------
# grid_name  xsize  inc(m)   first(m)     (square grids: y = x)
regrid_to() {
    local grid="$1" n="$2" inc="$3" first="$4"
    local gd="$out_dir/griddes_${grid}.txt"
    cat > "$gd" <<EOF
gridtype  = projection
gridsize  = $((n*n))
xsize     = $n
ysize     = $n
xname     = xc
yname     = yc
xunits    = "m"
yunits    = "m"
xfirst    = $first
xinc      = $inc
yfirst    = $first
yinc      = $inc
grid_mapping_name = polar_stereographic
straight_vertical_longitude_from_pole = 0.
latitude_of_projection_origin = -90.
standard_parallel = -71.
false_easting = 0.
false_northing = 0.
EOF
    local out="$out_dir/${grid}_RACMO24P_monclim_1981-2010.nc"
    echo "== remapcon -> $out =="
    cdo -O remapcon,"$gd" "$native" "$out"
}

regrid_to ANT-32KM 191 32000 -3040000
regrid_to ANT-16KM 381 16000 -3040000
regrid_to ANT-8KM  761  8000 -3040000

echo "== done =="
ls -1 "$out_dir"/*_RACMO24P_monclim_1981-2010.nc
