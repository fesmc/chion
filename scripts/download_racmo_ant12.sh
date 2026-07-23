#!/usr/bin/env bash
# ============================================================================
# Download the CORDEX ANT-12 RACMO2.4 (ERA5-forced) monthly variables that the
# chion Antarctica domain needs, for the 1981-2010 climatology period.
#
# Source: DMI THREDDS fileServer (CORDEX PolarRes, UU-IMAU RACMO24P-NN).
# Variables: tas (air temperature), pr (precip), rsds (SW down), clt (cloud).
# Period:    three decade files each -> 1981-2010  (~1.0 GB total).
#
# Files land in ~/data/racmo/<var>/, which scripts/build_ant12_climatology.sh
# then averages into the 12-month climatology and regrids onto the ANT grids.
# ============================================================================
set -u

base="https://cordex.dmi.dk/thredds/fileServer/esg_cordex/PolarRes/ANT-12/UU-IMAU/ERA5/evaluation/r1i1p1f1/RACMO24P-NN/v1-r1/mon"
ver="v20260116"
dest="$HOME/data/racmo"
vars="tas pr rsds clt"
decades="198101-199012 199101-200012 200101-201012"

for v in $vars; do
    mkdir -p "$dest/$v"
    for d in $decades; do
        fn="${v}_ANT-12_ERA5_evaluation_r1i1p1f1_UU-IMAU_RACMO24P-NN_v1-r1_mon_${d}.nc"
        out="$dest/$v/$fn"
        if [ -s "$out" ]; then echo "SKIP (exists): $v/$fn"; continue; fi
        echo "GET $v/$fn"
        curl -f -L --retry 3 --retry-delay 5 -s -S -o "$out" "$base/$v/$ver/$fn" \
            && echo "  OK  $(du -h "$out" | cut -f1)  $v/$fn" \
            || { echo "  FAIL $v/$fn"; rm -f "$out"; }
    done
done
echo "== done =="; du -sh "$dest"
