# Hand-off: cloud-proxy longwave recalibration (BESSI)

Status: open, decided (option a), not started. Written 2026-10-10 after the Chion.jl
9ec6cc7 sync (PR #1). BESSI stays an ongoing project; yelmox uses chion ITM by default.

## Problem

In yelmox (Greenland, CESM2-WACCM ssp585 2015-2300) the synced BESSI loses -0.80 m SLE
(ITM -0.34, graybody alpha_ice 0.45 -0.28). Attribution: about 2/3 of the extra loss is the
cloud-proxy longwave (`longwave_scheme = "cloud_proxy"`) over snow in the percolation zone
(1500-3000 m); the D40 thin-snow blend adds about 1/3; alpha_ice is minor.

Checked against MAR3.9-CESM2 (summer JJA, ice mask, yelmox-like inputs SW = 0.6 TOA):

| band [m] | proxy - MAR [W m-2] | graybody - MAR | dLWD/dT MAR / proxy / graybody [W m-2 K-1] |
|---|---|---|---|
| 0-1000 | +0 | -9 | 6.4 / 5.0 / 3.8 |
| 1500-2000 | +9 | -3 | 5.6 / 4.8 / 3.6 |
| 2000-2500 | +12 | -1 | 5.0 / 4.7 / 3.5 |
| 2500-3000 | +14 | +1 | 4.4 / 4.6 / 3.4 |
| >3000 | +16 | +3 | 4.3 / 4.4 / 3.3 |

- The proxy's temperature sensitivity is right (at or below MAR); graybody's is about 30 % too low.
- The defect is the summer LEVEL above 1500 m (+9 to +16 W m-2): cloudiness
  n = 1 - SW/(TOA*(lw_clear_sky_transmissivity + lw_clear_sky_transmissivity_per_km*z/km))
  depends only on elevation when SW = 0.6 TOA, and is too high aloft. Annual means agree within about 5 W m-2.
- So graybody is not the fix: it gets the present level right but under-responds to warming.

## Plan (option a)

1. Recalibrate the cloud-proxy `lw_*` parameters (mainly `lw_clear_sky_transmissivity`,
   `lw_clear_sky_transmissivity_per_km`, possibly `lw_emissivity_cloud_slope`) against MAR summer
   LWD by elevation band, for yelmox-like inputs (SW = trans_sw*TOA, MAR TT). Keep the
   temperature sensitivity (target: MAR dLWD/dT per band). Script to start from:
   `/work/ba1442/robinson/chion-runs/marproj/scripts/lw_check.jl` (Levante).
2. Recalibrate `alpha_ice` against MAR present day (GRL-8KM, yelmox fixed geometry) with the new
   `lw_*` (previous: 0.55 with default lw_*, graybody optimum 0.45).
3. Rerun yelmox ssp585 GRL-8KM (1-kyr spin-up + 2015-2300) and compare with ITM and the runs below.
4. Put the new values in yelmox's Greenland par files (`&chion_const`), not chion's defaults
   (they are Chion.jl's); add the issue to `docs/upstream_chionjl_issues.md` for Nils.
5. Open question to keep in mind: whether the D40 thin-snow blend (`swe_crit_albedo`) should stay
   on (-0.27 m of the excess).

## Data and earlier runs

- MAR3.9 ISMIP6 CESM2 1950-2100 monthly (SMB, ME, RU, SF, RF, TT, SWD, LWD, AL2, ...):
  `/work/ba1442/robinson/data/MAR-CESM2/v39/` (Levante).
- Standalone MAR-forced comparison (driver prescribing MAR LWD, slice runs, analysis):
  `/work/ba1442/robinson/chion-runs/marproj/` (`analysis_mar/summary_mar.md`, `analysis_mar/lw_check.md`).
- yelmox attribution runs: `/work/ba1442/robinson/yelmox-chion-sync/output/sync/attr/`
  (scripts `cmp_sync/submit_attr.sh`, par `cmp_sync/yelmox_esm_Greenland_attr.nml`).
  BESSI runs in yelmox use `surface_chion.time_equil=20` (converged; 100 yr needed only by ITM).

| yelmox ssp585 run | dV_sle 2300 [m] | ablation area 2300 [Mkm2] |
|---|---|---|
| ITM | -0.34 | 0.67 |
| BESSI pre-sync | -0.15 | 0.35 |
| BESSI cloud proxy, alpha_ice 0.55 | -0.80 | 1.12 |
| graybody, alpha_ice 0.45 | -0.28 | 0.54 |
| cloud proxy 0.55, blend off | -0.53 | 1.15 |
| cloud proxy, alpha_ice 0.50 | -0.90 | 1.10 |

Same-forcing standalone (MAR's own TT/SWD/LWD/SF/RF), dSMB 2081-2100 vs 1961-1990 as % of MAR:
new BESSI 65-85 % (own cloud-proxy LW 82 %), pre-sync 73 %, ITM 47 %. Refreezing is about half of
MAR's at present (no melt above about 1800 m in chion).
