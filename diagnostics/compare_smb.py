#!/usr/bin/env python3
# chion vs MAR surface-SMB comparison figure.
#
#   python3 diagnostics/compare_smb.py <chion_output.nc> [mar_file.nc] [out.png]
#
# Surface SMB (MAR definition) = precip - runoff - sublimation, taken as the
# final-year value from a spun-up chion run. Writes a 4-panel PNG: chion map,
# MAR map, difference, and a scatter with bias/RMSE/R2.

import sys
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

chion_file = sys.argv[1] if len(sys.argv) > 1 else "chion_grl16.nc"
mar_file = sys.argv[2] if len(sys.argv) > 2 else \
    "/Users/alrobi001/models/ice_data/Greenland/GRL-16KM/GRL-16KM_MARv3.11-ERA_monmean_1961-1990.nc"
out_png = sys.argv[3] if len(sys.argv) > 3 else "smb_chion_vs_mar.png"

m = nc.Dataset(mar_file)
mask = m["mask"][:]
ice = mask > 50
mar = np.asarray(m["smb"][:].mean(0) * 365.0)          # mm/yr
accum = np.asarray((m["sf"][:] + m["rf"][:]).mean(0) * 365.0)

d = nc.Dataset(chion_file)
ru = d["runoff"][:]
su = d["sublimation"][:]
chion = np.asarray(accum - (ru[-1] - ru[-2]) - (su[-1] - su[-2]))   # final-year surface SMB

x, y = chion[ice], mar[ice]
bias = np.mean(x - y)
rmse = np.sqrt(np.mean((x - y) ** 2))
r2 = 1 - np.sum((x - y) ** 2) / np.sum((y - np.mean(y)) ** 2)

chion_m = np.where(ice, chion, np.nan).T          # (yc,xc) for display
mar_m = np.where(ice, mar, np.nan).T
diff_m = chion_m - mar_m

fig, ax = plt.subplots(1, 4, figsize=(17, 6))
kw = dict(origin="lower", cmap="RdBu", vmin=-2000, vmax=2000)

im = ax[0].imshow(chion_m, **kw); ax[0].set_title("chion surface SMB")
ax[1].imshow(mar_m, **kw); ax[1].set_title("MAR v3.11 SMB")
imd = ax[2].imshow(diff_m, origin="lower", cmap="PuOr", vmin=-800, vmax=800)
ax[2].set_title("chion - MAR")
for a in ax[:3]:
    a.set_xticks([]); a.set_yticks([])
fig.colorbar(im, ax=ax[1], fraction=0.046, label="mm w.e. yr$^{-1}$")
fig.colorbar(imd, ax=ax[2], fraction=0.046, label="mm w.e. yr$^{-1}$")

ax[3].scatter(y, x, s=2, alpha=0.15, color="#333")
lim = [-4000, 3200]
ax[3].plot(lim, lim, "r-", lw=1)
ax[3].set_xlim(lim); ax[3].set_ylim(lim)
ax[3].set_xlabel("MAR SMB"); ax[3].set_ylabel("chion SMB")
ax[3].set_title(f"bias {bias:+.0f}   RMSE {rmse:.0f}   R$^2$ {r2:.2f}")
ax[3].set_aspect("equal")

fig.suptitle(f"{chion_file}   (ice cells, mm w.e. yr$^{{-1}}$)", y=0.98)
fig.tight_layout()
fig.savefig(out_png, dpi=110)
print(f"wrote {out_png}   bias={bias:+.1f}  RMSE={rmse:.1f}  R2={r2:.2f}")
