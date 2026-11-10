"""
Fit 3 OLS best-fit lines 
to the BG FWV vs gradSSH derived from ORAS5
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray as xr
import matplotlib.colors as colors
from matplotlib.patches import Rectangle

# Append custom module paths
PPTHN = None
if 'PPTHN' not in locals() or PPTHN is None:
  cwd = os.getcwd()
  parts = cwd.split(os.sep)
  if 'python' in parts:
    idx = parts.index('python')
    PPTHN = os.sep + os.path.join(*parts[:idx + 1])
  else:
    raise RuntimeError("Directory 'python' not found in current working directory path.")

sys.path.extend([
    os.path.join(PPTHN, 'MyPython', 'hycom_utils'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])

from mod_utils_fig import bottom_text, colorbar_horiz
import mod_time as mtime

use_deseaz = True

pthdata='/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data_misc/BeaufortGyre_FWanalysis_ORAS5'
flin = 'BG_1980-2025_FWV_and_SSHG.nc'
flin = 'BG_1980-2025_FWV_and_SSHG_deseas.nc'

dflin = os.path.join(pthdata, flin)

cff_fw = 1e-3
cff_ssh = 1.e7 
with xr.open_dataset(dflin) as dsin:
  YRS = dsin["year"].values.astype(int)
  MMS = dsin["month"].values.astype(int)
  if use_deseaz:
    gSSH = dsin["sshg_deseas"].values * cff_ssh           # ssh grad m/m
    FWV  = dsin["fwv_deseas"].values * cff_fw  # km 3
  else:
    gSSH = dsin["sshg"].values * cff_ssh           # ssh grad m/m
    FWV  = dsin["fwv_integrated_oras"].values * cff_fw  # km 3


# Define clusters of points for 3-trends:
CLSSH = {
  "1" : [0, 5],
  "2" : [4.5, 9.1],
  "3" : [5.8, 12]
}

CLFWV = {
  "1" : [15, 18.5],
  "2" : [16, 22],
  "3" : [21, 23]
}

from scipy.stats import linregress, t

FIT = {}

for irg in range(1, 4):

  ssh_bin = CLSSH[str(irg)]
  fwv_bin = CLFWV[str(irg)]

  # Select points within the region
  msk = (
      (gSSH >= ssh_bin[0]) & (gSSH <= ssh_bin[1]) &
      (FWV >= fwv_bin[0]) & (FWV <= fwv_bin[1])
  )

  x = gSSH[msk]
  y = FWV[msk]

  # OLS linear regression
  res = linregress(x, y)

  # 95% CI for slope
  dof = len(x) - 2
  tcrit = t.ppf(0.975, dof)

  slope_ci = (
      res.slope - tcrit * res.stderr,
      res.slope + tcrit * res.stderr
  )

  FIT[irg] = {
      'slope': res.slope,
      'intercept': res.intercept,
      'r2': res.rvalue**2,
      'slope_ci': slope_ci,
      'npts': len(x),
      'x': x,
      'y': y
  }

  print(f"Region {irg}: N={len(x)}, "
        f"slope={res.slope:.4f}, "
        f"intercept={res.intercept:.4f}, "
        f"  95% CI = [{slope_ci[0]:.4f}, {slope_ci[1]:.4f}] "
        f"R2={res.rvalue**2:.4f}")


nrecs = len(YRS)
rmin = np.min(YRS)
rmax = np.max(YRS)
norm = plt.Normalize(rmin, rmax)
clrmp = colors.Normalize(rmin, rmax)

plt.ion()

fig1 = plt.figure(1, figsize=(9, 9))
plt.clf()

ax1 = plt.axes([0.1, 0.15, 0.8, 0.8])

sc = ax1.scatter(gSSH, FWV, s=10, c=YRS, cmap='jet',
                 norm=norm, marker='o')

# Plot regions:
for irg in range(1,4):
  ssh_bin = CLSSH[str(irg)]
  fwv_bin = CLFWV[str(irg)]

  rect = Rectangle(
      (ssh_bin[0], fwv_bin[0]),
      ssh_bin[1] - ssh_bin[0],
      fwv_bin[1] - fwv_bin[0],
      linewidth=1.5,
      edgecolor='k',
      facecolor='none',
      linestyle='--'
  )

  ax1.add_patch(rect)  

  # Plot least-sq. trends
  for irg in range(1, 4):
    ssh_bin = CLSSH[str(irg)]
    res = FIT[irg]

    xx = np.linspace(ssh_bin[0], ssh_bin[1], 100)
    yy = res['slope'] * xx + res['intercept']

    ax1.plot(xx, yy, 'k-', lw=2)

    ax1.text(
        np.min(ssh_bin),
        np.max(CLFWV[str(irg)]),
        f"R{irg}",
        fontsize=12,
        ha='right'
    )



clb = colorbar_horiz(fig1, ax1, sc, rmin=rmin, rmax=rmax, extd='neither', cb_gap=0.07, decim=0)


ax1.set_xlabel(r'$\nabla SSH$ ($\times 10^{-7}$ m/m)')
ax1.set_ylabel('FWV (10$^3$ km$^3$)')
ax1.grid(True)
ax1.set_title(f"{flin} use deseas.={use_deseaz}")

btx = 'fit_OLS_FWV_gradSSH_BG.py'
bottom_text(btx, pos=[0.02, 0.02], fsz=8)


