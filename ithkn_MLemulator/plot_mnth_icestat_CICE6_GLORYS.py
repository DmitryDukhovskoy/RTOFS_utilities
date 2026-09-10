"""
  Plot monthly ice statistics (ice vol and ice thkn overall mean over ice-covered area)
  from JEDI / SOCA CICE6
  and GLORYS reanalysis

  SOCA CICE6:
  calc_mean_ithkn_SOCAcice6.py
  called by bash scripts

  GLORYS
  derived on PPAN
  derive_monthly_mean_ithkn.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.colors as colors
from yaml import safe_load
import argparse

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
    os.path.join(PPTHN, 'MyPython', 'draw_map'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])
import mod_time as mtime
from mod_utils_fig import bottom_text


parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="Year to plot, default=2025", default=2025, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    default='north', type=str)
args = parser.parse_args()

YR     = args.yr
regn   = args.regn

pthglr = "/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/GLORYS_ithkn_interp_UFSmesh025"
pthcice = "/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/ML_emulator_ithkn/predictors/mean_ithkn"

# GLORYS:
if YR <= 2025:
  fltmp = f"GLORYS_monthly_icevol_ithknmn_{regn}_1993_2025.npz"
else:
  fltmp = f"GLORYS_monthly_icevol_ithknmn_{regn}_2026_2026.npz"

dflout = os.path.join(pthglr, fltmp)
print(f"Loading ice vol and mean ice thickness --> {dflout}")
data = np.load(dflout)
ivolG = data['IVOL']
ithkG = data['ITHKM']
DNMB = data['DNMB']

YRS, MMS, DDS = mtime.datevec(DNMB[0])[:3]
YRE, MME, DDE = mtime.datevec(DNMB[-1])[:3]
years = np.arange(YRS,YRE+1)
iyr = np.where(years == YR)[0][0]

nrec = len(ivolG)
nyr = nrec // 12

# Extract year:
ivol_glr = ivolG.reshape(nyr,12)[iyr,:]
ithk_glr = ithkG.reshape(nyr,12)[iyr,:]

# SOCA CICE6 data:
months = []
ivol_cice = []
ithk_cice = []
for MM in range(1,13):
  flcice = f"cice6_mean_ithkn_{YR}{MM:02d}_north.npz"
  dflcice = os.path.join(pthcice, flcice)
  print(f"Loading ice vol and mean ice thickness --> {dflcice}")

  dcice = np.load(dflcice)
  ivol_cice.append(dcice['IVOL'] * 1e-9)  # m3 --> km3
  ithk_cice.append(dcice['ITHKN'])
  months.append(MM)  

clrg = [0.2, 0.7, 0.3]
clrc = [0.7, 0., 0.3]

plt.ion()

#mm = 7

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.07, 0.55, 0.75, 0.4])

ln1, = ax1.plot(months, ivol_glr, '-o', color=clrg, lw=2, ms=7, label="GLORYS")
ln2, = ax1.plot(months, ivol_cice, '-o', color=clrc, lw=2, ms=7, label="CICE6")

LGD = [ln1, ln2]

ax1.set_title(f'Ice Vol, km3, {YR}')
ax1.set_xticks(months)
ax1.grid('on')

ax2 = plt.axes([0.07, 0.1, 0.75, 0.4])
ax2.plot(months, ithk_glr, '-o', color=clrg,  lw=2, ms=7)
ax2.plot(months, ithk_cice, '-o', color=clrc, lw=2, ms=7)

ax2.set_title(f'Mean ice thikn (over ice-covered area), m, {YR}')
ax2.set_xticks(months)
ax2.grid('on')

ax3 = plt.axes([0.83, 0.1, 0.17, 0.4])
lgd = plt.legend(handles=LGD, loc='lower left')
ax3.axis('off')


btx = 'plot_mnth_icestat_CICE6_GLORYS.py'
bottom_text(btx, pos=[0.02,0.04], fsz=8)








