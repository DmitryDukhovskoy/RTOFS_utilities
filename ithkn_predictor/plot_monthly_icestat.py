"""
  Plot monthly ice statistics

  derived in:
  derive_monthly_mean_ithkn.py

The daily data is on uda:
/uda/Global_Ocean_Physics_Reanalysis/global/daily/siconc/

/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick

ERA5 T2m fields:
  Downloaded every 7-day daily mean 2m SAT for specified region from ERA5 website
  https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics?tab=download

OR use 1-hr fields on PPAN --> derive daily mean
/archive/uda/ERA5/Hourly_Data_On_Single_Levels/reanalysis/global/1hr-timestep/annual_file-range/Temperature_and_Pressure/T_2m
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray as xr
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, cm
from yaml import safe_load
import argparse
from pathlib import Path

#ROOT = Path(__file__).resolve().parent

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
import mod_icepredict as micepr
from mod_utils_fig import bottom_text

#from MyPython.mod_cice6_utils import change_base_template, flname_replace_date

parser = argparse.ArgumentParser()
parser.add_argument("--ys", help="Year start, default=1993", default=1993, type=int)
parser.add_argument("--ye", help="Year end, default=2025", default=2025, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'], 
                    default='north', type=str)
parser.add_argument("--mm", help="Month to show statistics for",
                    required=True, type=int,
                    nargs="+")
args = parser.parse_args()

YS     = args.ys
YE     = args.ye
regn   = args.regn
MONTHS = args.mm


# Only central Arctic 
regions = {
    "north": ("Arctic", 65.0),
    "south": ("Antarctic", -60.0),
}
regn_name, lat0 = regions[regn]


syst_info = os.uname()
machine = syst_info.nodename

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

DIRS = {
  "pthithkn" : config_predictor["linregr"]["pthithkn"],
  "pthiconc" : config_predictor["linregr"]["pthiconc"],
  "pthsst"   : config_predictor["linregr"]["pthsst"],
  "pthssh"   : config_predictor["linregr"]["pthssh"],
  "ptht2m"   : config_predictor["linregr"]["ptht2m"].format(regn_name=regn_name),
  "pthout"   : config_predictor["linregr"]["pthout"],
  }


pthout = DIRS["pthout"]
fltmp = f"GLORYS_monthly_icevol_ithknmn_{regn}_{YS}_{YE}.npz"
dflout = os.path.join(pthout, fltmp)
print(f"Loading ice vol and mean ice thickness --> {dflout}")
data = np.load(dflout)
IVOL = data['IVOL']
ITHK = data['ITHKM']
DNMB = data['DNMB']

YRS, MMS, DDS = mtime.datevec(DNMB[0])[:3]
YRE, MME, DDE = mtime.datevec(DNMB[-1])[:3]
years = np.arange(YRS,YRE+1)

nrec = len(IVOL)
nyr = nrec // 12

IVOL2D = IVOL.reshape(nyr,12)
ITHK2D = ITHK.reshape(nyr,12)

CLRS = micepr.sens_tests_colors()

plt.ion()

#mm = 7

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.07, 0.55, 0.75, 0.4])

LGD = []
for imm, mm in enumerate(MONTHS):
  clr = CLRS[imm,:]
  ln1, = ax1.plot(years, IVOL2D[:,mm-1], color=clr, label=f"MM={mm:02d}")
  LGD.append(ln1)

ax1.set_title(f'Ice Vol, km3, Months')
ax1.grid('on')

ax2 = plt.axes([0.07, 0.1, 0.75, 0.4])
for imm, mm in enumerate(MONTHS):
  clr = CLRS[imm,:]
  ax2.plot(years, ITHK2D[:,mm-1], color=clr)

ax2.set_title(f'Mean ice thikn (over ice), m, Months')
ax2.grid('on')

ax3 = plt.axes([0.83, 0.1, 0.17, 0.4])
lgd = plt.legend(handles=LGD, loc='lower left')
ax3.axis('off')


btx = f' @{machine}: plot_monthly_icestat.py'
bottom_text(btx, pos=[0.02,0.04], fsz=8)

