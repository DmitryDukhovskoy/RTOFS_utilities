"""
  Plot predicted ice thickness
  Saved fields in test_linregr_ithkn.py
  or
  predict_linregr_ithkn_Ndays.py
  or similar script that dumps
  ithkn fields at specified GLORYS grid points
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray as xr
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, cm
import argparse
from pathlib import Path
from yaml import safe_load

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
import mod_glorys as mglr
import mod_icepredict as micepr
from mod_utils_fig import bottom_text
import mod_colormaps as mclrmps

parser = argparse.ArgumentParser()
parser.add_argument("--rdate", help="Prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    default='north', type=str)
parser.add_argument("--model", help="Model name", 
                    choices=['clim','ols1','ols2','rf1','rf2','rf3','gbr1','gbr2','gbr3'],
                    required=True, type=str)
parser.add_argument("--iconc", help="Ice conc field used as a predictor",
                    choices=['glorys','amsr2'],
                    type=str,
                    default="glorys")
args = parser.parse_args()

model   = args.model
rdate    = args.rdate
regn     = args.regn
iconc_fld = args.iconc
syst_info = os.uname()
machine = syst_info.nodename

# Training linregr params:
MODEL_NAMES = micepr.models_info()
model_name = MODEL_NAMES[model]
YS = 1993
YE = 2025

if model == 'ols1':
  YE = 2002

regions = {
    "north": ("Arctic", 65.0),
    "south": ("Antarctic", -60.0),
}
regn_name, lat0 = regions[regn]

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)


# Load prediction and grid points:
pthfcst = os.path.join(config_predictor["linregr"]["pthfcst"],f"{model_name}")
flfcst = f"{model_name}_ithkn_fcast_{rdate}.npz"
if not iconc_fld == "glorys":
  flfcst = f"{model_name}_ithkn_fcast_{iconc_fld}_{rdate}.npz" 

#flfcst = f"{model_name}_{YS}_{YE}_ithkn_fcast_{rdate}.npz"
dflfcst = os.path.join(pthfcst, flfcst)
print(f"Loading fcst {dflfcst}")
data_fcst = np.load(dflfcst)
Ithkn = data_fcst['Yfcst']
JG    = data_fcst['JG']
IG    = data_fcst['IG']

# Get GLORYS grid
dnmb0 = mtime.rdate2datenum(rdate)
YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]

# Find file:
pthice = os.path.join(config_predictor["linregr"]["pthithkn"], f"{YR0}")
dflice = mglr.find_file(rdate, pthice)
assert dflice is not None, f"GLORYS file not found for {rdate} in {pthice}"

with xr.open_dataset(dflice) as dsice:
  A2d = dsice['sithick'].isel(time=0).data.squeeze()
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values

hlon, hlat = np.meshgrid(LON, LAT)

LMsk = None
pthssh = f"/uda/Global_Ocean_Physics_Reanalysis/global/daily/zos/{YR0}"
dflssh = mglr.find_file(rdate, pthssh)
with xr.open_dataset(dflssh) as dszos:
  SSH = dszos['zos'].isel(time=0).values.squeeze()

LMsk = np.where(np.isfinite(SSH),1,0)

# Replace glorys with predicted ithkn
AP = A2d * np.nan
AP[JG,IG] = Ithkn

if LMsk is not None:
  Jocn = (LMsk == 1) & (~np.isfinite(AP))   # open ocean
  AP[Jocn] = 0.

# Plotting
clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.2, 0.2, 0.2])
clrmp.set_under(color=[1,1,1])

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.15, 0.8, 0.8])

#m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l', ax=ax1)
if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l', ax=ax1)
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)

  # Subset region
  JJ = np.where(hlat[:, 0] <= -50)[0]

elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l', ax=ax1)
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)

  # Subset region
  JJ = np.where(hlat[:, 0] >= 50)[0]

hlat_s = hlat[JJ, :]
hlon_s = hlon[JJ, :]
AP_s   = AP[JJ, :]


xh, yh = m(hlon_s, hlat_s)

if regn == 'north':
  m.drawparallels(np.arange(60, 90, 10), labels=[0,0,0,0])
elif regn == 'south':
  m.drawparallels(np.arange(-80, -50, 10), labels=[0,0,0,0])

m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,0])
m.drawcoastlines()

img = m.pcolormesh(xh,yh, AP_s, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(f"ithkn prediction {model_name} {YR0}/{MM0:02d}/{DD0:02d}")

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ticks = np.linspace(rmin, rmax, 11)
clb.set_ticks(ticks)
clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
clb.ax.tick_params(direction='in', length=12)

btx = f' @{machine}: plot_ithkn_prediction.py'
bottom_text(btx, pos=[0.02,0.02], fsz=8)




