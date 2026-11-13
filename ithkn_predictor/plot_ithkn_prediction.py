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
parser.add_argument(
  "--regn", 
  help="Region to process",
  choices=['north','south'],
  default='north'
)
parser.add_argument(
  "--model", 
  help="ML model, e.g.: ols1, rf11, gbr20, ... error message if prediction does not exist", 
  required=True
)
parser.add_argument(
  "--iconc", help="Ice conc field used as a predictor",
  choices=['glorys','amsr2','nsidc'],
  default="glorys"
)
args = parser.parse_args()

ml_model  = args.model
rdate     = args.rdate
regn      = args.regn
iconc_fld = args.iconc
syst_info = os.uname()
machine   = syst_info.nodename

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)


# Split model into model short name and number:
ml_model = ml_model.lower()

if ml_model.startswith("gbr") and ml_model[3:].isdigit():
  ml_short = "gbr"
  ml_nmb = int(ml_model[3:])
elif ml_model.startswith("rf") and ml_model[2:].isdigit():
  ml_short = "rf"
  ml_nmb = int(ml_model[2:])
else:
  raise ValueError(f"Model name / number cannot be identified: {ml_model}")

# Check if specified model number exists:
rf_specs = config_predictor['ml_specs'][ml_short]

# Keep only numeric model IDs, exclude pthmodel and other settings
valid_models = [
  int(key) for key in rf_specs
  if key.isdigit()
]

if ml_nmb not in valid_models:
  raise ValueError(
    f"Invalid model '{ml_short}' number: {ml_nmb},  "
    f"Choose one of: {valid_models}"
  )

# Training linregr params:
ml_full = micepr.full_model_name(ml_model)
model_name = f"{ml_full}_{regn}"

pthinfo = config_predictor["params"]["pthinfo"]
flij   = config_predictor["params"]["flij"].format(regn=regn)
flgrid = config_predictor["params"]["flgrid"]
dflij   = os.path.join(pthinfo, flij)
dflgrid = os.path.join(pthinfo, flgrid)

# Load prediction and grid points:
pthfcst = os.path.join(config_predictor["linregr"]["pthfcst"],f"{model_name}")
#flfcst = f"{model_name}_ithkn_fcast_{rdate}.npz"
flfcst = f"{model_name}_ithkn_fcast_{iconc_fld}_{rdate}.npz" 

#flfcst = f"{model_name}_{YS}_{YE}_ithkn_fcast_{rdate}.npz"
dflfcst = os.path.join(pthfcst, flfcst)
print(f"Loading fcst {dflfcst}")
data_fcst = np.load(dflfcst)
Ithkn = data_fcst['Yfcst']
JG    = data_fcst['JG']
IG    = data_fcst['IG']

# GLORYS grid:
A = np.load(dflgrid)
hlon = A["LON"]
hlat = A["LAT"]
LMsk = A["LMsk"]

# Replace glorys with predicted ithkn
AP = LMsk * np.nan
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


dnmb = mtime.rdate2datenum(rdate)
YR0, MM0, DD0 = mtime.datevec(dnmb)[:3]

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
ax1.set_title(f"ithkn {model_name} iconc={iconc_fld}\n{YR0}/{MM0:02d}/{DD0:02d}")

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ticks = np.linspace(rmin, rmax, 11)
clb.set_ticks(ticks)
clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
clb.ax.tick_params(direction='in', length=12)

btx = f' @{machine}: plot_ithkn_prediction.py'
bottom_text(btx, pos=[0.02,0.02], fsz=8)




