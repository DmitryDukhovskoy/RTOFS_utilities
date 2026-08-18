"""
  Calc RMSE of ice thickness from stat model vs GLORYS
  Stat model:

  Predictions should be done using:
  predict_linregr_ithkn_Ndays.py
  or similar script that dumps
  ithkn fields at specified GLORYS grid points

  sdate and edate should be within predicted time
  check saved predictions:
  /work/Dmitry.Dukhovskoy/anls_output/GLORYS_anls/ice_linregr/{model_name}
 
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
from datetime import datetime, timedelta
import matplotlib.dates as mdates

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
from mod_mom6 import dx_dy
import mod_icepredict as micepr
from mod_utils_fig import bottom_text
import mod_colormaps as mclrmps
import calendar

parser = argparse.ArgumentParser()
parser.add_argument("--sdate", help="Start prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--edate", help="End prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--model", help="Model names to analyze",
                    choices=['clim','ols1','ols2','rf1','rf2','rf3','gbr1','gbr2','gbr3'],
                    required=True, 
                    type=str,
                    nargs="+")
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--ys", help="Year start, default=1993", default=1993, type=int)
parser.add_argument("--ye", help="Year end, default=2025", default=2025, type=int)
args = parser.parse_args()

MODELS = args.model
sdate  = args.sdate
edate  = args.edate
regn   = args.regn

ndays_era = 7   # freq. of saved era5 fields


syst_info = os.uname()
machine = syst_info.nodename

MODEL_NAMES = [micepr.models_info()[m] for m in MODELS]

print(f"Models requested: {MODEL_NAMES}")

# Requested start / end time - those may change
# based on saved ERA5 time
dnmbS0 = mtime.rdate2datenum(sdate)
YRS0, MMS0, DDS0 = mtime.datevec(dnmbS0)[:3]
dnmbE0 = mtime.rdate2datenum(edate)
YRE0, MME0, DDE0 = mtime.datevec(dnmbE0)[:3]

regions = {
    "north": ("Arctic", 65.0),
    "south": ("Antarctic", -60.0),
}
regn_name, lat0 = regions[regn]

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)


DIRS = {
  "pthithkn" : config_predictor["linregr"]["pthithkn"],
  "pthiconc" : config_predictor["linregr"]["pthiconc"],
  "pthsst"   : config_predictor["linregr"]["pthsst"],
  "pthssh"   : config_predictor["linregr"]["pthssh"],
  "pthui"    : config_predictor["linregr"]["pthui"],
  "pthvi"    : config_predictor["linregr"]["pthvi"],
  "ptht2m"   : config_predictor["linregr"]["ptht2m"].format(regn_name=regn_name),
  "pthgmapi" : config_predictor["linregr"]["pthgmapi"],
  "pthout"   : config_predictor["linregr"]["pthout"],
  "pthfcst"  : config_predictor["linregr"]["pthout"],
  }


# Construct time array of available ERA5 SAT fields
# within requested time window for prediction ithkn
ptht2m = DIRS["ptht2m"]
DNMB_full = micepr.derive_time(YRS0, YRE0, ptht2m, regn_name, ndays_era)

ixS = np.argmin(abs(DNMB_full - dnmbS0))
ixE = np.argmin(abs(DNMB_full - dnmbE0))

assert abs(DNMB_full[ixS] - dnmbS0) <= ndays_era, f"Check ixS={ixS} in ERA DNMB"
assert abs(DNMB_full[ixE] - dnmbE0) <= ndays_era, f"Check ixE={ixE} in ERA DNMB"

DNMB_fcst = DNMB_full[ixS:ixE+1].astype(int)
nfcst = len(DNMB_fcst)

def read_glorys_grid(dnmb0, config_predictor):
  # Get GLORYS grid
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  rdate = YR0*10000 + MM0*100 + DD0

  # Find file:
  pthice = os.path.join(config_predictor["linregr"]["pthithkn"], f"{YR0}")
  dflice = mglr.find_file(rdate, pthice)
  assert dflice is not None, f"GLORYS file not found for {YR0}/{MM0}/{DD0} in {pthice}"

  with xr.open_dataset(dflice) as dsice:
    #A2d = dsice['sithick'].isel(time=0).data.squeeze()
    LON = dsice['longitude'].values
    LAT = dsice['latitude'].values

  hlon, hlat = np.meshgrid(LON, LAT)
  DX, DY = dx_dy(hlon, hlat)
  Acell = DX*DY

  LMsk = None
  pthssh = f"/uda/Global_Ocean_Physics_Reanalysis/global/daily/zos/{YR0}"
  dflssh = mglr.find_file(rdate, pthssh)
  with xr.open_dataset(dflssh) as dszos:
    SSH = dszos['zos'].isel(time=0).values.squeeze()

  LMsk = np.where(np.isfinite(SSH),1,0)

  return LMsk, Acell, hlon, hlat

def read_prediction(dnmb0, model_name, LMsk):
  YR, MM, DD = mtime.datevec(dnmb0, round_hrs=True)[:3]
  rdate = YR*10000 + MM*100 + DD

  # Load prediction and grid points:
  pthfcst = os.path.join(config_predictor["linregr"]["pthfcst"],f"{model_name}")
  flfcst = f"{model_name}_ithkn_fcast_{rdate}.npz"
  #flfcst = f"{model_name}_{YS}_{YE}_ithkn_fcast_{rdate}.npz"
  dflfcst = os.path.join(pthfcst, flfcst)
  print(f"Loading fcst {dflfcst}")
  data_fcst = np.load(dflfcst, allow_pickle=True)
  Ithkn = data_fcst['Yfcst']
  JG    = data_fcst['JG']
  IG    = data_fcst['IG']

  # Replace glorys with predicted ithkn
  AP = LMsk * np.nan
  AP[JG,IG] = Ithkn

  if LMsk is not None:
    Jocn = (LMsk == 1) & (~np.isfinite(AP))   # open ocean
    AP[Jocn] = 0.

  return AP, IG, JG

def read_glorys(dnmb0, LMsk):
  YR, MM, DD = mtime.datevec(dnmb0, round_hrs=True)[:3]
  pthice = f"/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick/{YR}"

  rdate = YR*10000 + MM*100 + DD
  dflice = mglr.find_file(rdate, pthice)
  if not os.path.isfile(dflice):
    raise RuntimeError(f"File not found: {dflice}")

  with xr.open_dataset(dflice) as dsice:
    A2d = dsice['sithick'].isel(time=0).data.squeeze()
 
  if LMsk is not None:
    Jocn = (LMsk == 1) & (~np.isfinite(A2d))   # open ocean
    A2d[Jocn] = 0.

  return A2d 

def subset_ithkn_climatology(dnmb0, regn, LMsk, HIG):
  pthdata = '/archive/Dmitry.Dukhovskoy/data'
  pthice    = os.path.join(pthdata, 'ithkn_clim_combined')
  fliceout  = 'ithkn_mnthclim_cryo_avhrr_ices_1440x1080_north.nc'
  dfliceout = os.path.join(pthice,fliceout)
  YR, MM, DD = mtime.datevec(dnmb0, round_hrs=True)[:3]
  with xr.open_dataset(dfliceout) as dsice:
    hice = dsice['ice_thkn'].isel(time=MM-1).values.squeeze()

  # Gmapi indices:
  fgmapi = f"gmapi_closenghb_mesh025_to_GLORYS_{regn}.nc"
  pthindx = '/archive/Dmitry.Dukhovskoy/data/remap_indx'
  dflout = os.path.join(pthindx, fgmapi)
  with xr.open_dataset(dflout) as ds:
    IM25 = ds['mesh025_indx'].values
    JM25 = ds['mesh025_jndx'].values
    IGL  = ds['glorys_indx'].values
    JGL  = ds['glorys_jndx'].values

  assert LMsk.shape == HIG.shape, \
    f"LMsk={LMsk.shape}, HIG={HIG.shape}"

  HI = np.zeros_like(LMsk, dtype=float)
  HI[JGL,IGL]  = hice[JM25, IM25]

  # Compare only where GLORYS ice > 0
  # Climatology fields have been extrapolated
  # to the coast in the Arctic 
  HI = np.where(HIG <= 1.e-6, 0., HI)

  return HI

LMsk, Acell, hlon, hlat = read_glorys_grid(DNMB_fcst[0], config_predictor)


nrecs = len(DNMB_fcst)
nmod  = len(MODELS)
RMSE  = np.zeros((nrecs,nmod)) * np.nan
for ictr, model_name in enumerate(MODEL_NAMES): 
  # Load linregr info:
  #pthout = config_predictor["linregr"]["pthout"]
  #flinfo = f"{model_name}_info.npz"
  #dflinfo = os.path.join(pthout, flinfo)
  print(f"Processing {model_name}")

  # GLORYS subset of grid points, where prediction is tested:
  if model_name == 'clim':
    if regn == 'north':
      JG, IG = np.where((LMsk == 1) & (hlat > lat0))
    elif regn == 'south':
      JG, IG = np.where((LMsk == 1) & (hlat < lat0))
  else:
    _, IG, JG = read_prediction(DNMB_fcst[0], model_name, LMsk)

  for irec, dnmb0 in enumerate(DNMB_fcst):
    HIG = read_glorys(dnmb0, LMsk)
    if model_name == 'clim':
      # Read ice thickness CryoSat clim on UFS mesh025
      HIP = subset_ithkn_climatology(dnmb0, regn, LMsk, HIG)
    else:
      HIP, _, _ = read_prediction(dnmb0, model_name, LMsk)

    Atot = np.nansum(Acell[JG,IG])
    ErrSq = np.nansum(Acell[JG,IG] * (HIP[JG,IG] - HIG[JG,IG])**2) / Atot
    print(f"RMSE = {np.sqrt(ErrSq)}")
    RMSE[irec, ictr] = np.sqrt(ErrSq)

    
CLRS = micepr.sens_tests_colors()
YR0, MM0, DD0 = mtime.datevec(DNMB_fcst[0])[:3]
dnmbJ1 = mtime.datenum([YR0,1,1])
fday = (DNMB_fcst - DNMB_fcst[0]).astype(int)

# Create dates for axis 
start_date = datetime(YR0, MM0, DD0)
dates = [start_date + timedelta(days=int(d)) for d in fday]


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.55, 0.8, 0.4])

LGD = []
for imm in range(len(MODEL_NAMES)):
  mdl = MODEL_NAMES[imm]
  clr = CLRS[imm,:]
  ln1, = ax1.plot(dates, RMSE[:,imm], color=clr, label=f"{mdl}", linewidth=2)
  LGD.append(ln1)

# Major ticks every month
ax1.xaxis.set_major_locator(mdates.MonthLocator())

# Labels like Sep/2025, Oct/2025, ...
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
#fig1.autofmt_xdate()   # rotate/align labels if needed

#ax1.set_xticklabels(month_labels)
ax1.set_xlim(min(dates), max(dates))

ax1.set_title(f'RMSE, IceThkn')
ax1.grid('on')


ax3 = plt.axes([0.1, 0.3, 0.25, 0.4])
lgd = plt.legend(handles=LGD, loc='lower left')
ax3.axis('off')


btx = f' @{machine}: calc_rmse_ithkn_prediction.py'
bottom_text(btx, pos=[0.1,0.25], fsz=8)

