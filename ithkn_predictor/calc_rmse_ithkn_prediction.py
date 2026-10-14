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

# Choices: clim - plot climatology prediction
#          ML predictions: check existing forcasts
#          rf[1,....21] - Random Forest
#          gbr[1,...21] - Gradient Boosting Decision Tree
parser = argparse.ArgumentParser()
parser.add_argument("--sdate", help="Start prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--edate", help="End prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--model", help="Model names to analyze",
                    required=True, 
                    nargs="+")
parser.add_argument("--iconc", help="Ice conc field used as a predictor",
                    choices=['glorys','amsr2', 'nsidc'],
                    default='glorys')
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
args = parser.parse_args()

MODELS = args.model
sdate  = args.sdate
edate  = args.edate
regn   = args.regn
iconc_fld = args.iconc

syst_info = os.uname()
machine = syst_info.nodename

#MODEL_NAMES = [micepr.models_info()[m] for m in MODELS]

print(f"Models requested: {MODELS}")

# Requested start / end time - those may change
# based on saved ERA5 time
dnmbS0 = mtime.rdate2datenum(sdate)
YRS0, MMS0, DDS0 = mtime.datevec(dnmbS0)[:3]
dnmbE0 = mtime.rdate2datenum(edate)
YRE0, MME0, DDE0 = mtime.datevec(dnmbE0)[:3]

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

# Load parameters:
regn_name    = config_predictor["regn"][regn]["name"]
lat0         = config_predictor["regn"][regn]["lat_bnd"]
davrg        = config_predictor["params"]["davrg"]
tstep_era    = config_predictor["params"]["tstep"]     # ERA5 daily fields, skipping days

# Get time steps, grid points:
pthinfo = config_predictor["params"]["pthinfo"]
fltime = config_predictor["params"]["fltime"].format(regn=regn, tstep=tstep_era)
flij   = config_predictor["params"]["flij"].format(regn=regn)
flgrid = config_predictor["params"]["flgrid"]
dfltime = os.path.join(pthinfo, fltime)
dflij   = os.path.join(pthinfo, flij)
dflgrid = os.path.join(pthinfo, flgrid)

assert os.path.isfile(dfltime), f"Time steps file is missing: {dfltime}"
assert os.path.isfile(dflij), f"Subsample grid points file is missing: {dflij}"

# GLORYS grid:
A = np.load(dflgrid)
hlon = A["LON"]
hlat = A["LAT"]
LMsk = A["LMsk"]

# 0 <= lon < 360
hlon = (hlon + 360) % 360

DX, DY = dx_dy(hlon, hlat)
Acell = DX*DY

DIRS = {
  "pthithkn" : config_predictor["linregr"]["pthithkn"],
  "pthiconc" : config_predictor["linregr"]["pthiconc"],
  "pthsst"   : config_predictor["linregr"]["pthsst"],
  "pthssh"   : config_predictor["linregr"]["pthssh"],
  "pthui"    : config_predictor["linregr"]["pthui"],
  "pthvi"    : config_predictor["linregr"]["pthvi"],
  "ptht2m"   : config_predictor["linregr"]["ptht2m_1hr"],
  "pthgmapi" : config_predictor["linregr"]["pthgmapi"],
  "pthout"   : config_predictor["linregr"]["pthout"],
  "pthfcst"  : config_predictor["ml"]["pthfcst"],
  }

def read_prediction(dnmb0, model_name, LMsk, pthfcst, flbase):
  YR, MM, DD = mtime.datevec(dnmb0, round_hrs=True)[:3]
  rdate = YR*10000 + MM*100 + DD

  #flfcst = f"{model_name}_ithkn_fcast_{iconc_fld}_{rdate}.npz"
  flfcst = f"{flbase}{rdate}.npz"
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

  if regn == 'north':
    fliceout  = 'ithkn_mnthclim_cryo_avhrr_ices_1440x1080_north.nc'
  elif regn == 'south':
    fliceout = 'CryoSat_hice_mnthclim_2011_2020_mesh025_1440x1080_south.nc'

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

def parse_model_name(ml_model):
  if ml_model.startswith("gbr") and ml_model[3:].isdigit():
    ml_short = "gbr"
    ml_nmb = int(ml_model[3:])
  elif ml_model.startswith("rf") and ml_model[2:].isdigit():
    ml_short = "rf"
    ml_nmb = int(ml_model[2:])
  elif ml_model == "clim":
    ml_short = "clim"
    ml_nmb = None
  else:
    raise ValueError(f"Model name / number cannot be identified: {ml_model}")

  return ml_short, ml_nmb

def construct_time_fcasts(pthdump, flbase, sdate, edate):
  """
  Find start - end days 
  Requested time interval is truncated to existing f/cast time in the f/cast dir
  """
  path_fcst = Path(pthdump)
  files = sorted(path_fcst.glob(f"{flbase}*.npz"))   
  if not files:
    return None

  # Extract dates from f/cast filenames
  RDATES = [f.stem.split("_")[-1] for f in files]
  fcast_dnmbs = mtime.rdate2datenum(RDATES)

  # Requested time interval
  dnmbS = mtime.rdate2datenum(sdate)
  dnmbE = mtime.rdate2datenum(edate)
  assert dnmbE >= dnmbS, f"Start date {sdate} is after end date {edate}"

  # Available forecast time range
  dnmb_min = np.min(fcast_dnmbs)
  dnmb_max = np.max(fcast_dnmbs)

  # Ensure existing f/casts are within the requested time
  yr1, mm1, dd1 = mtime.datevec(dnmb_min)[:3]
  yr2, mm2, dd2 = mtime.datevec(dnmb_max)[:3]

  # Check if the requested time interval overlaps with the existing f/cast:
  if dnmbE < dnmb_min or dnmbS > dnmb_max:
    raise ValueError(
      f"Available forecasts {int(yr1):04d}{int(mm1):02d}{int(dd1):02d} - "
      f"{int(yr2):04d}{int(mm2):02d}{int(dd2):02d}\n"
      f"do not cover the requested interval {sdate}-{edate}"
    )

  # Find indices closest to the requested dates
  idxS = max(np.argmin(abs(fcast_dnmbs - dnmbS)), 0)
  idxE = np.argmin(abs(fcast_dnmbs - dnmbE))

  return fcast_dnmbs[idxS:idxE+1]

nmod  = len(MODELS)
RMSE  = []
FTIME = []
for ictr, ml_model in enumerate(MODELS): 
  if ml_model == "clim":
    # process clim after ML models
    continue

  ml_short, ml_nmb = parse_model_name(ml_model)
  ml_name = micepr.full_model_name(ml_model)
  model_name = f"{ml_name}_{regn}"
  print(f"Processing {model_name}")

  # Construct time array for existing f/casts within specified time period
  pthdump = os.path.join(DIRS["pthfcst"],f"{model_name}")
  flbase = flfcst = f"{model_name}_ithkn_fcast_{iconc_fld}_"
  
  DNMB_fcst = construct_time_fcasts(pthdump, flbase, sdate, edate)
  if DNMB_fcst is None:
    raise RuntimeError(f"No forecast found in {pthdump}")

  # Keep rmse and time as list to allow different receord lengths
  nrecs = len(DNMB_fcst)
  rmse_run  = np.zeros((nrecs))

  # Check the lenght of the f/cast time array, should be the same for all f/casts:
  #assert len(DNMB_fcst) == nrecs, f"Check f/cast time, expected nrecs={nrecs}"

  # Subset of grid points, where prediction is performed
  _, IG, JG = read_prediction(DNMB_fcst[0], model_name, LMsk, pthdump, flbase)
  Atot = np.nansum(Acell[JG, IG])
  if Atot <= 0:
    raise ValueError("Total grid-cell area must be positive")

  for irec, dnmb0 in enumerate(DNMB_fcst):
    HIG = read_glorys(dnmb0, LMsk)
    HIP, _, _ = read_prediction(dnmb0, model_name, LMsk, pthdump, flbase)

    ErrSq = np.nansum(Acell[JG,IG] * (HIP[JG,IG] - HIG[JG,IG])**2) / Atot
    print(f"RMSE = {np.sqrt(ErrSq)}")
    rmse_run[irec] = np.sqrt(ErrSq)

  RMSE.append(rmse_run)
  FTIME.append(DNMB_fcst)

# Process climatology f/cast:
if 'clim' in MODELS:
  ictr = MODELS.index('clim')
  print("Processing climatology")
  if regn == 'north':
    JG, IG = np.where((LMsk == 1) & (hlat > lat0))
  elif regn == 'south':
    JG, IG = np.where((LMsk == 1) & (hlat < lat0))

  nrecs = len(DNMB_fcst)
  rmse_run  = np.zeros((nrecs))
  for irec, dnmb0 in enumerate(DNMB_fcst):
    HIG = read_glorys(dnmb0, LMsk)

    # Read ice thickness CryoSat clim on UFS mesh025
    HIP = subset_ithkn_climatology(dnmb0, regn, LMsk, HIG)

    Atot = np.nansum(Acell[JG,IG])
    ErrSq = np.nansum(Acell[JG,IG] * (HIP[JG,IG] - HIG[JG,IG])**2) / Atot
    print(f"RMSE = {np.sqrt(ErrSq)}")
    rmse_run[irec] = np.sqrt(ErrSq)

  RMSE.append(rmse_run)
  FTIME.append(DNMB_fcst)
    
CLRS = micepr.sens_tests_colors()

def create_time_array(DNMB_fcst):
  YR0, MM0, DD0 = mtime.datevec(DNMB_fcst[0])[:3]
  fday = (DNMB_fcst - DNMB_fcst[0]).astype(int)

  # Create dates for axis 
  start_date = datetime(YR0, MM0, DD0)
  dates = [start_date + timedelta(days=int(d)) for d in fday]

  return dates

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.55, 0.8, 0.4])

LGD = []
for imm in range(len(MODELS)):
  mdl = MODELS[imm]
  clr = CLRS[imm,:]

  rmse_run = RMSE[imm]
  dates = create_time_array(FTIME[imm])

  ln1, = ax1.plot(dates, rmse_run, color=clr, label=f"{mdl}", linewidth=2)
  LGD.append(ln1)

# Major ticks every month
ax1.xaxis.set_major_locator(mdates.MonthLocator())

# Labels like Sep/2025, Oct/2025, ...
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
#fig1.autofmt_xdate()   # rotate/align labels if needed

#ax1.set_xticklabels(month_labels)
ax1.set_xlim(min(dates), max(dates))

ax1.set_title(f'RMSE, IceThkn, {regn}')
ax1.grid('on')


ax3 = plt.axes([0.1, 0.3, 0.25, 0.4])
lgd = plt.legend(handles=LGD, loc='lower left')
ax3.axis('off')


btx = f' @{machine}: calc_rmse_ithkn_prediction.py'
bottom_text(btx, pos=[0.1,0.25], fsz=8)

