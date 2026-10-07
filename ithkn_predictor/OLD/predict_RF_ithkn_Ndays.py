"""
  Random Forest
  model trained in
  train_RF_ithkn.py

  Run N days predictions
  predictions will be performed from sdate to edate
  using available ERA5 daily SAT fields

  Use every N-days of daily data
  and subsample high-res. GLORYS fields: there is lots of spatial
  correlation, no need to do every grid point

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
import argparse
from yaml import safe_load
import joblib
import json

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

parser = argparse.ArgumentParser()
parser.add_argument("--rfmod", help="random forest model number",
                    choices=[1,2], required=True, type=int)
parser.add_argument("--sdate", help="Start prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--edate", help="End prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--save", help="Save predicted ice thickness in npz (0=no, 1=yes)",
                    choices=[0,1], default=0, type=int)
args = parser.parse_args()

rfmod     = args.rfmod
sdate     = args.sdate
edate     = args.edate
regn      = args.regn
save_fcst = args.save == 1

if not save_fcst:
  print(" \n!!!!   PREDICTIONS WILL NOT BE SAVED !!!\n")

f_plt = False
iconc_min = 0.05  # Discard too low ice conc. predictions
sst_max = 5.      # Discard ice in too warm ocean


# Requested start / end time - those may change
# based on saved ERA5 time
dnmbS0 = mtime.rdate2datenum(sdate)
YRS0, MMS0, DDS0 = mtime.datevec(dnmbS0)[:3]
dnmbE0 = mtime.rdate2datenum(edate)
YRE0, MME0, DDE0 = mtime.datevec(dnmbE0)[:3]

model_name = f"RF_model{rfmod:02d}_{regn}"


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
  "pthrf"    : config_predictor["linregr"]["pthrf"],
  "pthfcst"  : config_predictor["linregr"]["pthfcst"],
  }


# Load model parameters and RF model:
pthrf = config_predictor["linregr"]["pthrf"]
model_file = os.path.join(pthrf, model_name + ".pkl")
print(f"Reading RF object from {model_file}")
rf = joblib.load(model_file)

# Training period:
info_file = os.path.join(pthrf, model_name + "_info.json")
with open(info_file, "r") as f:
  info = json.load(f)

YS           = info["training_yrS"]
YE           = info["training_yrE"]
dxy          = info["dxy"]             # ice length scale
Tfrz         = info["Tfrz"]
sqrt_frzdays = info["sqrt_frzdays"]
intgr_time   = info["intgr_time"]
ndays_era    = info["ndays_era"] 


print(f"Training period: {YS}-{YE}")

# RF Predictors:
pred_file = os.path.join(pthrf, model_name + "_predictors_trainidx.npz")
data_pred = np.load(pred_file)
PRED_NAMES = data_pred["PRED_NAMES"]

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

# GLORYS grid:
# Read GLORYS grid:
pthice = os.path.join(DIRS["pthithkn"],f"{YRS0}")

# Find file:
dflglr = mglr.find_file(sdate, pthice)
assert dflglr is not None, f"GLORYS file not found for {sdate} in {pthice}"

with xr.open_dataset(dflglr) as dsice:
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values

# 0 <= lon < 360
LON = (LON + 360) % 360
hlon, hlat = np.meshgrid(LON, LAT)

# Land mask:
LMsk = None
pthssh = os.path.join(DIRS["pthssh"],f"{YRS0}")
dflssh = mglr.find_file(sdate, pthssh)
with xr.open_dataset(dflssh) as dszos:
  SSH = dszos['zos'].isel(time=0).values.squeeze()

LMsk = np.where(np.isfinite(SSH),1,0)

# Points where ice thickness is predicted:
# Define domain:
if regn == 'north':
  DOMAIN = (hlat > lat0) & (LMsk == 1)
elif regn == 'south':
  DOMAIN = (hlat < lat0) & (LMsk == 1)

JG, IG = np.where(DOMAIN)

standz = False  # No standartization for RF !!!
print(f"Start forecasts, N forecasts: {nfcst}")
for dnmb0 in DNMB_fcst:
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  rdate = YR0*10000 + MM0*100 + DD0
  print(f"Day {YR0}/{MM0}/{DD0}")

  # Derive Predictors
  # hlon, hlat - grid where prediction is done (GLORYS)
  # IG, JG - grid point indices where prediction is done
  # DIRS - dictionary with input / output directories
  # PREAD_NAMES - list of predictors in the right order
  # sqrt_frzdays - True: use sqrt of integrated freezing degree days
  # sst_max - max sst threshold for possible sea ice
  # standz - True: standardize predictors (False for RF)
  # regn - region of prediction
  # YS, YE - training period, start/end years
  # ndays_era - time step in ERA5 atm. fields subsets
  # intgr_time - for freeze degree days, integration period, days
  # dxy - ice length scale, used for calc. ice predictors and gird point subset
  PRED, Iconc, SST = micepr.construct_predictors_day(
           hlon, hlat, IG, JG, dnmb0, DIRS, PRED_NAMES, 
           sqrt_frzdays, sst_max, standz, regn, YS, YE,
           ndays_era, intgr_time, Tfrz, 
           dxy=dxy
           )
  # Construct design / predictor matrix:
  AA = np.column_stack(PRED)

  # Prediction:
  Yfcst = rf.predict(AA)

  #Ierr = np.where(Yfcst < 0)[0]
  Yfcst[Yfcst < 0] = 0
  Yfcst[Iconc < iconc_min] = 0
  Yfcst[SST > sst_max] = 0

  if save_fcst:
    pthdump = os.path.join(DIRS["pthfcst"],f"{model_name}")
    os.makedirs(pthdump, exist_ok=True)
    flfcst = f"{model_name}_ithkn_fcast_{rdate}.npz"
    dflfcst = os.path.join(pthdump, flfcst)
    print(f"Saving fcst --> {dflfcst}")
    np.savez(dflfcst,
             Yfcst=Yfcst,
             JG=JG,
             IG=IG)


if f_plt:

  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
  ax1.contour(LMsk, [0.99], linestyles='solid', colors=[(0.5,0.8,1)])
  #ax1.scatter(IG, JG, s=5, color=(0.8,0.4,0), marker='.')

  # Plot ERA5 SAT at the sample locations:
  
  ir0 = 0
  sc = ax1.scatter(
    IG, JG,
    c=Yfcst,
    cmap='jet',
    s=20,          # marker size
    vmin=0,        # optional color scale limits
    vmax=4
 )

  if regn == 'north':
    yl1 = 1650
    yl2 = 2040

  ax1.set_ylim([yl1,yl2])
  ax1.set_title(f"Predicted ithkn, {YR0}/{MM0:02d}/{DD0:02d}")

  plt.colorbar(sc, ax=ax1, label='ithkn, m')



