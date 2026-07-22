"""
  As a first step in developing ice thickness predictor
  train a linear reg. model to test different predictors
  and check if the taks is feasable, e.g. if it can be treated
  as linear problem

  Test linear regr. derived in train_linregr_ithkn.py

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
from pathlib import Path
from yaml import safe_load
from statsmodels.iolib.smpickle import load_pickle

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
parser.add_argument("--rdate", help="Prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
args = parser.parse_args()

rdate  = args.rdate
regn  = args.regn

dnmb0 = mtime.rdate2datenum(rdate)
YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
DNMB = np.asarray([dnmb0], dtype=int)

# Training linregr params:
# Model OLS_model1:
model_name = "OLS_model1"
dxy = 50  # correlation spatial scale, km - distance btw sampled grd pnts
YS = 1993  # start of the training window
YE = 2002  # end of the training window
intgr_time = 90  # Time for freeze degree days accumulation, back from current time
Tfrz = -1.85    # ocea freezing T
ndays_era = 7   # freq. of saved era5 fields

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
  "ithkntmp" : config_predictor["linregr"]["ithkntmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "iconctmp" : config_predictor["linregr"]["iconctmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "ssttmp"   : config_predictor["linregr"]["ssttmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "divutmp"  : config_predictor["linregr"]["divutmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "sattmp"   : config_predictor["linregr"]["sattmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "dfrztmp"  : config_predictor["linregr"]["dfrztmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "iconc"    : "iconctmp",
  "sst"      : "ssttmp",
  "divu"     : "divutmp",
  "frzdays"  : "dfrztmp",
  "sat"      : "sattmp",
  "ithkn"    : "ithkntmp",
  }

# Load linregr info:
pthout = DIRS["pthout"]
flinfo = f"{model_name}_info.npz"
dflinfo = os.path.join(pthout, flinfo)

print(f"Reading: mean. stdev, predict. names --> {dflinfo}")
data = np.load(dflinfo)
PRED_NAMES = data["PRED_NAMES"]
PRED_MEAN  = data["PRED_MEAN"]
PRED_STDEV = data["PRED_STDEV"]
nparams = len(PRED_NAMES) + 1  # for intersept

# Load regr. results / parameters
flstat = f"{model_name}.pkl"
dflstat = os.path.join(pthout, flstat)
print(f"Reading stat model results {dflstat}")
results = load_pickle(dflstat)

COEF = results.params
assert len(COEF) == nparams, f"Expected N parameters {nparams} mismatches COEF {len(COEF)}"

# Derive predictors

# GLORYS grid:
# Read GLORYS grid:
pthice = os.path.join(DIRS["pthithkn"],f"{YR0}")

# Find file:
dflglr = mglr.find_file(rdate, pthice)
assert dflglr is not None, f"GLORYS file not found for {rdate} in {pthice}"

with xr.open_dataset(dflglr) as dsice:
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values

# 0 <= lon < 360
LON = (LON + 360) % 360
hlon, hlat = np.meshgrid(LON, LAT)

# Land mask:
LMsk = None
pthssh = os.path.join(DIRS["pthssh"],f"{YR0}")
dflssh = mglr.find_file(rdate, pthssh)
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
Xcrd, Ycrd, Zcrd = micepr.construct_coord_sphere(hlon, hlat, IG, JG, len(DNMB), order_fast="time")

# Time:
cosD, sinD = micepr.construct_ydays(DNMB, len(IG), order_fast="time")

# Ice conc
YR = YR0
print("\nDeriving GLORYS iconc")
pthice = os.path.join(DIRS["pthiconc"],f"{YR}")
dflice = mglr.find_file(rdate, pthice)
Iconc = micepr.subset_glorys_iconc(dflice, IG, JG)

# SST
print("\nDeriving GLORYS sst")
pthice = os.path.join(DIRS["pthsst"],f"{YR}")
dflice = mglr.find_file(rdate, pthice)
SST = micepr.subset_glorys_sst(dflice, IG, JG)

# divU ice
print("\nDeriving GLORYS divu ice")
pthu = os.path.join(DIRS["pthui"],f"{YR}")
dfui = mglr.find_file(rdate, pthu)
pthv = os.path.join(DIRS["pthvi"],f"{YR}")
dfvi = mglr.find_file(rdate, pthv)
divU = micepr.subset_glorys_divu(dfui, dfvi, IG, JG, hlon, hlat, dxy)

# Freeze days
print("\nDeriving GLORYS Freeze degree days")
ptht2m = DIRS['ptht2m']
pthgmapi = DIRS["pthgmapi"]
flout = f"gmapi_ERA5_to_GLORYS_{regn}.nc"
dfgmapi = os.path.join(pthgmapi, flout)

frzdays = micepr.subset_era_frzdays(IG, JG, dnmb0, intgr_time, ndays_era, 
                                    ptht2m, dfgmapi, regn, Tfrz=Tfrz)

# SAT
print("\nDeriving GLORYS SAT")
flt2m = f"era5_2mTemp_daily7day_Arctic_{YR}.nc"
ptht2m = DIRS['ptht2m']
dflt2m = os.path.join(ptht2m, flt2m)

SAT = micepr.subset_era_sat(dflt2m, IG, JG, dfgmapi, dnmb0)

# Eliminate ice in the warm ocean:
Iconc[SST>5.] = 0.


pred_dict = {
    'cosD'    : cosD,
    'sinD'    : sinD,
    'Xcrd'    : Xcrd,
    'Ycrd'    : Ycrd,
    'Zcrd'    : Zcrd_stdz,
    'iconc'   : Iconc_stdz,
    'sst'     : SST_stdz,
    'divu'    : divU_stdz,
    'frzdays' : frzdays_stdz,
    'SAT'     : SAT_stdz,
}

# Standardize:
for ik, pr in enumerate(PRED_NAMES): 
  mu   = PREAD_MEAN[ik]
  sgm  = PREAD_STDEV[ik]
  if np.isnan(mu) or np.isnan(sgm)::
    continue

  arr = pred_dict[pr]
  match pr:
    case 'Zcrd':
      arr = (Zcrd - mu) / sgm

    case 'iconc':
      arr = (Iconc - mu) / sgm


# COnstruct design / predictor matrix:
# The order of the predictor should match the order in PRED_NAMES
# saved by lin. regr. model output
PRED_STDZ = []

PRED_STDZ = [pred_dict[p] for p in PRED_NAMES]




f_check = False
if f_check:
  # Land mask:
  LMsk = None
  pthssh = os.path.join(DIRS["pthssh"],f"{YR0}")
  dflssh = mglr.find_file(rdate, pthssh)
  with xr.open_dataset(dflssh) as dszos:
    SSH = dszos['zos'].isel(time=0).values.squeeze()

  LMsk = np.where(np.isfinite(SSH),1,0)

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
    c=Iconc,
    cmap='jet',
    s=20,          # marker size
    vmin=0,        # optional color scale limits
    vmax=1
 )

  plt.colorbar(sc, ax=ax1, label='SAT, degC')



