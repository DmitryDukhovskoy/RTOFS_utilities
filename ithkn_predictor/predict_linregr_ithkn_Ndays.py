"""
  Run N days predictions
  predictions will be performed from sdate to edate
  using available ERA5 daily SAT fields

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
parser.add_argument("--sdate", help="Start prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--edate", help="End prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--save", help="Save predicted ice thickness in npz (0=no, 1=yes)",
                    choices=[0,1], default=0, type=int)
args = parser.parse_args()

sdate     = args.sdate
edate     = args.edate
regn      = args.regn
save_fcst = args.save == 1

f_plt = False
iconc_min = 0.05  # Discard too low ice conc. predictions
sst_max = 5.      # Discard ice in too warm ocean

# Requested start / end time - those may change
# based on saved ERA5 time
dnmbS0 = mtime.rdate2datenum(sdate)
YRS0, MMS0, DDS0 = mtime.datevec(dnmbS0)[:3]
dnmbE0 = mtime.rdate2datenum(edate)
YRE0, MME0, DDE0 = mtime.datevec(dnmbE0)[:3]

# Training linregr params:
# Model OLS_model1:
model_name = "OLS_model1"

regions = {
    "north": ("Arctic", 65.0),
    "south": ("Antarctic", -60.0),
}
regn_name, lat0 = regions[regn]

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

# Load linregr info:
pthout = config_predictor["linregr"]["pthout"]
flinfo = f"{model_name}_info.npz"
dflinfo = os.path.join(pthout, flinfo)

print(f"Reading: mean. stdev, predict. names --> {dflinfo}")
data = np.load(dflinfo)
PRED_NAMES   = data["PRED_NAMES"]
PRED_MEAN    = data["PRED_MEAN"]
PRED_STDEV   = data["PRED_STDEV"]
Tfrz         = data["Tfrz"].item()
sqrt_frzdays = data["sqrt_frzdays"].item()
intgr_time   = data["intgr_time"].item()
ndays_era    = data["ndays_era"].item() 
dxy          = data["dxy"].item()
YS           = data["YS"].item()  # Year start for model training
YE           = data["YE"].item()  # year end for model training
nparams = len(PRED_NAMES) + 1  # for intersept

# Load regr. results / parameters
flstat = f"{model_name}.pkl"
dflstat = os.path.join(pthout, flstat)
print(f"Reading stat model results {dflstat}")
linregr_results = load_pickle(dflstat)

COEF = linregr_results.params
assert len(COEF) == nparams, f"Expected N parameters {nparams} mismatches COEF {len(COEF)}"

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

def construct_predictors_day(hlon, hlat, IG, JG, dnmb0, DIRS):
  """
    Construct predictors for 1 day forecast
  """
  DNMB = np.asarray([dnmb0])
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  rdate = YR0*10000 + MM0*100 + DD0 
  Xcrd, Ycrd, Zcrd = micepr.construct_coord_sphere(hlon, hlat, IG, JG, len(DNMB), order_fast="time")

  # Time:
  cosD, sinD = micepr.construct_ydays(DNMB, len(IG), order_fast="time")

  # Ice conc
  YR = YR0
  print("\nDeriving GLORYS iconc")
  pthice = os.path.join(DIRS["pthiconc"],f"{YR}")
  dflice = mglr.find_file(rdate, pthice)
  if dflice is None:
      raise FileNotFoundError(f"Not found {dflice}")
  Iconc = micepr.subset_glorys_iconc(dflice, IG, JG)

  # SST
  print("\nDeriving GLORYS sst")
  pthice = os.path.join(DIRS["pthsst"],f"{YR}")
  dflice = mglr.find_file(rdate, pthice)
  if dflice is None:
      raise FileNotFoundError(f"Not found {dflice}")
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
  if sqrt_frzdays:
    frzdays = np.sqrt(frzdays)


  # SAT
  print("\nDeriving GLORYS SAT")
  flt2m = f"era5_2mTemp_daily{ndays_era}day_{regn_name}_{YR}.nc"
  ptht2m = DIRS['ptht2m']
  dflt2m = os.path.join(ptht2m, flt2m)

  SAT = micepr.subset_era_sat(dflt2m, IG, JG, dfgmapi, dnmb0)

  # Eliminate ice in the warm ocean:
  Iconc[SST > sst_max] = 0.


  # Transform and Standardize predictors
  raw = {
      "Zcrd": Zcrd,
      "iconc": Iconc,
      "sst": SST,
      "divu": divU,
      "frzdays": frzdays,
      "sat": SAT,
  }

  # Dict. with standardized arrays:
  std_arr = {}
  for name, mu, sigma in zip(PRED_NAMES, PRED_MEAN, PRED_STDEV):
    if np.isnan(mu):
      continue
    std_arr[name] = (raw[name] - mu) / sigma

  # Construct predictors dictionary where each
  # predictor is linked to the standardized or raw array
  # The order of the predictor should match the order in PRED_NAMES
  # saved by lin. regr. model output
  pred_dict = {
      'cosD'    : cosD,
      'sinD'    : sinD,
      'Xcrd'    : Xcrd,
      'Ycrd'    : Ycrd,
      'Zcrd'    : std_arr['Zcrd'],
      'iconc'   : std_arr['iconc'],
      'sst'     : std_arr['sst'],
      'divu'    : std_arr['divu'],
      'frzdays' : std_arr['frzdays'],
      'sat'     : std_arr['sat'],
  }

  #PRED_STDZ = [pred_dict[p] for p in PRED_NAMES]

  # Combine all standardized predictors into a list:
  PRED_STDZ = []
  for predict in PRED_NAMES:
    PRED_STDZ.append(pred_dict[predict])

  return PRED_STDZ, Iconc, SST


print(f"Start forecasts, N forecasts: {nfcst}")

for dnmb0 in DNMB_fcst:
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  rdate = YR0*10000 + MM0*100 + DD0
  print(f"Day {YR0}/{MM0}/{DD0}")
  PRED_STDZ, Iconc, SST = construct_predictors_day(hlon, hlat, IG, JG, dnmb0, DIRS)
  # Construct design / predictor matrix:
  # Include intercept:
  Astdz = np.column_stack((np.ones(len(JG)), *PRED_STDZ))

  # Retrieve lin. regr coefficients:
  X = linregr_results.params

  assert len(X) == Astdz.shape[1], f"N of coeff does not match A shape {A.shape}"

  # Check means of the predictors - not necess. be ~0
  # because use mean and std from trained data !
  print("Means and stdev of the standzd predictors (except: cosD, sinD, Xcrd, Ycrd)")
  for ipred, (name, x) in enumerate(zip(PRED_NAMES, PRED_STDZ)):
    print(f"x{ipred+1} {name}:   {x.mean()},   {x.std()}")


  # Prediction:
  Yfcst = Astdz @ X

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



