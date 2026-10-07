"""
Version 2
Arctic / Antarctic regions

  General prediction using ML model
  and any input fields

  Specify ice concentration fields used as an input for ML emulator

  Random Forest or 
  Hist Gradient Boost Regressor (decision trees) predictor

  Run N days predictions
  predictions will be performed from sdate to edate
  using tstep - skip days time for derived daily ERA5 fields

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
importlib.reload(micepr)

# See mod_icepredict.py for more info about the models
#    Random Forest (all 1993-2025):
#   Added predictor: integrated Heat Degree Days (IHDD)
#      RF10: Ntrees 25  Max leaf: 15
#      RF11: Ntrees 50  Max leaf: 5
#      RF12: Ntrees 100 Max leaf: 20
#
#    Hist. Grad Boost Regressor (decision tree)
#   Added predictor: integrated Heat Degree Days (IHDD)
#      GBR10: max_leaf = 50, max_iter=1200, max_depth=10
#      GBR11: max_leaf = 63, max_iter=1200, max_depth=15 
#
parser = argparse.ArgumentParser()
parser.add_argument("--model", help="ML model to use",
                   choices=['rf6','rf7','rf8','rf9','rf10','rf11','rf12','gbr10','gbr11'],
                   required=True)
parser.add_argument("--sdate", help="Start prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--edate", help="End prediction date YYYYMMDD, skip if edate=sdate", 
                    type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--iconc", help="Ice conc field used as a predictor",
                    choices=['glorys','amsr2', 'nsidc'],
                    type=str,
                    required=True)
parser.add_argument("--save", help="Save predicted ice thickness in npz (0=no, 1=yes)",
                    choices=[0,1], default=0, type=int)
args = parser.parse_args()

ml_model  = args.model
iconc_fld = args.iconc
sdate     = args.sdate
edate     = args.edate if args.edate is not None else args.sdate
regn      = args.regn
save_fcst = args.save == 1

if not save_fcst:
  print(" \n!!!!   PREDICTIONS WILL NOT BE SAVED !!!\n")

f_plt = False
emul = "ml"
prdfld = "glorys"
iconc_min = 0.05  # Discard too low ice conc. predictions
sst_max = 5.      # Discard ice in too warm ocean
# No standartization for RF or grad. boosting
standz = False


# Requested start / end time - those may change
# based on saved ERA5 time
dnmbS0 = mtime.rdate2datenum(sdate)
YRS0, MMS0, DDS0 = mtime.datevec(dnmbS0)[:3]
dnmbE0 = mtime.rdate2datenum(edate)
YRE0, MME0, DDE0 = mtime.datevec(dnmbE0)[:3]

# Parameters used in generating predictors
fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

# Load parameters:
regn_name    = config_predictor["regn"][regn]["name"]
lat0         = config_predictor["regn"][regn]["lat_bnd"]
davrg        = config_predictor["params"]["davrg"]
tstep_era    = config_predictor["params"]["tstep"]     # ERA5 daily fields, skipping days

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)


model_name = micepr.models_info()[ml_model]

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


# Load model parameters and ML model:
pthmodel = config_predictor[emul][ml_model]["pthmodel"]
model_file = os.path.join(pthmodel, model_name + ".pkl")
print(f"Reading {ml_model} object from {model_file}")
mlem = joblib.load(model_file)

# Training period:
info_file = os.path.join(pthmodel, model_name + "_info.json")
with open(info_file, "r") as f:
  info = json.load(f)

YS           = info["training_yrS"]
YE           = info["training_yrE"]
dxy          = info["dxy"]             # ice length scale
Tfrz         = info["Tfrz"]
sqrt_frzdays = info["sqrt_frzdays"]
intgr_time   = info["intgr_time"]
#tstep_era   = info["tstep_era"] 

print(f"Training period: {YS}-{YE}")

DIRS = {
  "pthithkn"  : config_predictor[emul][prdfld]["pthithkn"],
  "pthsst"    : config_predictor[emul][prdfld]["pthsst"],
  "pthssh"    : config_predictor[emul][prdfld]["pthssh"],
  "pthui"     : config_predictor[emul][prdfld]["pthui"],
  "pthvi"     : config_predictor[emul][prdfld]["pthvi"],
  "ptht2m"    : config_predictor["linregr"]["ptht2m_1hr"],
  "pthgmapi"  : config_predictor[emul][prdfld]["pthgmapi"],
  "pthfcst"   : config_predictor[emul][prdfld]["pthfcst"],
  "pthout"    : config_predictor[emul][prdfld]["pthout"],
  "pthiconc"  : config_predictor[emul][iconc_fld]["pthiconc"],
  "pthmodel"  : config_predictor[emul][ml_model]["pthmodel"],
  "regn"      : regn,
  "dxy"       : dxy,
  "YS"        : YS,
  "YE"        : YE,
  "Tfrz"      : Tfrz,
  "sst_max"   : sst_max,
  "tstep_era" : tstep_era,
  "iconc_fld" : iconc_fld,
  "standz"    : standz,
  "sqrt_frzdays" : sqrt_frzdays,
  "intgr_time"   : intgr_time,
  }

# ML Predictors:
pred_file = os.path.join(pthmodel, model_name + "_predictors_trainidx.npz")
data_pred = np.load(pred_file)
PRED_NAMES = data_pred["PRED_NAMES"]

# Construct time array of ERA5 SAT fields with tstep_era skipping
# within requested time window for prediction ithkn
DNMB_full = micepr.derive_time(YS, YE, tstep_era)

ixS = np.argmin(abs(DNMB_full - dnmbS0))
ixE = np.argmin(abs(DNMB_full - dnmbE0))

assert abs(DNMB_full[ixS] - dnmbS0) <= tstep_era, f"Check ixS={ixS} in ERA DNMB"
assert abs(DNMB_full[ixE] - dnmbE0) <= tstep_era, f"Check ixE={ixE} in ERA DNMB"

DNMB_fcst = DNMB_full[ixS:ixE+1].astype(int)
nfcst = len(DNMB_fcst)

# Points where ice thickness is predicted:
# Define domain:
if regn == 'north':
  DOMAIN = (hlat > lat0) & (LMsk == 1)
elif regn == 'south':
  DOMAIN = (hlat < lat0) & (LMsk == 1)

JG, IG = np.where(DOMAIN)


print(f"Start forecasts, N forecasts: {nfcst}")
for dnmb0 in DNMB_fcst:
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  rdate = YR0*10000 + MM0*100 + DD0
  print(f"Day {YR0}/{MM0}/{DD0}")

  if emul == "ml":
    DIRS["fliconc"] = None # Search for the appropriate file in the directory
  else:
    DIRS["fliconc"] = fliconc.format(
        YR=f"{YR0}",
        MM=f"{MM0:02d}",
        DD=f"{DD0:02d}"
    ) 
 
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
  # tstep - time step in ERA5 atm. fields subsets
  # intgr_time - for freeze degree days, integration period, days
  # dxy - ice length scale, used for calc. ice predictors and gird point subset
  PRED, Iconc, SST = micepr.construct_predictors_day(
           hlon, hlat, IG, JG, dnmb0, DIRS, PRED_NAMES 
           )
  # Construct design / predictor matrix:
  AA = np.column_stack(PRED)

  # Prediction:
  Yfcst = mlem.predict(AA)

  #Ierr = np.where(Yfcst < 0)[0]
  Yfcst[Yfcst < 0] = 0
  Yfcst[Iconc < iconc_min] = 0
  Yfcst[SST > sst_max] = 0

  if save_fcst:
    pthdump = os.path.join(DIRS["pthfcst"],f"{model_name}")
    os.makedirs(pthdump, exist_ok=True)
    flfcst = f"{model_name}_ithkn_fcast_{iconc_fld}_{rdate}.npz"
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



