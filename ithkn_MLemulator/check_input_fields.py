NOT FINISHED
"""
  Check if all input fields are ready for ML prediction
"""
import numpy as np
import mod_time as mtime
import xarray as xr
import os
import argparse
import sys
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
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])
from mod_utils_fig import bottom_text, colorbar_horiz
import mod_time as mtime
import mod_mom6 as mmom6
import mod_ml_emulator as mml

# See mod_icepredict.py for more info about the models
#    Random Forest (all 1993-2025):
#      RF1:  Ntrees 100
#      RF2:  Ntrees 200
#    Hist. Grad Boost Regressor (decision tree)
#      GBR1: max_leaf = 8,  max_iter=500,  max_depth=8
#      GBR2: max_leaf = 31, max_iter=1000, max_depth=10
#      GBR3: max_leaf = 63, max_iter=1500, max_depth=15
parser = argparse.ArgumentParser()
parser.add_argument("--model", help="ML model to use",
                   choices=['rf1','rf2','gbr2','gbr3'],
                   type=str,
                   required=True)
parser.add_argument("--rdate", help="Start prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
args = parser.parse_args()

ml_model  = args.model
rdate     = args.rdate
regn      = args.regn

fyaml = "paths_ML.yaml"
with open(fyaml) as ff:
  config_ml = safe_load(ff)

# Grid
pthgrid    = config_ml["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

# Requested start / end time - those may change
dnmbS = mtime.rdate2datenum(rdate)
YRS, MMS, DDS = mtime.datevec(dnmbS)[:3]

# Load model parameters and RF model:
model_name = config_ml["ml"][ml_model]["mlname"].format(regn=regn)
pthmodel = config_ml["ml"][ml_model]["pthmodel"]
model_file = os.path.join(pthmodel, model_name + ".pkl")
print(f"Reading {ml_model} object from {model_file}")
mlem = joblib.load(model_file)


# ML Predictors / Input fields used for training:
pred_file = os.path.join(pthmodel, model_name + "_predictors_trainidx.npz")
data_pred = np.load(pred_file)
PRED_NAMES = data_pred["PRED_NAMES"]

if 'Xcrd' in PRED_NAMES or 'Ycrd' in PRED_NAMES or 'Zcrd' in PRED_NAMES:
  if os.paths.isfile(dfgrid_mom):
    print("Xcrd, Ycrd, Zcrd:    ok")
  else:
    print(f"Xcrd, Ycrd, Zcrd:  missing {dfgrid_mom}")

if 'mnithkn' in PRED_NAMES:
  pthglr = MLYAML["PRED"]["pthglr"]
  if YR0 < 2026:
    fltmp = f"GLORYS_monthly_icevol_ithknmn_{regn}_1993_2025.npz"
  elif YR0 == 2026:
    fltmp = f"GLORYS_monthly_icevol_ithknmn_north_{YR0}_{YR0}.npz"

  dflmni = os.path.join(pthglr, fltmp)
  if os.path.isfile(dflmni): 
    print("mnithkn:   ok")
  else
    print(f"mnithkn:  File is missing: {dflmni}")

if 'sst' in PRED_NAMES:
  # SOCA MOM6 rdate, 6hr before the f/cast time for IAU:
  fhr = 6  # IAU time window
  dnmb_mom = dnmb0 - 6/24
  YRm, MMm, DDm, HRm = mtime.datevec(dnmb_mom)[:4]
  mom_rdate = int(YRm)*10000 + int(MMm)*100 + int(DDm)
  pthsst = MLYAML["PRED"]["pthsst"].format(mom_rdate=mom_rdate)
  flsst  = MLYAML["PRED"]["flsst"].format(hr=int(HRm), fhr=fhr)
  dflsst = os.path.join(pthsst, flsst)

  if os.path.isfile(dflsst):
    print('SST:   ok')
  else:
    print(f"SST:  File is missing: {dflsst}")

