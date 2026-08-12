"""
  Train random forest using predictors tested in linregr model

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
import argparse
from yaml import safe_load
import time
from datetime import datetime
import json
import joblib
import statsmodels.api as statm
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

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

#from MyPython.mod_cice6_utils import change_base_template, flname_replace_date

parser = argparse.ArgumentParser()
parser.add_argument("--rfmod", help="random forest model number",
                    choices=[1,2,3,4], required=True, type=int)
parser.add_argument("--regn", help="Region to process, default=north", 
                    choices=['north','south'], 
                    default='north', 
                    type=str)
parser.add_argument("--save", help="Save model object (=1)",
                    choices=[0,1],
                    required=True,
                    type=int)
args = parser.parse_args()

regn  = args.regn
rfmod = args.rfmod
f_save = args.save == 1


# Parameters used in generating predictors
dxy   = 50 
YS    = 1993
YE    = 2025

sqrt_frzdays = True   # use sqrt(integrated freeze days) to better fit Zubov relation
intgr_time = 90  # Time for freeze degree days accumulation, back from current time
Tfrz = -1.85    # ocea freezing T
ndays_era = 7   # freq. of saved era5 fields

regions = {
    "north": ("Arctic", 65.0),
    "south": ("Antarctic", -60.0),
}
regn_name, lat0 = regions[regn]

# Static predictors:
#   yday   - year day represented as cos(yday) + sin(yday)
#   gcoord - geogr. coord. in spherical coordinates
# Dynamic predictors:
#   mnithkn - monthly mean ice thickness (over ice area!), to account for interann. trend
#   iconc   - ice concentration
#   sst     - ocean SST
#   divu    - area-mean ice divergence
#   frzdays - integrated freeze degree days
#   sat     - atm. surf. temp
#
PRED = ["yday", "gcoord", "mnithkn", "iconc", "sst", "divu", "frzdays", "sat"]

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

# Read time array and J,I sample grid points:
# Time array should match ERA5 extracted fields
pthout = DIRS["pthout"]
flithkn = DIRS["ithkntmp"]
dflithkn = os.path.join(pthout, flithkn)

# Sample locations on GLORYS grid and time (date numbers):
print(f"Loading saved {dflithkn}")
assert os.path.isfile(dflithkn), f"Missing tmp file {dflithkn}\n First, create all predictors derive_*py"
data = np.load(dflithkn)
JG   = data["JG"]
IG   = data["IG"]
DNMB = data["DNMB"]

# Read GLORYS grid:
pthice = os.path.join(DIRS["pthithkn"],f"{YS}")
rdate = f"{YS*10000+100+1}"

# Find file:
dflglr = mglr.find_file(rdate, pthice)

with xr.open_dataset(dflglr) as dsice:
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values

# 0 <= lon < 360
LON = (LON + 360) % 360
hlon, hlat = np.meshgrid(LON, LAT)

# Land mask:
LMsk = None
pthssh = os.path.join(DIRS["pthssh"],f"{YS}")
dflssh = mglr.find_file(rdate, pthssh)
with xr.open_dataset(dflssh) as dszos:
  SSH = dszos['zos'].isel(time=0).values.squeeze()

LMsk = np.where(np.isfinite(SSH),1,0)


def load_predictor(predict, DIRS, DNMB):
  """
    These are created in derive_{predictor}.py
  """
  pthout  = DIRS["pthout"]
  fldname = DIRS[predict]
  fltmp   = DIRS[fldname]
  dfltmp  = os.path.join(pthout, fltmp)

  print(f"Loading {predict}: {dfltmp}")
  data = np.load(dfltmp)
  F1d  = data["YY"]
  DNMB_check = data["DNMB"]
  print(f"{predict} N records = {len(DNMB_check)}, expected={len(DNMB)}")
  dtmp = np.floor(np.abs(DNMB - DNMB_check))
  assert np.max(dtmp) == 0, f"{predict}: Check DNMB - dates do not match processed days"

  return F1d

def update_lists_predictor(predict, sqrt_frzdays, PRED_LIST, PRED_STDZ, PRED_MEAN, 
                           PRED_STDEV, PRED_NAMES, DIRS, DNMB):
  FLD2d = load_predictor(predict, DIRS, DNMB)

  if predict == 'frzdays' and sqrt_frzdays == True:
    assert np.min(FLD2d) >= 0, f"Negative freeze degree days found"
    FLD2d = np.sqrt(FLD2d)

  FLD1d = FLD2d.ravel(order='C')
  print(f"{predict} array size={len(FLD1d)}")

  PRED_LIST.append(FLD1d)
  mu = FLD1d.mean()
  stdev = FLD1d.std()

  if stdev == 0:
    raise ValueError(f"Predictor '{predict}' is constant.")

  PRED_STDZ.append((FLD1d - mu) / stdev)

  PRED_MEAN.append(mu)
  PRED_STDEV.append(stdev)

  PRED_NAMES.append(predict)

  return 


def save_random_forest(rf, outdir, model_name="random_forest",
                       predictor_names=None,
                       train_idx=None, test_idx=None,
                       metadata=None):
  """
  Save trained Random Forest model and related information.
  rf : sklearn RandomForestRegressor
      Trained Random Forest estimator.
  outdir : str
      Directory where files will be saved.
  model_name : str
      Prefix for output files.
  predictor_names : list, optional
      Names of predictors in the same order as columns in X.
  train_idx, test_idx : array-like, optional
      Indices of training and testing samples.
  metadata : dict, optional
      Additional information to save (e.g., date range, data source).
  """

  os.makedirs(outdir, exist_ok=True)

  # Save trained model
  model_file = os.path.join(outdir, model_name + ".pkl")
  print(f"Saving trained model ---> {model_file}")
  joblib.dump(rf, model_file)

  # Save estimator parameters
  params_file = os.path.join(outdir, model_name + "_params.json")
  print(f"Saving model parameters ---> {params_file}")
  with open(params_file, "w") as f:
    json.dump(rf.get_params(), f, indent=4, default=str)

  # Save predictor names
  if predictor_names is not None:
    if train_idx is None:
      train_idx = []
    if test_idx is None:
      test_idx = []

    pred_file = os.path.join(outdir, model_name + "_predictors_trainidx.npz")
    print(f"Saving predictors, train/test indx ---> {pred_file}")
    np.savez(
      pred_file,
      PRED_NAMES=PRED_NAMES,
      train_idx=train_idx,
      test_idx=test_idx
      )

  # Save metadata
  info = {
      "saved_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
      "sklearn_model": type(rf).__name__,
  }

  if metadata is not None:
      info.update(metadata)

  info_file = os.path.join(outdir, model_name + "_info.json")
  print(f"Saving model info ---> {info_file}")
  with open(info_file, "w") as f:
      json.dump(info, f, indent=4)


PRED_LIST  = [] 
PRED_STDZ  = []
PRED_MEAN  = []
PRED_STDEV = []
PRED_NAMES = []

for predict in PRED:
  # Construct static predictors:
  if predict == "yday":
    cosD, sinD = micepr.construct_ydays(DNMB, len(IG), order_fast="time") 

    PRED_LIST.extend([cosD, sinD])
    PRED_STDZ.extend([cosD, sinD]) # Do not standardize cosD, sinD
    PRED_MEAN.extend([np.nan, np.nan])
    PRED_STDEV.extend([np.nan, np.nan])
    PRED_NAMES.extend(["cosD", "sinD"])

  elif predict == "gcoord":
    Xcrd, Ycrd, Zcrd = micepr.construct_coord_sphere(hlon, hlat, IG, JG, len(DNMB), order_fast="time")

    # Note: random forest does not need standardized response/ predctor var
    # Standardize only Zcrd because it is >0.9 (mean is not ~0)
    # Xcrd and Ycrd means are close to 0 and within [-1, 1]
    Zstd = (Zcrd - Zcrd.mean()) / Zcrd.std()
    PRED_LIST.extend([Xcrd, Ycrd, Zcrd])
    PRED_STDZ.extend([Xcrd, Ycrd, Zstd])
    PRED_MEAN.extend([np.nan, np.nan, Zcrd.mean()])
    PRED_STDEV.extend([np.nan, np.nan, Zcrd.std()])
    PRED_NAMES.extend(["Xcrd", "Ycrd", "Zcrd"])

  elif predict == "mnithkn":
    pthout = DIRS["pthout"]
    fltmp = f"GLORYS_monthly_icevol_ithknmn_{regn}_{YS}_{YE}.npz"
    dflmni = os.path.join(pthout, fltmp)
    assert os.path.isfile(dflmni), f"File is missing: {dflmni}"

    mnithkn = micepr.construct_mean_ithkn(len(IG), DNMB, dflmni, order_fast="time", Mavrg = 3)
    mnithkn_std = (mnithkn - mnithkn.mean()) / mnithkn.std()
    PRED_LIST.extend([mnithkn])
    PRED_STDZ.extend([mnithkn_std])
    PRED_MEAN.extend([mnithkn.mean()])
    PRED_STDEV.extend([mnithkn.std()])
    PRED_NAMES.extend(["mnithkn"])

  else:
    # Construct dynamic predictors
    # Flatten 2D arrays into 1D ("row-major order (C): row1: col1, ..., colN, row2: col1, ..., colN, ...)
    print(f"Constructing {predict}")
    update_lists_predictor(
        predict,
        sqrt_frzdays,
        PRED_LIST,
        PRED_STDZ,
        PRED_MEAN,
        PRED_STDEV,
        PRED_NAMES,
        DIRS,
        DNMB,
    )

# Response variable:
hice_max = 5.
Ithkn = load_predictor("ithkn", DIRS, DNMB)
Ithkn[Ithkn > hice_max] = hice_max
Y = Ithkn.ravel(order='C')

# Predictor 2D array, predictor matrix or design matrix
# Not standardized - not needed
AA = np.column_stack(PRED_LIST)

# Rand. Forest does not need standardized var.
# Standardized predictors (no intercept !):
#Astdz = np.column_stack(PRED_STDZ)


# Model parameters:
test_sz = 0.2
Ntrees = None    # N of trees
MinLeaf = 5
RandState = 42

if rfmod == 1:
  Ntrees = 100
elif rfmod == 2:
  Ntrees = 200
elif rfmod == 3:
  Ntrees = 5
elif rfmod == 4:
  Ntrees = 2

print("Start training RF")
# Split data
print(f"  Splitting data --> Train ({(1-test_sz)*100:.1f}%) / Test ({test_sz*100:.1f}%)")
INDX = np.arange(len(Y))
# Do no use this as this will split the data mixing grid points from the same
# date into tests/train grid points
#Xtrain, Xtest, Ytrain, Ytest, train_idx, test_idx = train_test_split(
#    AA, Y, INDX, test_size=test_sz, random_state=RandState
#    )

nrec = len(DNMB)
nloc = len(IG)

# Split train/test by eliminating all grid points from the same day entirely
RandGen = np.random.default_rng(RandState)
ndays_test = int(test_sz * nrec)
test_days  = RandGen.choice(nrec, size=ndays_test, replace=False)
train_days = np.setdiff1d(np.arange(nrec), test_days)

# Find indices of flattend 1D data arrays: loc1_time1, loc1_time2, ..., loc2_time1, loc2_time2, ...
# The records in day d: d, d+nrec, d+2*nrec, ..., d+(nloc-1)*nrec
train_idx = (
  np.asarray(train_days)[:, None] +   # col. array of train day indices [[d1],[d2], [d3], ...]
  np.arange(nloc)[None,:]*nrec       # row array of location_indx*nrec
  ).ravel()

test_idx = (
  np.asarray(test_days)[:, None] + 
  np.arange(nloc)[None,:]*nrec
  ).ravel()

Xtrain = AA[train_idx, :]
Ytrain = Y[train_idx]

Xtest = AA[test_idx, :]
Ytest = Y[test_idx]


# Check:
print("N training:", len(train_idx))
print("N testing :", len(test_idx))

print("Unique training days:", len(np.unique(train_days)))
print("Unique testing days :", len(np.unique(test_days)))

common_days = np.intersect1d(
    np.unique(train_days),
    np.unique(test_days)
)

print("Days appearing in BOTH:", len(common_days))


# Create forest
rf = RandomForestRegressor(
    n_estimators=Ntrees,      # number of trees
    max_depth=None,        # grow trees until stopping criteria - no stopping in RF
    min_samples_leaf=MinLeaf,
    n_jobs=-1,             # use all CPUs
    random_state=RandState,
    verbose=2
)

t0 = time.time()
# Train
print("Training ...")
rf.fit(Xtrain, Ytrain)

telaps = time.time() - t0

print(f"{Ntrees} trees elapsed time: {telaps/60:.1f} min")


# Predict
Ypred = rf.predict(Xtest)

# Score
r2 = rf.score(Xtest, Ytest)
print(f"Test prediction:  R2 = {r2:.8f}")

f_modelinfo = False
if f_modelinfo:
  total_nodes = 0
  for i, tree in enumerate(rf.estimators_):
    n = tree.tree_.node_count
    total_nodes += n
    print(f"Tree {i}: {n:,} nodes")

  print(f"Total nodes = {total_nodes:,}")

  depths = [t.tree_.max_depth for t in rf.estimators_]
  print("Mean depth =", np.mean(depths))
  print("Max depth  =", np.max(depths))


metadata={
    "training_yrS": YS,
    "training_yrE": YE,
    "Tfrz": Tfrz,
    "sqrt_frzdays": sqrt_frzdays,
    "ndays_era": ndays_era,
    "intgr_time": intgr_time,
    "dxy": dxy,
    "target": "ice thickness",
    "n_trees": Ntrees,
    "test_size": test_sz,
    "MinLeaf":  MinLeaf,
    "RandState": RandState,
}


if f_save:
  model_name = f"RF_model{rfmod:02d}_{regn}"
  pthrf = DIRS['pthrf']
  print(f"Saving stat model --> {pthrf}")

  save_random_forest(rf, pthrf, model_name=model_name,
                       predictor_names=PRED_NAMES,
                       train_idx=train_idx, test_idx=test_idx,
                       metadata=metadata)

f_plt = False
if f_plt:
  # Check training / testing data
  # Reshape train / test indices: location x time to check 
  # what points / time are used
  nrec = len(DNMB)
  nloc = len(IG)
  ilocTR, itimeTR = np.unravel_index(train_idx, (nloc, nrec))
  ilocTS, itimeTS = np.unravel_index(test_idx, (nloc, nrec))

  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
  ax1.scatter(itimeTR,ilocTR,  s=5, color=(0.5,0.6,0.9), marker='.') # blue - training data points
  ax1.scatter(itimeTS,ilocTS,  s=5, color=(0.9,0.4,0.), marker='.')

  ax1.set_xlabel('Time')
  ax1.set_ylabel('Location')


