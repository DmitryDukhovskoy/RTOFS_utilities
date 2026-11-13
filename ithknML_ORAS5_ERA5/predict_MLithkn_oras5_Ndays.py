"""
General prediction using ML model
and any input fields

Input fields (predictors) are from ORAS5 and ERA5

Specify ice concentration fields used as an input for ML emulator

Hist Gradient Boost Regressor (decision trees) predictor

Run N days predictions
predictions will be performed from sdate to edate
using ORAS5 and ERA5 fields

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
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])
from mod_utils_fig import bottom_text, colorbar_horiz
import mod_time as mtime
import mod_mom6 as mmom6
import mod_ml_emulator as mml
importlib.reload(mml)

def write_MLithkn_nc(dfliceout, time_dnmb, A3d, regn, nctitle):
  yr1, mm1, dd1 = mtime.datevec(time_dnmb[0])[:3]
  nrecs, jdim, idim = A3d.shape

  # Days wrt to the Jan 1st:
  time_days = time_dnmb - mtime.datenum([yr1, 1, 1])
  darr_cice = xr.DataArray(A3d, dims=("time","jdim","idim"),\
                     coords={"time": time_days,\
                             "jdim": np.arange(jdim),\
                             "idim": np.arange(idim)})

  dset = xr.Dataset({"ice_thkn": darr_cice})
  dset['ice_thkn'].attrs['long_name'] = 'ice thickness'
  dset["time"].attrs = {
       "long_name": f"days since {yr1}/01/01",
       "units" : "meters",
  }

  # Add global attributes:
  dset.attrs['title']    = nctitle
  dset.attrs['source']   = "predict_MLithkn_oras5_Ndays.py"
  dset.attrs['region']   = regn

  print(f'Dumping ML ice thickness  --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')

def write_iconc_nc(dfliceout, time_dnmb, A3d, regn, nctitle):
  yr1, mm1, dd1 = mtime.datevec(time_dnmb[0])[:3]
  nrecs, jdim, idim = A3d.shape

  # Days wrt to the Jan 1st:
  time_days = time_dnmb - mtime.datenum([yr1, 1, 1])
  darr_cice = xr.DataArray(A3d, dims=("time","jdim","idim"),\
                     coords={"time": time_days,\
                             "jdim": np.arange(jdim),\
                             "idim": np.arange(idim)})

  dset = xr.Dataset({"ice_conc": darr_cice})
  dset['ice_conc'].attrs['long_name'] = 'ice partial area / fraction'
  dset["time"].attrs = {
       "long_name": f"days since {yr1}/01/01",
       "units" : "meters",
  }

  # Add global attributes:
  dset.attrs['title']    = nctitle
  dset.attrs['source']   = "predict_MLithkn_oras5_Ndays.py"
  dset.attrs['region']   = regn

  print(f'Dumping input ice concentration  --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')

# ====
# See config_ML.yaml for more info about the models
#    Random Forest (all 1993-2025):
#      rf21: btmtopo, IHDD + other predictors, 100 trees, 20 MinLeaf
#    Hist. Grad Boost Regressor (decision tree)
#      gbr21: btmtopo, IHDD + other predictors, max_leaf = 63, mad_depth=15
#
parser = argparse.ArgumentParser()
parser.add_argument(
  "--model",
  help="ML model to use, see config_ML.yaml for options, e.g. gbr21",
  required=True
)
parser.add_argument("--sdate", help="Start prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--edate", help="End prediction date YYYYMMDD, skip if edate=sdate", type=int)
parser.add_argument(
  "--regn", 
  help="Region to process", 
  choices=['north','south'],
  required=True
)
parser.add_argument(
  "--save",
  help="Save predicted ithkn: 0=no, 1=npz (prediction grid pnts), 2=netcdf (whole domain, default))",
  choices=[0, 1, 2],
  default=2,
  type=int
)
parser.add_argument(
  "--sinput", 
  help="=1: Save created input fields for inspection, =0: no (default)",
  choices=[0,1], 
  default=0,
  type=int
)
parser.add_argument(
  "--saveiconc",
  help="=1: Save input iceonc in netcdf, =0: no (default)",
  choices=[0, 1],
  default=0,
  type=int
)

args = parser.parse_args()

ml_model  = args.model
sdate     = args.sdate
edate     = args.edate if args.edate is not None else args.sdate
regn      = args.regn
save_input = args.sinput == 1    # save created input (predictor) fields
save_npz   = args.save == 1      # save numpy binary
save_nc    = args.save == 2      # save netcdf
save_iconc = args.saveiconc == 1

if save_npz:
  print("Output will be saved to numpy binary npz")
elif save_nc:
  print("Output will be saved to netcdf")
else:
  print(" \n!!!!   PREDICTIONS WILL NOT BE SAVED !!!\n")

f_plt = False
standz = False    # standardized predictors only for lin. regr
iconc_min = 0.05  # Discard too low ice conc. predictions
sst_max = 5.      # Discard ice in too warm ocean

fyaml = "config_MLithkn_oras5.yaml"
with open(fyaml) as ff:
  config_ml = safe_load(ff)

# Define ML model to use, check if it has been trained:
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
rf_specs = config_ml['ml_specs'][ml_short]

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

# Model params
# Training period:
YS        = config_ml["params"]["ys"]
YE        = config_ml["params"]["ye"]
lat0      = config_ml["regn"][regn]["lat_bnd"]
dxy       = config_ml["params"]["dxy"]             # ice length scale
divavrg   = config_ml["params"]["davrg"]           # divU time averaging, days
Tfrz      = config_ml["params"]["Tfrz"]
sqrt_frz  = config_ml["params"]["sqrt_frz"]
intgr_frz = config_ml["params"]["intgr_frz"]
Hbtm_min  = config_ml["params"]["Hbtm_min"]   # for btmtopo predictor, normalization

# Requested start / end time 
dnmbS = mtime.rdate2datenum(sdate)
YRS, MMS, DDS = mtime.datevec(dnmbS)[:3]
dnmbE = mtime.rdate2datenum(edate)
YRE, MME, DDE = mtime.datevec(dnmbE)[:3]

# Load model parameters and RF model:
ml_long = mml.full_model_name(ml_model)
model_name = f"{ml_long}_{regn}"

pthmodel = config_ml["ml_specs"][ml_short]["pthmodel"]

model_file = os.path.join(pthmodel, model_name + ".pkl")
print(f"Reading {ml_model} object from {model_file}")
mlem = joblib.load(model_file)

# Training period:
info_file = os.path.join(pthmodel, model_name + "_info.json")
with open(info_file, "r") as f:
  info = json.load(f)

# check that the model is the right one
assert YS == info["training_yrS"], \
  f"Training start year YS={YS} does not match saved value {info['training_yrS']} in {info_file}"
assert YE == info["training_yrE"], \
  f"Training end year YE={YE} does not match saved value {info['training_yrE']} in {info_file}"
assert np.isclose(dxy, info["dxy"]), \
  f"Length scale dxy={dxy:.2f} does not match saved value {info['dxy']:.2f} in {info_file}"
assert Tfrz == info["Tfrz"], \
  f"Tfrz={dxy:.3f} does not match saved value {info['Tfrz']:.3f} in {info_file}"
assert sqrt_frz == info["sqrt_frzdays"],\
  f"Square root of IFDD = {sqrt_frzdays} does not match saved value {info['sqrt_frzdays']} in {info_file}"
assert intgr_frz  == info["intgr_time"],\
  f"IFDD/IHDD integr. time={intgr_frz} does not match saved value {info['intgr_time']} in {info_file}"

# ML Predictors / Input fields used for training:
pred_file = os.path.join(pthmodel, model_name + "_predictors_trainidx.npz")
data_pred = np.load(pred_file)
PRED_NAMES = data_pred["PRED_NAMES"]
npred = len(PRED_NAMES)

print("Predictor Names:")
for f in PRED_NAMES:
  print(f"  {f}")

# Prediction grid points:
# Get MOM6 grid:
pthgrid    = config_ml["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xr.open_dataset(dftopo_mom) as dstopo:
  depth = dstopo['depth'].data.squeeze()

# Convert all positive values -> land (100) and ocean (<0):
HH = np.where(depth < 1.e-20, 100., -depth)

jdm, idm = HH.shape
LMsk = HH < 0

# Set mask for region and ML lat bounds:
assert HH.shape == hlat.shape, f"Shape mismatch: HH={HH.shape}, hlat={hlat.shape}"

if regn == 'north':
  DOMAIN = (hlat > lat0) & (LMsk)
elif regn == 'south':
  DOMAIN = (hlat < lat0) & (LMsk)

Jsub, Isub = np.where(DOMAIN)

# Find gmapi indices for mapping ORAS5 --> mesh025:
# Mapping index pairs mesh025 - ORAS5:
pthgmapi = config_ml["ORAS5"]["pthgmapi"]
flgmapi  = config_ml["ORAS5"]["flgmapi"]
dflgmapi = os.path.join(pthgmapi, flgmapi)
with xr.open_dataset(dflgmapi) as ds:
  LONOR = ds["oras_longit"].values
  LATOR = ds["oras_latit"].values
  IORAS = ds["oras_indx"].values.astype(int)
  JORAS = ds["oras_jndx"].values.astype(int)
  I025  = ds["mesh025_indx"].values.astype(int)
  J025  = ds["mesh025_jndx"].values.astype(int)

LONOR = (LONOR + 360) % 360

# Create a look-up table for matching mesh025 -> ORAS5 index pairs
# Subsample mesh025 grid points to be processed
JORsub = np.empty(len(Jsub), dtype=int)
IORsub = np.empty(len(Isub), dtype=int)

# Build a lookup table for fast finding corresponding ORAS5 indices
print("Building ORAS5 - mesh025 remapping lookup table ...")
gmapi = {(jg,ig):(je,ie)
         for jg,ig,je,ie in zip(J025, I025, JORAS, IORAS)}

for k, (jj, ii) in enumerate(zip(Jsub, Isub)):
  if k > 0 and k % 50000 == 0:
    prc = k/len(Jsub)*100.
    print(f"  {prc:.2f}% processed")
  key = (int(jj), int(ii))
  if key not in gmapi:
    raise ValueError(f"Missing gmapi entry for UFS mesh025 index {key}")
  JORsub[k], IORsub[k] = gmapi[key]

check_gmapi = True
if check_gmapi:
  mml.check_oras_mesh025_coord(LATOR, LONOR, IORsub, JORsub, hlon, hlat, Isub, Jsub)

# gmapi indices ERA5 --> mesh025
dflgmapi = os.path.join(pthgmapi, 'gmapi_ERA5_to_UFSmesh025_global.nc')
with xr.open_dataset(dflgmapi) as ds:
  LONE = ds["era_longit"].values
  LATE = ds["era_latit"].values
  IERA = ds["era_indx"].values.astype(int)
  JERA = ds["era_jndx"].values.astype(int)
  I025 = ds["mesh025_indx"].values.astype(int)
  J025 = ds["mesh025_jndx"].values.astype(int)

JEsub = np.empty(len(Jsub), dtype=int)
IEsub = np.empty(len(Isub), dtype=int)

# Build a lookup table for fast finding corresponding ERA indices
print("Building ERA5 remapping lookup table ...")
gmapi = {(jg,ig):(je,ie)
         for jg,ig,je,ie in zip(J025, I025, JERA, IERA)}

for k, (jj, ii) in enumerate(zip(Jsub, Isub)):
  if k > 0 and k % 50000 == 0:
    prc = k/len(Jsub)*100.
    print(f"  {prc:.2f}% processed")
  key = (int(jj), int(ii))
  if key not in gmapi:
    raise ValueError(f"Missing gmapi entry for UFS mesh025 index {key}")
  JEsub[k], IEsub[k] = gmapi[key]

check_gmapi = True
if check_gmapi:
  mml.check_era_mesh025_coord(LATE, LONE, IEsub, JEsub, hlon, hlat, Isub, Jsub)




# Combine all ML information:
MLPARAM = {
  "topo"      : HH,
  "LON"       : hlon,
  "LAT"       : hlat,
  "LMsk"      : LMsk,
  "IG"        : Isub,       # mesh025 analysis grid points
  "JG"        : Jsub,       # mesh025 analysis grid points
  "regn"      : regn,
  "dxy"       : dxy,
  "YS"        : YS,
  "YE"        : YE,
  "Tfrz"      : Tfrz,
  "sst_max"   : sst_max,
  "sqrt_frzdays" : sqrt_frz,
  "intgr_frz"   : intgr_frz,
  "dxy"       : dxy,
  "Hbtm_min"  : Hbtm_min,
  "JORsub"    : JORsub,
  "IORsub"    : IORsub,
  "JERsub"   : JEsub,
  "IERsub"   : IEsub,
  }

DNMB_fcst = np.arange(dnmbS, dnmbE+1)
nfcst = len(DNMB_fcst)

print(f"Start forecasts, N forecasts: {nfcst}")
for dnmb0 in DNMB_fcst:
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  rdate = YR0*10000 + MM0*100 + DD0
  print(f"Day {YR0}/{MM0}/{DD0}")

  # Derive Predictors (Input fields) 
  # hlon, hlat - grid where prediction is done (mesh025)
  # IG, JG - grid point indices where prediction is done
  # DIRS - dictionary with input / output directories
  # PREAD_NAMES - list of predictors in the right order
  # sqrt_frzdays - True: use sqrt of integrated freezing degree days
  # sst_max - max sst threshold for possible sea ice
  # standz - True: standardize predictors (False for RF)
  # regn - region of prediction
  # intgr_frz - for freeze degree days, integration period, days
  # dxy - ice length scale, used for calc. ice predictors and gird point subset
  PRED, Iconc, SST = mml.construct_predictors_day(dnmb0, fyaml, MLPARAM, PRED_NAMES, regn)

  # Construct design / predictor matrix:
  AA = np.column_stack(PRED)

  # Prediction:
  Yfcst = mlem.predict(AA)

  print(f"ML prediction ice thickness min/max: {np.min(Yfcst):.2f} / {np.max(Yfcst):.2f}")

  #Ierr = np.where(Yfcst < 0)[0]
  Yfcst[Yfcst < 0] = 0
  Yfcst[Iconc < iconc_min] = 0
  Yfcst[SST > sst_max] = 0

  if save_npz:
    pthdump = config_ml["PRED"]["pthdump"]
    os.makedirs(pthdump, exist_ok=True)
    flfcst = f"{ml_model}_ithkn_ORAS5_ERA5_mesh025_{rdate}_{regn}.npz"
    dflfcst = os.path.join(pthdump, flfcst)
    print(f"Saving ML fcst --> {dflfcst}")
    np.savez(dflfcst, Yfcst=Yfcst, JG=Jsub, IG=Isub, idim=idm, jdim=jdm)

  elif save_nc:
    pthdump = config_ml["PRED"]["pthdump"]
    os.makedirs(pthdump, exist_ok=True)
    flfcst = f"{ml_model}_ithkn_ORAS5_ERA5_mesh025_{rdate}_{regn}.nc"
    dflfcst = os.path.join(pthdump, flfcst)
    print(f"Saving ML fcst --> {dflfcst}")

    A2d = np.zeros_like(HH) * np.nan
    A2d[LMsk] = 0
    A2d[Jsub, Isub] = Yfcst

    A3d = np.expand_dims(A2d, axis=0)
    time_dnmb = np.asarray([dnmb0])
    nctitle = f"ML {ml_model} predicted ithkn, input fields derived from ORAS5 & ERA5"

    write_MLithkn_nc(dflfcst, time_dnmb, A3d, regn, nctitle)

  if save_iconc:
    pthdump = config_ml["PRED"]["pthprd"]
    os.makedirs(pthdump, exist_ok=True)
    fliconc = f"input_iconc_ORAS5_mesh025_{rdate}_{regn}.nc"
    dfliconc = os.path.join(pthdump, fliconc)
    print(f"Saving iconc --> {dfliconc}")
       
    A2d = np.zeros_like(HH) * np.nan
    A2d[LMsk] = 0
    A2d[Jsub, Isub] = Iconc
  
    A3d = np.expand_dims(A2d, axis=0)
    time_dnmb = np.asarray([dnmb0])
    nctitle = f"Input / predictor ice concentration fields from ORAS5 on mesh025"
  
    write_iconc_nc(dfliconc, time_dnmb, A3d, regn, nctitle)
   
  if save_input:
    # Save input fields for inspection
    pthdump = config_ml["PRED"]["pthprd"]
    fldump = f"{ml_model}_predictors_ORAS5_ERA5_mesh025_{rdate}_{regn}.npz"
    dfldump = os.path.join(pthdump, fldump)
    print(f"Saving predictor fields --> {dfldump}")
    np.savez(dfldump, PredFlds=AA, JG=Jsub, IG=Isub, idim=idm, jdim=jdm)
  
print("All Done")

