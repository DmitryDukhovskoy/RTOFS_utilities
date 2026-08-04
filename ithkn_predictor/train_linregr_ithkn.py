"""
  As a first step in developing ice thickness predictor
  train a linear reg. model to test different predictors
  and check if the taks is feasable, e.g. if it can be treated
  as linear problem

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
from mpl_toolkits.basemap import Basemap, cm
import argparse
from pathlib import Path
from yaml import safe_load
import statsmodels.api as statm

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
parser.add_argument("--dxy", help=f"Min dist (km) between data points (~corr.scale), to skip close i,j points", 
                    type=int, required=True)
parser.add_argument("--ys", help="Year start, default=1993", default=1993, type=int)
parser.add_argument("--ye", help="Year end, default=2025", default=2025, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'], 
                    required=True, type=str)
args = parser.parse_args()

dxy   = args.dxy    
YS    = args.ys
YE    = args.ye
regn  = args.regn

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

# Read time array ad J,I sample grid points:
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
    These are created in derive_*.py
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
A = np.column_stack((np.ones(len(Y)), *PRED_LIST))

# Standardized predictors:
# Include intercept:
Astdz = np.column_stack((np.ones(len(Y)), *PRED_STDZ))

# Check means of the predictors - should be ~0:
print("Means and stdev of the standardized predictors (except: cosD, sinD, Xcrd, Ycrd)")
for ipred, (name, x) in enumerate(zip(PRED_NAMES, PRED_STDZ)):
  print(f"x{ipred+1} {name}:   {x.mean()},   {x.std()}")


# Simple lin. regr:
# divU = 7
#Asmp = np.column_stack((np.ones(len(Y)), PRED_STDZ[7]))
#msimp = statm.OLS(Y, Asmp)
#ressimp = msimp.fit()
#print(ressimp.summary())

# No intercept:
#Astdz_nointrcp = np.column_stack(PRED_STDZ)
#model_nointrcp = statm.OLS(Y, Astdz_nointrcp)
#results_nointrcp = model_nointrcp.fit()

model = statm.OLS(Y, Astdz)
results = model.fit()

print(results.summary())

"""
coef = results.params          # coefficients
stderr = results.bse           # standard errors
pvalues = results.pvalues
tvalues = results.tvalues
confint = results.conf_int()
r2 = results.rsquared
r2adj = results.rsquared_adj

yfit = results.fittedvalues
resid = results.resid
"""

f_save = True
if f_save:
  model_name = "OLS_model1"
  flstat = f"{model_name}_{YS}_{YE}.pkl"
  dflstat = os.path.join(pthout, flstat)
  print(f"Saving stat model --> {dflstat}")
  results.save(dflstat)

  # Predictor names, mean and stdev and model constants
  flinfo = f"{model_name}_{YS}_{YE}_info.npz"
  dflinfo = os.path.join(pthout, flinfo)
  print(f"Saving: mean. stdev, predict. names --> {dflinfo}")
  np.savez(dflinfo, 
            PRED_NAMES=PRED_NAMES,
            PRED_MEAN=PRED_MEAN,
            PRED_STDEV=PRED_STDEV,
            Tfrz=Tfrz,
            sqrt_frzdays=sqrt_frzdays,
            intgr_time=intgr_time,
            ndays_era=ndays_era,
            dxy=dxy,
            YS=YS,
            YE=YE
           )


plt.ion()

fld_plot = 'frzdays'
idx = PRED_NAMES.index(fld_plot) + 1 # offset one's for intercept
#PR = A[:,idx]
PR = Astdz[:,idx]

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
ax1.scatter(PR, Y, s=5, color=(0.5,0.6,0.9), marker='.')
ax1.set_xlabel(fld_plot)
ax1.set_ylabel('ithkn')

#ax1.contour(LMsk, [0.99], linestyles='solid', colors=[(0.5,0.8,1)])
#ax1.scatter(IG, JG, s=5, color=(0.8,0.4,0), marker='.')


