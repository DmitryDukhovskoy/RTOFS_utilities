"""
  Prepare GDAS atm surface temperature for 
  creating predictor fields for ML ithkn emulators
  - Calculate daily mean SAT (T2m)
  - subset North and South polar regions

  First step in deriving other T2m - based predictors
  integr. freeze degree days
  or SAT 

  ML models developed on PPAN

  GFSv17 status with HPSS / WCOSS directories:
  https://docs.google.com/spreadsheets/d/1N3isKTVmE4ITdiULDLP5lK1RoZOzkNlHFN-NFHwrH6o/edit?gid=492588212#gid=492588212

  Fetch GDAS fields from HPSS:
  atm GDAS: scripts/DATA_HPSS/get_GDASatm.sh
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray as xr
from yaml import safe_load
from mpl_toolkits.basemap import Basemap, cm
import argparse

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


from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_sis2_relax as msisrlx
import mod_regmom as mrmom 


parser = argparse.ArgumentParser()
parser.add_argument("--rdate",  
             help="GDAS date YYYYMMDD to calc daily average on GDAS grid", 
             type=int,
             required=True)
parser.add_argument("--hr", help="use GDAS hour(s) for daily avrg., list, default=[0,6,12,18]",
                    type=int,
                    nargs="+")
parser.add_argument("--flipN", 
           help="=1: Flip South-North GDAS grid to have North at the top, =0 - keep original GDAS grid",
           default=1,
           choices=[0,1],
           type=int)

args       = parser.parse_args()
gdas_date  = args.rdate
HRS        = args.hr if args.hr is not None else [0,6,12,18]
flip_north = args.flipN == 1

    
fyaml = 'paths_ML.yaml'
with open(fyaml) as ff:
  pths_ml = safe_load(ff)

def read_gdas(dflgdas, k2c=True, flip_north=True):
  assert os.path.isfile(dflgdas), f"Not found: {dflgdas}"
  with xr.open_dataset(dflgdas) as ds:
    A2d = ds["tmp2m"].isel(time=0).values

  if k2c:
    A2d = A2d - 273.15   # K --> Celsius

  if flip_north:
    # Flip array to have N. at the top:
    A2d = np.flipud(A2d)

  return A2d

def get_gdas_coord(dflgdas, flip_north=True):
  assert os.path.isfile(dflgdas), f"Not found: {dflgdas}"
  with xr.open_dataset(dflgdas) as ds:
    LON = ds["lon"].values
    LAT = ds["lat"].values

  if flip_north:
    # Flip array to have N. at the top:
    LON = np.flipud(LON)
    LAT = np.flipud(LAT)

  return LON, LAT

def T2m_daily_mean(pthgdas, HRS, flip_north):
  """
    Derive daily mean T2m from GDAS
    using HRS hours
  """
  T2m_mean = None
  ik = 0
  for hrz in HRS:
    flname = f"gdas.t{hrz:02d}z.sfc.f000.nc"
    print(f"Reading {flname}")
    dflgdas = os.path.join(pthgdas, flname)
    T2m = read_gdas(dflgdas, flip_north=flip_north)

    ik += 1
    if T2m_mean is None:
      T2m_mean = T2m.copy()
    else:
      T2m_mean += T2m

  T2m_mean /= ik

  return T2m_mean

dnmb0 = mtime.rdate2datenum(gdas_date)
YR, MM, DD = mtime.datevec(dnmb0)[:3]
  
pthdata = pths_ml["GDAS"]["pthdata"]  # root dir for processed data
pthgdas = pths_ml["GDAS"]["pthgdas"].format(YR=YR, MM=MM, DD=DD)  # unprocessed data from HPSS
pthdaily = pths_ml["GDAS"]["pthdaily"].format(YR=YR)  # daily T2m, GDAS grid

# Derive daily mean:
T2m_day = T2m_daily_mean(pthgdas, HRS, flip_north)

# Save T2m daily field, not interpolated
#tmp_dir = os.path.join(pthdata, 'GDAS_T2m','glob','daily_tmp')
os.makedirs(pthdaily, exist_ok=True)
if flip_north:
  tmp_file = f"GDAS_flipN_T2m_{gdas_date}.npy"
else:
  tmp_file = f"GDAS_notflipN_T2m_{gdas_date}.npy"

dfl = os.path.join(pthdaily, tmp_file)

print(f"Saving daily T2m {gdas_date} --> {dfl}")
np.save(dfl, T2m_day)


