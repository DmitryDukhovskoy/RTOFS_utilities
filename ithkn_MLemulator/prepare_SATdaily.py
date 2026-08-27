"""
  Prepare GDAS atm surface temperature for 
  creating predictor fields for ML ithkn emulators
  - Calculate daily mean SAT
  - subset North and South polar regions
  - Interpolate onto mesh025 grid

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

init_date = 20250701
init_hr = 0

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", 
                    choices=['north','south','global'],
                    type=str, required=True)
parser.add_argument("--rdate",  help="GDAS date YYYYMMDD", type=int, required=True)
parser.add_argument("--hr", help="GDAS hour(s) for daily avrg., list, default=0,6,12,18",
                    type=int,
                    nargs="+")

args = parser.parse_args()
regn      = args.regn if args.regn else None
gdas_date = args.rdate
HRS       = args.hr if args.hr is not None else [0,6,12,18]


pthgdas = f"/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/GFSv17/gdas.{gdas_date}"

def read_gdas(dflgdas, k2c=True, flip_north=True):
  assert os.path.isfile(dflgdas), f"Not found: {dflgdas}"
  with xr.open_dataset(dflgdas) as ds:
    A2d = ds["tmp2m"].isel(time=0).values

  if k2c:
    A2d += -273.15   # K --> Celsius

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


T2m_mean = None
for ik, hrz in enumerate(HRS):
  flname = f"gdas.t{hrz:02d}z.sfc.f000.nc"
  print(f"Reading {flname}")
  dflgdas = os.path.join(pthgdas, flname)
  T2m = read_gdas(dflgdas)

  if ik == 0:
    LON, LAT = get_gdas_coord(dflgdas)

  if T2m_mean is None:
    T2m_mean = T2m.copy()
  else:
    T2m_mean += T2m

T2m_mean /= ik

# Interpolate to mesh025




