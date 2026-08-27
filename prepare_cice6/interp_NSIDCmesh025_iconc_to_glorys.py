"""
  Interpolate iconc field NSIDC on mesh025
  to GLORYS grid 

  For running ML estimator for generating ithkn field
  Note ML taining performed on PPAN machines

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
from datetime import datetime, timedelta
import matplotlib.dates as mdates

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
from mod_mom6 import dx_dy
import mod_icepredict as micepr
from mod_utils_fig import bottom_text
import mod_colormaps as mclrmps
import calendar

parser = argparse.ArgumentParser()
parser.add_argument("--sdate", help="Start prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--edate", help="End prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--model", help="Model names to analyze",
                    choices=['clim','ols1','ols2','rf1','rf2','rf3','gbr1','gbr2','gbr3'],
                    required=True,
                    type=str,
                    nargs="+")
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--ys", help="Year start, default=1993", default=1993, type=int)
parser.add_argument("--ye", help="Year end, default=2025", default=2025, type=int)
args = parser.parse_args()


