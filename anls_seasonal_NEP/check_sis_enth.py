import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
from netCDF4 import Dataset as ncFile
import importlib
import xarray
import yaml
from yaml import safe_load

PPTHN = []
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_anls_seas as manseas


iconc_new = 0.3016
iconc_old = 0.0
ithk_new  = 2214.57    # kg m-2
ithk_old  = 2327.25

# Find change in ice = dlt_iconc*ithk_old + dlt_ithk*iconc_new
# same for snow 
dlt_iconc = iconc_new - iconc_old
dlt_ithk  = ithk_new - ithk_old
dlt_ice   = dlt_iconc * ithk_old + dlt_ithk * iconc_new

# Better way: compute mass*fract_area 
dlt_ice   = (ithk_new*iconc_new - ithk_old*iconc_old)
# Using bluk ice S:
ice_salin = 3.35
dlt_salt = dlt_ice * ice_salin  # kg m-2 * g kg-1



