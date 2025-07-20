"""
  For seasonal forecasts with COBALT (NEP_BGC)
  need to create OB fields for TA (total alk) and DIC

  For the hindcasts, those were created from T/S fields
  using ESPER subroutines that provide time-varying (annual fields)
  of TA and DIC

  For seas. f/casts, create esper clim. 
"""
import numpy as np
import xarray

import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
#import pickle
import matplotlib.pyplot as plt
from yaml import safe_load

import mod_utils_ob as mutob
importlib.reload(mutob)

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
sys.path.append('./seasonal-workflow')
from boundary import Segment
import mod_time as mtime
import mod_mom6 as mmom6


def modulo(ds):
  ds['time'] = np.arange(0, 365, dtype='float')
  ds['time'].attrs['units'] = 'days since 0001-01-01'
  ds['time'].attrs['calendar'] = 'noleap'
  ds['time'].attrs['modulo'] = ' '
  ds['time'].attrs['cartesian_axis'] = 'T'
  return ds

pthbgc = '/work/Dmitry.Dukhovskoy/NEP_input/BGC_esper_seasfcast'
flbgc  = 'bgc_esper_annual_1993_2024.nc'
flbgc_out = 'bgc_esper_annual_clim.nc'




