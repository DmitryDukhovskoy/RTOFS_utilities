"""
  Redistribute  relax ice fields by thcikness categories
  redistribute evenly
  then excess move to thickest cat
  if not enough 
  move from thickest to thinnest

"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
import pickle
from yaml import safe_load

PPTHN = '/home/Dmitry.Dukhovskoy/python'
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc
from mod_utils_fig import bottom_text
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

YR0 = 1993
MM0 = 5
ifld = 'ithkn'  # ithkn, iarea
file_type = 'monthly'  # monthly, daily, ... or clim
                       # for climatologies, do not need padded time - data will be recycled
                       # for monthly, daily, etc. need -dt and +dt at the beginn/end 
YR1 = YR0
YR2 = YR1+1

""" 
ICAT in ARC:
 check: CatIce=10
distribute_ice2cats thkn cat 1 hLim= 0.000
distribute_ice2cats thkn cat 2 hLim= 0.100
distribute_ice2cats thkn cat 3 hLim= 0.300
distribute_ice2cats thkn cat 4 hLim= 0.700
distribute_ice2cats thkn cat 5 hLim= 1.100
distribute_ice2cats thkn cat 6 hLim= 1.500
distribute_ice2cats thkn cat 7 hLim= 2.000
distribute_ice2cats thkn cat 8 hLim= 2.500
distribute_ice2cats thkn cat 9 hLim= 3.000
distribute_ice2cats thkn cat 10 hLim= 3.500
distribute_ice2cats thkn cat 11 hLim= 4.000
"""
# ICAT in NEP:
ICAT = np.array([1.0e-10, 0.1, 0.3, 0.7, 1.1])
# ICAT in ARC:
#ICAT = np.array([1.0e-10, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5, 3.0, 3.5])

hice = 2.5
aice = 0.85
ncat = ICAT.shape[0]
vin = np.zeros((ncat))
hin = np.zeros((ncat))
ain = np.zeros((ncat))
dhm = np.zeros((ncat))

vin[:] = hice / ncat
ain[:] = aice / ncat
dhm[:ncat-1] = 0.5*np.diff(ICAT) + ICAT[:ncat-1]
dhm[-1] = ICAT[-1] + 50.
dlt_vin = vin - dhm*ain
vin_new = vin - dlt_vin
vin_new[-1] = vin[-1] + np.sum(dlt_vin[:ncat-1])




