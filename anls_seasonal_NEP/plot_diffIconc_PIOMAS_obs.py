"""
  Plot difference maps of monthly ice conc from PIOMAS reanalysis

  Specify months (calendar numbering!) to average statistics by seasons
  and compare against NSIDC conc fields

  If more than 1 year is given in keyargs than statistics are averaged over these years
  
   monthly mean and StDev  bottom T derived in calc_mnthlyTSbtm.py
  Save by years

  Usage: plot_diffIthkn_PIOMAS_obs.py --YAS=1993 --YAE=1993 --MAS=10 --MAE=10
  use --help for more information on keywargs

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import pickle
from copy import copy
import matplotlib.colors as colors
from yaml import safe_load
import argparse
import pickle

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
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')
import mod_swstate as msw
import conversions as gsw

from mod_utils_fig import bottom_text
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_rtofs as mrtofs
importlib.reload(mutob)
importlib.reload(manseas)

parser = argparse.ArgumentParser()
parser.add_argument("--YAS", help="Calendar (! not init.) year to start stat. averaging", type=int)
parser.add_argument("--YAE", help="Calendar year to end averaging of statistics", type=int)
parser.add_argument("--MAS", help="Calendar month to start averaging of statistics", type=int)
parser.add_argument("--MAE", help="Calendar month to end averaging of statistics", type=int)
args = parser.parse_args()

# Default Averaging time period:
YAS = 1993   # Start: f/cast init. year to use for monthly averaging
YAE = YAS   # End
MAS = 10
MAE = MAS

if args.YAS:
  YAS = args.YAS
if args.YAE:
  YAE = args.YAE
else:
  YAE=YAS
if args.MAS:
  MAS = args.MAS
if args.MAE:
  MAE = args.MAE
else:
  MAE = MAS

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

run_name   = 'seasonal_fcst_daily'
pthsis  = gridfls['MOM6_NEP'][run_name]['pthsis']

icc = 0
for YRA in range(YAS,YAE+1):
  YR1 = YRA
  YR2 = YR1+1
  flout = f'PIOMAS_ithkn_iconc_{YR1}_{YR2}_monthly.nc'
  diclim = os.path.join(pthsis, flout)
  ds_rlx = xarray.open_dataset(diclim)
  Time = ds_rlx['time'].data
  TM = mmisc.convert_nptime_to_datenum(Time)
  for MMA in range(MAS,MAE+1):
    dnmb0 = mtime.datenum([YRA,MMA,15,12])
    D = abs(TM-dnmb0)
    itime = np.argmin(D)
    dv0 = mtime.datevec(TM[itime])
    assert dv0[0]==YRA, f'Requested YR={YRA}, year in rlx file={dv0[0]}'
    assert dv0[1]==MMA, f'Requested month={MMA}, month in rlx file={dv0[1]}'

    Cice = ds_rlx['iarea'].isel(time=itime).data
    if icc == 0:
      CI = Cice
    else:
      CI = CI + Cice

    icc += 1

CI = CI.squeeze()/icc

# Get NSIDC concentration averaged over same time period/seasons
CMobs,Xobs,Yobs = manseas.avrg_cice_NSIDC(YAS, YAE, MAS, MAE)

# Get ice edge contour in the Bering Sea
# get rid of ice in unneeded part of the domain
CMobs[:,150:] = np.nan
CMobs[:200,:] = np.nan

CNTR = manseas.derive_ice_contour(CMobs, nmin=10)




