"""
  Change time value of the saved restart files 
  from spinup to initialize hindcast runs

  Also need to
  modify coupler.res - Current model time

"""
import numpy as np
import os
from pathlib import Path
import importlib
import sys
import xarray
import argparse
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
sys.path.append('./seasonal-workflow')
import mod_time as mtime
import mod_mom6 as mmom6

# Original restarts:
pthinp = '/archive/Dmitry.Dukhovskoy/fre/NEP/hindcast_bgc/NEPbgc_nudged_spinup/restart/1996'

# Modified restarts:
pthout = '/work/Dmitry.Dukhovskoy/NEP_input/restart_bgc'

# Modified dates:
YRR = 1993
MMR = 1
DDR = 1
dstamp = f'{YRR}{MMR:02d}{DDR:02d}'

# Count day number since 0001-01-01
dnmbR = mtime.datenum([YRR,MMR,DDR]) - 1 


def change_time(dfin,dfout,dnmbR):
  ds = xarray.open_dataset(dfin, decode_times=False)
  told = ds['Time'].data[0]
  # Replace the 'Time' coordinate with new value
  ds = ds.assign_coords(Time=[dnmbR])
  #ds['Time'].values = [dnmbR]  <-- this won't work since Time is dimension coordinate
  # Add attributes to the Time coordinate
  ds['Time'].attrs['long_name'] = 'Time'
  ds['Time'].attrs['units'] = 'days'
  ds['Time'].attrs['cartesian_axis'] = 'T'

  # Define encoding for Time to prevent _FillValue
  encoding = {
      'Time': {
          '_FillValue': None  # Disable _FillValue for Time
      }
  }

  print(f'Old restart={told:.1f}, New restart={dnmbR:.1f}, saving --> {dfout}')
  ds.to_netcdf(dfout, encoding=encoding)


flrst = 'ice_model.res.nc'
flnew = f'ice_model_{dstamp}.res.nc'
dfin  = os.path.join(pthinp,flrst)
dfout = os.path.join(pthout,flnew)
change_time(dfin,dfout,dnmbR)

for imm in range(7):
  if imm == 0:
    flrst = 'MOM.res.nc'
    flnew = f'MOM_{dstamp}.res.nc'
  else:
    flrst = f'MOM.res_{imm}.nc'
    flnew = f'MOM_{dstamp}.res_{imm}.nc'

  dfin  = os.path.join(pthinp,flrst)
  dfout = os.path.join(pthout,flnew)
  change_time(dfin,dfout,dnmbR)


