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
dstamp_old='19960101'
pthinp = '/archive/Dmitry.Dukhovskoy/fre/NEP/hindcast_bgc/NEPbgc_nudged_spinup/'+\
         f'restart/restdate_{dstamp_old}'

# Modified restarts:
pthout = '/work/Dmitry.Dukhovskoy/NEP_input/restart_bgc'
os.makedirs(pthout, exist_ok=True)

# Modified dates:
YRR = 1993
MMR = 1
DDR = 1
dstamp = f'{YRR}{MMR:02d}{DDR:02d}'
yr_rest = int(dstamp[:4])
mm_rest = int(dstamp[4:6])
dd_rest = int(dstamp[6:8])

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


flrst = f'ice_model_{dstamp_old}.res.nc'
flnew = f'ice_model_{dstamp}.res.nc'
dfin  = os.path.join(pthinp,flrst)
dfout = os.path.join(pthout,flnew)
change_time(dfin,dfout,dnmbR)

for imm in range(8):
  if imm == 0:
    flrst = f'MOM_{dstamp_old}.res.nc'
    flnew = f'MOM_{dstamp}.res.nc'
  else:
    flrst = f'MOM_{dstamp_old}.res_{imm}.nc'
    flnew = f'MOM_{dstamp}.res_{imm}.nc'

  dfin  = os.path.join(pthinp,flrst)
  dfout = os.path.join(pthout,flnew)
  change_time(dfin,dfout,dnmbR)

# COBALT restarts do not have Time as a variable
import shutil

for sfx in ['ice_cobalt', 'ocean_cobalt_airsea_flux']:
  flcob_in = f'{sfx}_{dstamp_old}.res.nc'
  dfinp = os.path.join(pthinp,flcob_in)
  flcob_out = f'{sfx}_{dstamp}.res.nc'
  dfout = os.path.join(pthout,flcob_out)

  print(f'Copying {dfinp} --> {dfout}')
  shutil.copy(dfinp,dfout)

# Update coupler.res:
dfnm_in = os.path.join(pthinp,  f'coupler_{dstamp_old}.res')
dfnm_out = os.path.join(pthout, f'coupler_{dstamp}.res')
print(f'Udpating {dfnm_in} --> {dfnm_out}')

new_date = [yr_rest, mm_rest, dd_rest, 0, 0, 0]

with open(dfnm_in, 'r') as infile, open(dfnm_out, 'w') as outfile:
  for line in infile:
    stripped = line.strip()

    # Find lines that start with a year
    if stripped and stripped.split()[0].isdigit() and len(stripped.split()[0]) == 4:
      # Replace old restart with the new date
      parts = line.split()
      parts[:6] = [f"{dmm:5d}" for dmm in new_date]  
      new_line = "  ".join(parts)
      outfile.write(f"{new_line}\n")
    else:
      # Leave all other lines untouched
      outfile.write(line)

print('All done')


