"""
  Derive GLORYS daily climatology fields

location for monthly GLORYS means for NEP region:
/archive/e1n/datasets/GLORYS/monthly_means/

location for padded monthly GLORYS means, concatenated by year used to generate NEP clim nudging files:
/archive/e1n/datasets/GLORYS/monthly_climatologies/

location for monthly GLORYS means, regridded to NEP for nudging as individual months:
/archive/e1n/mom6/NEP/sponge/monthly_sponge_files/
/archive/e1n/mom6/NEP/sponge/glorys/nep_10k

location for padded monthly GLORYS means, regridded to NEP for nudging and concatenated by year:
/archive/e1n/mom6/NEP/sponge/clims/

The last directory contains the files used for nudging the solution to GLORYS. 

Daily GLORYS reanalysis for NEP domain prepared by Liz:
/archive/e1n/datasets/GLORYS/YYYY/nep_10

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import matplotlib.colors as colors
import argparse
import time

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
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
importlib.reload(mutob)


parser = argparse.ArgumentParser()
parser.add_argument("--yrs", help=f"start year to compute clim: 1993, ..., 2024", type=int)
parser.add_argument("--yre", help=f"end year to compute clim ", type=int)
args = parser.parse_args()

# Default:
YRS = 1993
YRE = 2024
if args.yrs:
  YRS= args.yrs
if args.yre:
  YRE = args.yre

VARS = ['thetao','so','uo','vo']

def save_netcdf(ds_out,dflout):
  print(f"Saving to {dflout}")
  ds_out.to_netcdf(
    dflout,
    format='NETCDF4',
    engine='netcdf4',
  )


# Compute daily climatologies
for jday in range(1,366):
  dnmb0 = mtime.jday2dnmb(1993,jday)
  yr0,mm0,dd0 = mtime.datevec(dnmb0)[:3]
  iyr = 0
  for YR in range(YRS,YRE+1):
    timeS = time.time()
    pthglorys = f'/archive/e1n/datasets/GLORYS/{YR}/nep_10/filled'
    dnmbS = int(mtime.datenum([YR,1,1]))
    dnmbE = int(mtime.datenum([YR,12,31]))
    iyr += 1

    flin = f'GLORYS_REANALYSIS_NEP_{YR}-{mm0:02d}-{dd0:02d}.nc'
    dflin = os.path.join(pthglorys,flin)
    print(f'Processing {dflin}')

    # Clim fields:
    pthclm = '/archive/Dmitry.Dukhovskoy/datasets_NEP/GLORYS_clim'
    flclm = f'GLORYS_CLIM_{YRS}-{YRE}_{mm0:02d}{dd0:02d}.nc'
    dfclm = os.path.join(pthclm,flclm)

    ds_glorys = xarray.open_dataset(dflin)
    ds_glorys = ds_glorys.drop_vars('time')
    if iyr == 1:
      ds_clm = ds_glorys.copy()

    for varnm in VARS:
      print(f'   Processing {varnm} ...')
      ds_clm[varnm] = ds_clm[varnm] + ds_glorys[varnm]
     
    timeE = time.time()
    print(f'Elapsed time: {(timeE-timeS)*1./60.:.3f} min ')
 
  # Save daily fields:
  for varnm in VARS:
    ds_clm[varnm] /= float(iyr)

  save_netcdf(ds_clm,dfclm)



