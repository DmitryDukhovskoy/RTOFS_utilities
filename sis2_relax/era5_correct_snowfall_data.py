"""
  ERA5 sf (snowfall) padded forcing files for 2001 cause ice continuity 
  blow up (negative ice thickness)
  This is due to the min values being close to zero but negative: -8.673617379884035e-19
  instead of 0.0 
  In SIS2, snowfall is converted to the snowice with negative thickness

  ERA5 atm. forcing files:
  ERA5_u10_2000_padded.nc
  ERA5_v10_2000_padded.nc
  ERA5_lp_2000_padded.nc
  ERA5_msl_2000_padded.nc
  ERA5_sf_2000_padded.nc
  ERA5_sphum_2000_padded.nc
  ERA5_ssrd_2000_padded.nc
  ERA5_strd_2000_padded.nc
  ERA5_t2m_2000_padded.nc
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
from copy import copy
import matplotlib.colors as colors
from yaml import safe_load
import argparse
import pandas as pd

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
import mod_interp1D as mint1d
#importlib.reload(mutob)

parser = argparse.ArgumentParser()
parser.add_argument("--yrs", help="Start year for checking/correcting snowfall ERA5", type=int)
parser.add_argument("--yre", help="End year for checking/correcting snowfall ERA5", type=int)
#parser.add_argument("--varnm", help="field to change: u10,v10,lp,msl,sf,sphum,ssrd,strd,t2m", type=str)
args = parser.parse_args()

varnm = 'sf'
f_save = True

if args.yrs:
  YRS=args.yrs
  YRE=YRS
if args.yre:
  YRE=args.yre

import time
# Read data:
pthera = '/archive/e1n/mom6/NEP/atmos_forcing/era5_padded'
for YRR in range(YRS,YRE+1):
  print(f'Checking ERA5 sf for {YRR}')
  ptherafld = os.path.join(pthera,varnm)
  flera = f'ERA5_{varnm}_{YRR}_padded.nc'

  dflera = os.path.join(ptherafld,flera)

  print(f'Processing {YRR} {varnm}')
  print(f'Opening {dflera}')
  dset = xarray.open_dataset(dflera)
  darray = dset[varnm]
  Time = dset['time'].data
  tmP = pd.to_datetime(Time)
  nrec = len(tmP)
  years  = tmP.year.to_numpy()
  months = tmP.month.to_numpy()
  days   = tmP.day.to_numpy()
  hours  = tmP.hour.to_numpy()
  found_neg = False
  for itime in range(nrec):
    timeS = time.time()
    A2d = darray.values[itime,:,:]
    #A2d = dset[varnm].isel(time=itime).data
    sfmin = np.min(A2d)
    #Jneg,Ineg = np.where(A2d<0.0)
    if sfmin < 0.0:
      yy,mm,dd,hh = years[itime],months[itime],days[itime],hours[itime]
      print(f'Found neg. {varnm}, {yy}/{mm}/{dd}:{hh:02d}hr, replacing with 0.0')
      A2d = np.where(A2d<0., 0.0, A2d)
      darray.values[itime,:,:] = A2d
      found_neg = True
      timeE = (time.time()-timeS)
      print(f'Ellapsed time: {timeE:.3f} sec')

  if not found_neg:
    print(f' No negative values in {YRR}, ok ...')
    continue

  dset[varnm] = darray
  dset.attrs["info"]=f"Corrected near-zero negative snow fall values "
  dset.attrs["code"]="/home/Dmitry.Dukhovskoy/python/sis2_relax/era5_correct_snowfall_data.py"

  if f_save:
    pthoutp = f'/work/Dmitry.Dukhovskoy/NEP_input/ERA5_padded_changed/{YRR}'
    os.makedirs(pthoutp, exist_ok=True)

    flnew = f'ERA5_{varnm}_{YRR}_corrected_padded.nc'

    dflout = os.path.join(pthoutp, flnew)
    print(f"Saving to {dflout}")
    dset.to_netcdf(
      dflout,
      format='NETCDF3_64BIT',
      engine='netcdf4',
    )

  dset.close()

