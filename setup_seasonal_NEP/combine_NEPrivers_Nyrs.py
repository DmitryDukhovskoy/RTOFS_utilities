"""
  Prepare river data for NEP 1yr forecasts
  daily GLOFAS river fields
  Combine several years into 1 to allow for the forecasts that are initialized 
  any time during start year

  Need to overlap last year of N-yr block with the 1st year of the next N-yr block

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

PPTHN = '/home/Dmitry.Dukhovskoy/python/'
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

Nyrs = 5   # N of year in 1 file
YR1 = 1996
YR2 = 2021  # last year for generating runoff data for simulation, may be changed to have 5 years
yriver_end = 2019 # last year of year data

# List of years:
riv_years = np.zeros((1,2))
yrS = YR1
icc = 0
while yrS < YR2:
  yrE = yrS + Nyrs-1
  dmm = np.array([[yrS, yrE]])
  icc += 1
  if icc == 1:
    riv_years = dmm 
  else:
    riv_years = np.append(riv_years, dmm, axis=0)
  yrS = yrE - 1


fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)

flriv = config['filesystem']['yearly_river_files']

# Note all runoff have +1 padded days for interpolation in time
#rivers = xarray.open_mfdataset( input_files,
#        preprocess=lambda x: x.isel(time=slice(None, 365)) )# skip padded days & leap yrs

nblocks = riv_years.shape[0]
for ibl in range(nblocks):
  yr1, yr2 = riv_years[ibl,:]
  years = np.arange(yr1, yr2+1)
  input_files = [flriv.format(year=y) for y in years]
  nfls = len(input_files)

  print(f"Creating runoff for {yr1}-{yr2}")

  ds_riv = xarray.Dataset()
  for ifl in range(nfls):
    yr0 = yr1 + ifl
    print(f"-->  {yr0} ...")
    if yr0%4 == 0:
      nyrdays = 366
    else: 
      nyrdays = 365
    if yr0 <= yriver_end:
      ds = xarray.open_dataset(input_files[ifl]) 
      ds0 = ds.copy()
#      if yr0 < yr2:
      R  = ds['runoff'].data
      dmm = ds['time'].data
      strd1 = str(dmm[0])
      strd2 = str(dmm[-1])
      yr1_input = int(strd1[:4])
      mo1_input = int(strd1[5:7])
      dd1_input = int(strd1[8:10])
      yr2_input = int(strd2[:4])
      mo2_input = int(strd2[5:7])
      dd2_input = int(strd2[8:10])
      assert yr1_input == yr0, f"Input year {yr1_input} does not agree with init YR {yr0}"

      jday1 = mtime.date2jday([yr1_input, mo1_input, dd1_input])
      jday2 = mtime.date2jday([yr2_input, mo2_input, dd2_input])
      TM = mtime.npdatetime_year(yr1_input, yr2_input, day_start=jday1, day_end=jday2, tprecis='D') 
    else:
    # After 2019 - no data, use last year
    # add + 1 day - padded data
    # for leap year - add +1 day at the end
    # Assumed that last runoff year (2019) is NOT leap year !
      ds = ds0.copy()
      TM0 = mtime.npdatetime_year(yr0, yrE=yr0+1, day_end=1, tprecis='D')
      R0  = ds['runoff'].data
      if len(TM0) > R0.shape[0]:
        dmm = np.expand_dims(R0[-1,:], axis=0)
        R0 = np.append(R0, dmm, axis=0)

      if not len(TM0) == len(R0):
        raise Exception(f"Fake Time and runoff arrays for year {yr0} not equal")

      R  = R0.copy()
      TM = TM0.copy()

# Skip padded days at the end if not the last year
    ds = ds.drop_vars('runoff')
    ds = ds.drop_vars('time')
    ds = ds.assign_coords(time=TM[:-1])
    ds = ds.assign(runoff=(['time','y','x'], R[:-1,:,:]))
    if ifl == 0:
      ds_riv = ds.copy()
    else:
      ds_riv = xarray.concat([ds_riv, ds], dim="time")

# Drop unneeded dims in lon, lat
  ds_riv['area'] = ds_riv['area'].isel(time=1).drop_vars('time')
  ds_riv['lat'] = ds_riv['lat'].isel(time=1).drop_vars('time')
  ds_riv['lon'] = ds_riv['lon'].isel(time=1).drop_vars('time')


  # Add attributes for time var
#  ds_riv['time'].attrs['units'] = f'days since 1950-01-01'
#  ds_riv['time'].attrs['calendar'] = 'gregorian'
  ds_riv['time'].attrs['cartesian_axis'] = 'T'
  ds_riv['runoff'].attrs['units'] = 'kg m-2 s-1'

  idm = ds_riv.sizes['x']
  jdm = ds_riv.sizes['y']
  output_dir = config['filesystem']['river_clim_path']
  dfriv_out = os.path.join(output_dir, f'glofas_runoff_NEP_{jdm}x{idm}_daily_{yr1}-{yr2}.nc')
  print(f"Writing river daily --> {dfriv_out}")

  ds_riv.to_netcdf(
      dfriv_out, 
      format='NETCDF3_64BIT',
      engine='netcdf4',
      encoding={'time': {'dtype': 'float64', 'calendar': 'gregorian'}},
      unlimited_dims='time'
  )

f_chck = False
if f_chck:
  dsR = xarray.open_dataset(dfriv_out)
  j0 = 337
  i0 = 260

#  Ryr = vardata.isel(time=slice(None,365), y=j0, x=i0).load()
  Ryr = vardata.isel(y=j0, x=i0).load().data
  nyrs = int(Ryr.shape[0]/365)
  Ryr = Ryr.reshape((nyrs, 365)).transpose()

  Rav = rivmn.data[:,j0,i0].squeeze()

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])
  ax1.plot(Ryr, color=[0.6, 0.6, 0.6])
  ax1.plot(Rav)

  ax1.set_title(f'River climatology and original GLOFAS runoff, i0={i0}, j0={j0}')


