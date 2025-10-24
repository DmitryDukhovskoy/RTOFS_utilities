"""
  On Gaea - change restart date 
  in MOM restart files
  for given date
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import datetime
import xarray as xr
import argparse
from pathlib import Path

#PPTHN = '/home/Dmitry.Dukhovskoy/python'
PPTHN = None
if PPTHN is None:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_mom6 as mom6util
import mod_regmom as mrgm

# Default values:
momhr = 21  # Input MOM restart, hour
icehr = 3   # Input CICE restart, hour

parser = argparse.ArgumentParser()
parser.add_argument("--rdate", help="Created restart date YYYYMMDD", required=True, type=int)
parser.add_argument("--rhr", help="Created restart hour =0,..,24, default=0", type=int)
parser.add_argument("--momdate", help="Input MOM restart date YYYYMMDD, if None - no MOM restart will be done", type=int)
parser.add_argument("--momhr", help=f"Input MOM restart hour, default={momhr}", type=int)
parser.add_argument("--icedate", help="Input CICE restart date YYYYMMDD, if None - non CICE restart will be done", type=int)
parser.add_argument("--icehr", help=f"Input CICE restart hour, default={icehr}", type=int)
args = parser.parse_args()

rdate   = args.rdate if args.rdate else None
rhr     = args.rhr if args.rhr is not None else None
momdate = args.momdate if args.momdate else None
momhr   = args.momhr if args.momhr is not None else momhr
icedate = args.icedate if args.icedate else None
icehr   = args.icehr if args.icehr is not None else icehr

pthrst  = Path('/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_fields')
pthrnew = '/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_fields/new'
dateout = f"{rdate}.{rhr:02d}0000"

assert os.path.isdir(pthrst), f"Does not exist {pthrst}"

# New restart date and time:  
dnmb_rest = mtime.rdate2datenum(rdate*100+rhr)
yr_rest, mm_rest, dd_rest, hr_rest, _ = mtime.datevec(dnmb_rest, round_hrs=True)

# Get input and output date numbers:
if momdate is not None:
  dnmb_mom = mtime.rdate2datenum(momdate*100 + momhr)
  yr_mom, mm_mom, dd_mom, hr_mom, _ = mtime.datevec(dnmb_mom, round_hrs=True)

  assert yr_mom*10000+mm_mom*100+dd_mom == momdate and hr_mom == momhr, \
          f"converted dates or time do not match input {momdate}:{momhr}"

  # Time stamps for file naming
  mom_datein = f"{momdate}.{momhr:02d}0000"

  mom_dayoffset = 14  # day number in MOM restart is +14 days of what it should be for file datestamp

  # For MOM6, add ofsset to the New restart date and time:
  dnmb_restoff = dnmb_rest + mom_dayoffset


  # Process MOM6 restarts:
  momrest_list = sorted(pthrst.glob(f"{mom_datein}.MOM.res*.nc"))

  for flnm in momrest_list:
    print(f"Processing {flnm.name}")
    ds = xr.open_dataset(flnm, decode_times=False)
    
    # Make a deep copy to modify
    ds = ds.copy(deep=True)

    # Check if 'Time' is a coordinate and remove it to change for a new value
    # Note Time is an index coordinate (Time dimension)
    if 'Time' in ds.coords:
      ds = ds.assign_coords(Time=("Time", [dnmb_restoff]))

    flmom = flnm.stem
    sfx = '.'.join(flmom.split('.')[-2:])
    fmomrst_new = f"{sfx}.{dateout}.nc"
    dfl_out = os.path.join(pthrnew, fmomrst_new)

    # Save to new NetCDF
    print(f"Saving restart --> {dfl_out}")
    #ds.to_netcdf(dfl_out)
    # check the original format with ncdump -k MOM.res.nc
    # keep same format as original netcdf - MOM saved in netCDF-4
    ds.to_netcdf(dfl_out, encoding={var: {'_FillValue': None} for var in ds.data_vars})

    ds.close()

if icedate is not None:
  dnmb_ice = mtime.rdate2datenum(icedate*100 + icehr)
  yr_ice, mm_ice, dd_ice, hr_ice, _ = mtime.datevec(dnmb_ice, round_hrs=True)

  assert yr_ice*10000+mm_ice*100+dd_ice == icedate and hr_ice == icehr, \
          f"converted dates or time do not match input {icedate}:{icehr}"

  ice_datein = f"{icedate}.{icehr:02d}0000"
  flice_rest = f"{ice_datein}.cice_model.res.nc"
  dflice_rest = os.path.join(pthrst, flice_rest)
  print(f"\nProcessing {dflice_rest}")

  ds = xr.open_dataset(dflice_rest)
  ds.attrs['myear']  = np.int32(yr_rest)
  ds.attrs['mmonth'] = np.int32(mm_rest)
  ds.attrs['mday']   = np.int32(dd_rest)
  ds.attrs['msec']   = np.int32(hr_rest * 3600)

  # Get rid off automatically added _FillValue = NaN:
  #for var in ds.data_vars:
  #  if '_FillValue' in ds[var].encoding:
  #    del ds[var].encoding['_FillValue']

  dflice_rest_new = os.path.join(pthrnew,f"cice_model.res.{dateout}.nc")
  print(f"Saving cice restart ---> {dflice_rest_new}")

  #ds.to_netcdf(dflice_rest_new)
  # To get rid off _FillValue attribute automatically added:
  #ds.to_netcdf(dflice_rest_new, encoding={var: {'_FillValue': None} for var in ds.data_vars})
  # To keep same format as original netcdf:
  ds.to_netcdf(dflice_rest_new, encoding={var: {'_FillValue': None} for var in ds.data_vars}, format='NETCDF3_64BIT')
  ds.close()


