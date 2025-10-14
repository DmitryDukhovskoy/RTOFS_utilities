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
parser.add_argument("--momdate", help="Input MOM restart date YYYYMMDD", required=True, type=int)
parser.add_argument("--momhr", help=f"Input MOM restart hour, default={momhr}", type=int)
parser.add_argument("--icedate", help="Input CICE restart date YYYYMMDD", required=True, type=int)
parser.add_argument("--icehr", help=f"Input CICE restart hour, default={icehr}", type=int)
args = parser.parse_args()

rdate   = args.rdate if args.rdate else None
rhr     = args.rhr if args.rhr is not None else None
momdate = args.momdate if args.momdate else None
momhr   = args.momhr if args.momhr is not None else momhr
icedate = args.icedate if args.icedate else None
icehr   = args.icehr if args.icehr is not None else icehr

pthrst = Path('/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_fields')
assert os.path.isdir(pthrst), f"Does not exist {pthrst}"

# Get input and output date numbers:
dnmb_mom = mtime.rdate2datenum(momdate*100 + momhr)
yr_mom, mm_mom, dd_mom, hr_mom, _ = mtime.datevec(dnmb_mom, round_hrs=True)

assert yr_mom*10000+mm_mom*100+dd_mom == momdate and hr_mom == momhr, \
        f"converted dates or time do not match input {momdate}:{momhr}"

dnmb_ice = mtime.rdate2datenum(icedate*100 + icehr)
yr_ice, mm_ice, dd_ice, hr_ice, _ = mtime.datevec(dnmb_ice, round_hrs=True)

assert yr_ice*10000+mm_ice*100+dd_ice == icedate and hr_ice == icehr, \
        f"converted dates or time do not match input {icedate}:{icehr}"

# Time stamps for file naming
mom_datein = f"{momdate}.{momhr:02d}0000"
ice_datein = f"{icedate}.{icehr:02d}0000"
dateout    = f"{rdate}.{rhr:02d}0000"

mom_dayoffset = 14  # day number in MOM restart is +14 days of what it should be for file datestamp

# New restart date and time:
dnmb_rest = mtime.rdate2datenum(rdate*100+rhr)
dnmb_restoff = dnmb_rest + mom_dayoffset


# Process MOM6 restarts:
momrest_list = sorted(pthrst.glob(f"{mom_datein}.MOM.res*.nc"))

for flnm in momrest_list:
  print(f"Processing {flnm.name}")
  ds = xr.open_dataset(flnm)
  
  # Make a deep copy to modify
  ds = ds.copy(deep=True)

  # Check if 'Time' is a coordinate and remove it to change for a new value
  if 'Time' in ds.coords:
    ds = ds.reset_coords('Time')

  ds['Time'][:] = dnmb0
  # Optionally set Time back as coordinate
  ds = ds.set_coords('Time')

  new_flnm = flnm.with_name(flnm.stem + '_mod.nc')

  # Save to new NetCDF
  ds.to_netcdf(new_flnm)

  ds.close()

  print(f"Saving restart --> {new_flnm}")
     


