"""
  On Gaea - change time variable 
  in  DATM_GFS.cpl.r.YYYY-MM-DD-00000.nc
  and save it to the file with requested restart date

  Original cpl.r. file can be created using cold start 

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
rmin = 0.  # restart minutes, 0,..., 59

parser = argparse.ArgumentParser()
parser.add_argument("--rdate", help="New restart date YYYYMMDD", required=True, type=int)
parser.add_argument("--sec", help="New Restart seconds =0,..,86400, default=0", type=int)
parser.add_argument("--rdold", help="Old restart date YYYYMMDD, default=rdate", type=int)
parser.add_argument("--secold", help="Old restart seconds in filename YYYY-MM-DD-sec.nc", type=int)
args = parser.parse_args()

rdate_new = args.rdate if args.rdate else None
sec_new   = args.sec if args.sec is not None else 0
rdate_old = args.rdold if args.rdold else rdate_new
sec_old   = args.secold if args.secold else 0

#pthrst  = Path('/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_fields')
pthrold = '/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/ufs_datm_mx025_cold/RESTART'
pthrnew = '/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_fields/new'

assert os.path.isdir(pthrold), f"Does not exist {pthrold}"
assert os.path.isdir(pthrnew), f"Does not exist {pthrnew}"

sfx_name = 'DATM_GFS'

# Old restart date and time (seconds):
rhr_old = sec_old //86400
dnmb_old = mtime.rdate2datenum(rdate_old*100+rhr_old)
yr_old, mm_old, dd_old, hr_old, _ = mtime.datevec(dnmb_old, round_hrs=True)
date_old = f"{yr_old}-{mm_old:02d}-{dd_old:02d}-{sec_old:05d}"
dflold = os.path.join(pthrold,f"{sfx_name}.cpl.r.{date_old}.nc")


# New restart date and time:  
rhr_new = sec_new // 86400
dnmb_rest = mtime.rdate2datenum(rdate_new*100+rhr_new)
yr_rest, mm_rest, dd_rest, hr_rest, _ = mtime.datevec(dnmb_rest, round_hrs=True)
dateout = f"{yr_rest}-{mm_rest:02d}-{dd_rest:02d}-{sec_new:05d}"
dflnew = os.path.join(pthrnew,f"{sfx_name}.cpl.r.{dateout}.nc")

assert os.path.isfile(dflold), f"Does not exist: {dflold}"

print(f"Processing {dflold}")

# Time in cpl.r is # of days since the restart day 0:00hr
ds = xr.open_dataset(dflold, decode_times=False)
ds = ds.copy(deep=True)
# Check if 'Time' is a coordinate and remove it to change for a new value
# Note Time is an index coordinate (Time dimension)
time_new = 0.0
if 'time' in ds.variables:
    ds['time'].values[:] = time_new

# Optional: set time_bnds if it exists
if 'time_bnds' in ds.variables:
    ds['time_bnds'].values[:] = np.array([[time_new, time_new]], dtype=ds['time_bnds'].dtype)

ds.attrs['info']  = f"Created from cold-run restart: {dflold}" 

print(f"Saving cice restart ---> {dflnew}")

# NetCDF - 3 with 64-bit offest does not support int64, need to manually downscale those:
# Fill value apply only for float value:
encoding = {}
for var in ds.data_vars:
  dtype = ds[var].dtype
  if np.issubdtype(dtype, np.floating):
    encoding[var] = {'_FillValue': 1.e+30}
  elif np.issubdtype(dtype, np.integer):
    encoding[var] = {'_FillValue': None}
  else:
    # For unsupported types (e.g., strings), omit FillValue
    encoding[var] = {}

#ds.to_netcdf(dflice_rest_new)
# To get rid off _FillValue attribute automatically added:
#ds.to_netcdf(dflice_rest_new, encoding={var: {'_FillValue': None} for var in ds.data_vars})
# To keep same format as original netcdf:
ds.to_netcdf(dflnew, encoding=encoding, format='NETCDF3_64BIT')
ds.close()


