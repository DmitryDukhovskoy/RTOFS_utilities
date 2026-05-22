"""
  Change time value of the saved restart files 
  In CICE IC files, prepared from GDAS SOCA, file name and actual time 
  in the netcdf file may differ

  Change restart time in the file to the desired restart time,
  change file name to YYYYMMDD.HR.cice_model.res.nc

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

parser = argparse.ArgumentParser()
parser.add_argument("--flinp", help="CICE restart original file to be modified", 
                    required=True, type=str)
parser.add_argument("--dateout", help="New restart date YYYYMMDD", type=int, required=True)
#parser.add_argument("--hrin", help="Old Restart hour, deafult=0", default=0, type=int)
parser.add_argument("--hrout", help="New Restart hour, default=0", default=0, type=int)
args = parser.parse_args()

flrst_old = args.flinp
#dstamp_old = args.datein
dstamp_new = args.dateout
#hr_old = args.hrin
hr_new = args.hrout

# Original restarts:
pthinp = "/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_sfs_C192mx025/ice"
dfold = os.path.join(pthinp, flrst_old)

# Get date in the original restat file
with xarray.open_dataset(dfold) as ds:
  year  = ds.attrs["myear"]
  month = ds.attrs["mmonth"]
  day   = ds.attrs["mday"]
  msec  = ds.attrs["msec"]

dstamp_old = year*10000 + month*100 + day
hr_old = msec // 3600
dnmb0 = mtime.rdate2datenum(dstamp_old*100 + hr_old)  # restart day nmb
yrP, mmP, ddP, hrP = mtime.datevec(dnmb0, round_hrs=True)[:4]

# Modified dates:
dnmb_new = mtime.rdate2datenum(dstamp_new*100 + hr_new)
yrN, mmN, ddN = mtime.datevec(dnmb_new, round_hrs=True)[:3]

#flrst_old = f"{yrP}{mmP:02d}{ddP:02d}.{nsec:06d}.cice_model.res.nc"
flrst_new = f"{yrN}{mmN:02d}{ddN:02d}.{hr_new:02d}.cice_model.res.nc"

# For information only:
print(f"Original CICE restart: {flrst_old} --> new: {flrst_new}")
print(f"Original time in the file: {dstamp_old} msec={msec}")
print(f"New time in the file: {dstamp_new} msec={hr_new*3600}")

# Modified restarts:
pthout = pthinp
os.makedirs(pthout, exist_ok=True)

dfnew  = os.path.join(pthout, flrst_new)

ds = xarray.open_dataset(dfold)
ds.attrs['myear'] = yrN
ds.attrs['mmonth'] = mmN
ds.attrs['mday'] = ddN
ds.attrs['msec'] = hr_new*3600

print(f"Saving CICE restart --> {dfnew}")
ds.to_netcdf(
             dfnew, 
             encoding={var: {'_FillValue': None} for var in ds.data_vars}, 
             format='NETCDF3_64BIT'
)

ds.close()

