"""
  Change time value of the saved restart files 

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
parser.add_argument("--datein", help="Old restart date YYYYMMDD", type=int, required=True)
parser.add_argument("--dateout", help="New restart date YYYYMMDD", type=int, required=True)
parser.add_argument("--hrin", help="Old Restart hour, deafult=0", default=0, type=int)
parser.add_argument("--hrout", help="New Restart hour, default=0", default=0, type=int)
parser.add_argument("--nmem", help="Ensemble member number", type=int, required=True)
parser.add_argument("--prfx", help="prefix in file name prfx_YYYYMMDD.XXXXXX.analysis.cice_model.res.nc",
                    type=str, default=None)
args = parser.parse_args()

dstamp_old = args.datein
dstamp_new = args.dateout
hr_old = args.hrin
hr_new = args.hrout
nmem   = args.nmem
prfx   = args.prfx

# Original restarts:
pthinp = "/gpfs/f6/sfs-emc/world-shared/Dmitry.Dukhovskoy/restart_fields/"+\
         f"cice6_restart_fixed/enkfgdas.{dstamp_old}/{hr_old:02d}/mem{nmem:03d}"

# Get date of the 3hr back from the init_date:
dnmb0 = mtime.rdate2datenum(dstamp_old*100 + hr_old)  # restart day nmb
dnmbP = dnmb0 - 1./8. # assuming previous date/time is 3 hr back
yrP, mmP, ddP, hrP = mtime.datevec(dnmbP, round_hrs=True)[:4]

# Modified dates:
dnmb_new = mtime.rdate2datenum(dstamp_new*100 + hr_new)
yrN, mmN, ddN = mtime.datevec(dnmb_new, round_hrs=True)[:3]
if prfx is None:
  flrst_old = f"{yrP}{mmP:02d}{ddP:02d}.{hrP:02d}0000.analysis.cice_model.res.nc"
  flrst_new = f"{yrN}{mmN:02d}{ddN:02d}.{hr_new:02d}.analysis.cice_model.res.nc"
else:
  flrst_old = f"{prfx}_{yrP}{mmP:02d}{ddP:02d}.{hrP:02d}0000.analysis.cice_model.res.nc"
  flrst_new = f"{prfx}_{yrN}{mmN:02d}{ddN:02d}.{hr_new:02d}.analysis.cice_model.res.nc"


# Modified restarts:
pthout = pthinp
os.makedirs(pthout, exist_ok=True)

dfold = os.path.join(pthinp, flrst_old)
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

