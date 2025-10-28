"""
  Modify snowfall rate in GFS atm. forcing fields
  unsed in datm UFS MOM6-CICE6 runs

  for CICE6 sensitivity runs

  Atm. fields prepared from GFS GDAS atm. output

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys 
import importlib
import datetime
import xarray 
import argparse
from pathlib import Path
from yaml import safe_load

# Append custom module paths
PPTHN = None
if 'PPTHN' not in locals() or PPTHN is None:
  cwd = os.getcwd()    
  parts = cwd.split(os.sep)
  if 'python' in parts:
    idx = parts.index('python')
    PPTHN = os.sep + os.path.join(*parts[:idx + 1])
  else:
    raise RuntimeError("Directory 'python' not found in current working directory path.")

sys.path.extend([
    os.path.join(PPTHN, 'MyPython', 'hycom_utils'),
    os.path.join(PPTHN, 'MyPython', 'draw_map'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])  
  
    
from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

init_date = 20250103
init_hr = 0
end_date = 20250118
end_hr = 18
  
parser = argparse.ArgumentParser()
parser.add_argument("--initd", help=f"init date, default={init_date}", type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--endd", help=f"end date of atm. forcing, defeault={end_date}", type=int)
parser.add_argument("--ehr", help=f"end hour on the endd, default={end_hr}", type=int)
parser.add_argument("--incr", help=f"snowfall increment factor: 2, 3, ...", type=int, required=True)
args = parser.parse_args()
  
init_date = args.initd if args.initd else init_date
init_hr   = args.ihr if args.ihr else init_hr
end_date  = args.endd if args.endd else end_date
end_hr    = args.ehr if args.ehr else end_hr
incr      = args.incr if args.incr else None

assert incr > 1, f"incr={incr} is invalid, incr should be > 1"
 
syst_info = os.uname()
machine = syst_info.nodename
  
if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
else:
  print("Unknown machine:", machine)
     
fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

pthrun = pths_ufs[node_nm]["MOM6"]["pthrun"]
pthatm = os.path.join(pthrun,"gfs_atm_forcing")

dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
yrI,mmI,ddI,hrI = mtime.datevec(dnmbI, round_hrs=True)[:4]
dateI_stamp = init_date*100+init_hr

dnmbE = mtime.rdate2datenum(end_date*100+end_hr)  # end day nmb
yrE,mmE,ddE,hrE = mtime.datevec(dnmbE, round_hrs=True)[:4]
dateE_stamp = end_date*100+end_hr

# Input file:
flin = f"gfs.{dateI_stamp}_{dateE_stamp}.merged.nc"
dflin = os.path.join(pthatm,flin)

# Output file:
flout = f"gfs.{dateI_stamp}_{dateE_stamp}.sf{incr}x.nc"
dflout = os.path.join(pthatm,flout)

varnm = 'fprecp'  # "surface snow precipitation rate" kg/(m2*s)

ds_atm = xarray.open_dataset(dflin)
ds_atm = ds_atm.copy(deep=True)
ds_atm[varnm][:] = ds_atm[varnm]*float(incr)

ds_atm.attrs['info']  = f"Modified {varnm} x {incr}"

print(f"Saving cice restart ---> {dflout}")
ds_atm.to_netcdf(dflout, format='NETCDF3_64BIT')
ds_atm.close()







