"""
  Driver to prepare CICE6 restart for specified
  setup of initial conditions

  Can use any restart adding missing fields to it, e.g.:
  edit cice6rest_files_SFS.yaml
  rdate_out:  requested restart date in the new restart
  rdate_in, rhr_in: can be blank, the date/time will be derived from the input restart name assuming one of the formats:
                    YYYYMMDD.XX[XXXX].restart_name.*.nc
                    some_sfx.YYYYMMDD.XX[XXXX].*.nc

  pth_in, pth_out: specify input/output restart dirctories

  flrst_in: intput restart file name
  flrst_out:  leave blank, output restart name will be created based on the flrst_tmp (template) and restart time
               by adding fields being corrected, e.g. cice_restart.20240701.00.iconc_ithkn.nc

  Example:
  create restart with corrected hsnow:
  in YAML file:
  rdate_in:  ~          # input restart date
  rhr_in:    ~          # input restart hour
  rdate_out: 20240701   # output restart date

  flrst_in: "20240701.064800.cice_restart.nc"
  flrst_out: ~
  flrst_tmp: "cice_restart.YYYYMMDD.HH"

  run prepare_cice6_restart_driver.py --iconc 0 --ithkn 0 --hsnow 1 --snitd 0 --regn global

  To create restart with iconc + ithkn + hsnow
  2 options:
  (1) use ice restart from flrst_in (with already inserted iconc and ithkn) and will add
    hsnow on top of this
  flrst_in: "cice_restart.20240701.00.iconc_ithkn.nc"   <--- use restart with already added iconc and ithkn
  flrst_out: ~
  flrst_tmp: "cice_restart.YYYYMMDD.HH"

  run prepare_cice6_restart_driver.py --iconc 0 --ithkn 0 --hsnow 1 --snitd 0 --regn global

  (2) create from original restart, note that intermediate restart files will also be created and saved
  flrst_in: "20240701.064800.cice_restart.nc"   <--- use original restart
  flrst_out: ~
  flrst_tmp: "cice_restart.YYYYMMDD.HH"

  run prepare_cice6_restart_driver.py --iconc 1 --ithkn 1 --hsnow 1 --snitd 0 --regn global


  Note that for snitd = 1, use ice_in with added options for snow redistribution, snow physics

"""
import os
import numpy as np
import importlib
import sys
import xarray
from yaml import safe_load
import argparse
import subprocess

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
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])

import mod_cice6_utils as mc6util

#rest_date = 20250103
#rest_date = 20250701 - not needed, check yaml 
rest_hr   = 0
rhr_out   = rest_hr
flrst_in  = None
flrst_out = None
#fyaml = 'cice6rest_files_SFSmem.yaml'  # to run Neils restart for SFS ensmb GFS / CPC cases
#fyaml = 'cice6rest_files_SFS_RTOFS.yaml'
fyaml = 'cice6rest_files_SFS_ML.yaml'
iconc  = 1  # insert iconc from NRT NSIDC
ithkn  = 0  # insert ithkn clim
hsnow  = 0  # insert snow depth clim
snitd  = 0  # snow redistribution over ice (ITDrdg - requires extra tracers in restart) 
sstmom = 0  # adjust MOM6 SST to bring closer to Tfreeze(S) under sea ice, will modify MOM.res.nc 

parser = argparse.ArgumentParser()
parser.add_argument("--iconc", type=int, 
                    choices=[0,1], help=f"insert ice conc NRT NSIDC, default={iconc}")
parser.add_argument("--ithkn", type=int, 
                    choices=[0,1], help=f"insert ice thickn climatology, default={ithkn}")
parser.add_argument("--hsnow", type=int, 
                    choices=[0,1], help=f"insert snow depth climatology, default={hsnow}")
parser.add_argument("--snitd", type=int, 
                    choices=[0,1], help=f"use snow redistribution over ice, default={snitd}")
parser.add_argument("--regn", help=f"region where restart is being updated, default=global", 
                    choices=['south','north','global'], default="global", type=str)
parser.add_argument("--sst", help=f"adjust SST in MOM restart under sea ice, default={sstmom}",
                    choices=[0,1], default=sstmom, type=int) 
parser.add_argument("--fyaml", 
                    help=f"YAML with local directories, filenames, restart dates, default={fyaml}",
                    default=fyaml,
                    type=str)
# For GFS ensembles with different restart files under mem000, mem001, ... dirs
# Not needed otherwise
parser.add_argument("--enmb",
      help=f"ensemble number for ensamble runs, requires YAML for enmb choice", 
      default=None, type=int)
args = parser.parse_args()

regn  = args.regn
iconc = args.iconc  if args.iconc  is not None else iconc
ithkn = args.ithkn  if args.ithkn  is not None else ithkn
hsnow = args.hsnow  if args.hsnow  is not None else hsnow
snitd = args.snitd  if args.snitd  is not None else snitd
fyaml = args.fyaml  if args.fyaml  is not None else fyaml 
enmb  = args.enmb 
sstmom = args.sst

print(f"Reading YAML with restart info: {fyaml}\n")
with open(fyaml) as ff:
  config_rest = safe_load(ff)

# Output Restart dates if missing - same as input
# input dates are deduced from restart input file
# or provided in YAML
def to_int_or_none(val):
  return None if val is None else int(val)

rdate_out = to_int_or_none(config_rest["restart_time"]["rdate_out"])
rhr_out   = to_int_or_none(config_rest["restart_time"]["rhr_out"])
rdate_in  = to_int_or_none(config_rest["restart_time"]["rdate_in"])
rhr_in    = to_int_or_none(config_rest["restart_time"]["rhr_in"])

# Input restart:
# Where original CICE6 restart file is located:
if enmb is None:
  pthrst_in = config_rest["cice_paths"]["pth_in"]
  flrst_in  = config_rest["rest_names"]["flrst_in"]
else:
  pthrst_in = config_rest["cice_mem_paths"]["pth_in"].format(enmb=f"{enmb:03d}")
  flrst_in  = config_rest["rest_names"]["mem"]["flrst_in"]


if flrst_in is None:
  raise RuntimeError("Input restart file name is missing in YAML")
else:
  flrst_in_start  = flrst_in

# Derive restart time from file if not specified in YAML
if rdate_in is None or rhr_in is None:
  with xarray.open_dataset(os.path.join(pthrst_in, flrst_in)) as ds:
    try:
      year  = ds.attrs["myear"]
      month = ds.attrs["mmonth"]
      day   = ds.attrs["mday"]
      sec   = ds.attrs["msec"]
    except KeyError as err:
      raise RuntimeError(f"Missing expected attribute in restart file: {err}")

    rdate_in = int(year*10000 + month*100 + day)
    rhr_in   = sec // 3600

if rdate_out is None:
  rdate_out = rdate_in
if rhr_out is None:
  rhr_out = rhr_in

# Output restart:
#pthrst_out = config_rest["cice_paths"]["pth_out"].format(regn=regn)
#flrst_out  = config_rest["rest_names"]["flrst_out"]
if enmb is None:
  pthrst_out = config_rest["cice_paths"]["pth_out"]
  flrst_out = config_rest["rest_names"]["flrst_out"]
else:
  pthrst_out = config_rest["cice_mem_paths"]["pth_out"].format(enmb=f"{enmb:03d}")
  flrst_out = config_rest["rest_names"]["mem"]["flrst_out"]

# If restart out is not specified, check if template has been provided 
# to use in constructing file name:
#if flrst_out_tmp is None:
#  flrst_out_tmp = 'cice_restart.YYYYMMDD.HH'  # default
if enmb is None:
  flrst_out_tmp = config_rest["rest_names"]["flrst_tmp"]
else:
  flrst_out_tmp = config_rest["rest_names"]["mem"]["flrst_tmp"]

flrst_out_start = flrst_out

# Temporary file names passed along the processes with suffixes being added
# after inserting fields .iconc or .hsnow etc.
# construct temporary file name 
# from restart input if provided otherwise use default file name
# use template file name if provided:
if flrst_out_tmp is not None:
  fltmp_base = mc6util.change_base_template(flrst_out_tmp, rdate_out, rhr_out, flrst_in)
else:
  # Template not provided, construct from input restart but replace date:
  if flrst_in_start is not None:
    flrst_time_new = mc6util.flname_replace_date(flrst_in_start, rdate_out, rhr_out)
    fltmp_base = os.path.splitext(flrst_time_new)[0]
  else: 
    fltmp_base = f'cice_restart.{rdate_out}.{rhr_out:02d}'  # default


os.makedirs(pthrst_out, exist_ok=True)

print(f"Restart input  directory:\n  {pthrst_in}")
print(f"Restart input  file:\n  {flrst_in_start}")
print(f"Restart output directory:\n  {pthrst_out}")
print(f"Restart output file (if None, will be constructed using {fltmp_base}):\n  {flrst_out_start}")
print(" ==== START INSERTION ====\n\n")

# Insert ice concentration and ice thickness if defined:
if ithkn == 1 or iconc == 1:
  # Create temporary restart if other changes are following:
  # Use it as the final name if flrst_out is not provided
  flrst_out = f'{fltmp_base}.iconc.nc'
  if ithkn == 1:
    flrst_out = f'{fltmp_base}.iconc_ithkn.nc'

  if hsnow == 0 and snitd == 0:
    # No other restart updates, final output name:
    if flrst_out_start is not None:
      flrst_out = flrst_out_start 
      
  print(f"Driver: running ithkn iconc --> CICE6 restart")
  print(f"{flrst_in} --> {flrst_out}")

  cmd_iconc = [
    "python", "insert_iconc_ithkn_cice6rest_global.py",
    "--rdate", str(rdate_in),
    "--rhr", str(rhr_in),
    "--rdate_out", str(rdate_out),
    "--rhr_out", str(rhr_out),
    "--pth_in", str(pthrst_in),
    "--flrst_in", str(flrst_in),
    "--pth_out", str(pthrst_out),
    "--flrst_out", str(flrst_out),
    "--ithkn", str(ithkn),
    "--regn", str(regn),
    "--fyaml", str(fyaml),
  ]

  subprocess.run(cmd_iconc, check=True)

# Insert snow depth:
if hsnow == 1:
  # Update input/output restart files:
  # Output directory keep the same
  if ithkn == 1 or iconc == 1:
    # Continued, use previous restart names:
    flrst_in   = flrst_out
    pthrst_in  = pthrst_out 
    fltmp_base = os.path.splitext(flrst_out)[0]
  else:
    flrst_in = flrst_in_start

  flrst_out = f'{fltmp_base}.hsnow.nc'
  if snitd == 0:
    # Final update:
    if flrst_out_start is not None:
      flrst_out = flrst_out_start
   
  print(f"Driver: inserting hsnow --> CICE6 restart")
  print(f"{flrst_in} --> {flrst_out}")

  cmd_hsnow = [
    "python", "insert_hsnow_cice6rest_global.py",
    "--rdate", str(rdate_in),
    "--rhr", str(rhr_in),
    "--rdate_out", str(rdate_out),
    "--rhr_out", str(rhr_out),
    "--pth_in", str(pthrst_in),
    "--flrst_in", str(flrst_in),
    "--pth_out", str(pthrst_out),
    "--flrst_out", str(flrst_out),
    "--regn", str(regn),
    "--fyaml", str(fyaml),
  ]
    
  subprocess.run(cmd_hsnow, check=True)

if snitd == 1:
  # Need to turn on flags in ice_in to use ITS and snow physics in CICE6
  if ithkn == 1 or iconc == 1 or hsnow == 1:
    # Continued:
    flrst_in = flrst_out
    pthrst_in = pthrst_out
    fltmp_base = os.path.splitext(flrst_out)[0]
  else:
    flrst_in = flrst_in_start
    
  flrst_out = f"{fltmp_base}.snphys.nc"
  # Final update:
  if flrst_out_start is not None:
    flrst_out = flrst_out_start

  print(f"Driver: adding snow ITD and physics variables --> CICE6 restart")
  print(f"{flrst_in} --> {flrst_out}")

  cmd_snitd = [
    "python", "add_snowvar_snowphys_cice6rest_global.py",
    "--rdate", str(rdate_in),
    "--rhr", str(rhr_in),
    "--rdate_out", str(rdate_out),
    "--rhr_out", str(rhr_out),
    "--pth_in", str(pthrst_in),
    "--flrst_in", str(flrst_in),
    "--pth_out", str(pthrst_out),
    "--flrst_out", str(flrst_out),
    "--flrst_out", str(flrst_out),
  ]

  subprocess.run(cmd_snitd, check=True)


if sstmom == 1: 
  # Adjust MOM6 SST  under sea ice
  nsec = int(rhr_in * 3600)
  flmom_in = f"{rdate_in}.{nsec:06d}.MOM.res.nc"
  print(f"Running MOM6 SST update, YAML file = {fyaml}, MOM orig restart = {flmom_in}")

  cmd_sst = [
    "python", "correct_surfT_mom6_restart.py",
    "--flmom_in", str(flmom_in),
    "--flice", str(flrst_out),
    "--fyaml", str(fyaml)
  ]

  subprocess.run(cmd_sst, check=True)



