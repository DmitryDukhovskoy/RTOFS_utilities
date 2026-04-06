"""
  Driver to prepare CICE6 restart for specified
  setup of initial conditions

  Can use any restart adding missing fields to it, e.g.:
run prepare_cice6_restart_driver.py --iconc 0 --ithkn 0 --hsnow 1 --snitd 0 --regn global --rdate_in 20250103 --pth_in {pthrst_in} --flrst_in cice_restart.20250103.00.iconc_ithkn.nc

will grab ice restart from flrst_in (with already inserted iconc and ithkn) and will add
hsnow on top of this

similarly can create from original restart (default name) creating iconc+ithkn then hsnow:
run prepare_cice6_restart_driver.py --iconc 1 --ithkn 1 --hsnow 1 --snitd 0 --regn global --rdate_in 20250103 --pth_in {pthrst_in} 


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
    os.path.join(PPTHN, 'MyPython', 'draw_map'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])

import mod_misc1 as mmisc
import mod_cice6_utils as mc6util
import mod_time as mtime

#rest_date = 20250103
rest_date = 20250704
rest_hr   = 0
rhr_out   = rest_hr
flrst_in  = None
flrst_out = None
#fyaml = 'cice6rest_files.yaml'   # YAML with paths, files, dates for restart
fyaml = 'cice6rest_files_SFS.yaml'
iconc  = 1  # insert iconc from NRT NSIDC
ithkn  = 0  # insert ithkn clim
hsnow  = 0  # insert snow depth clim
snphys = 0  # snow physics on / off
snitd  = 0  # snow redistribution over ice 

parser = argparse.ArgumentParser()
parser.add_argument("--iconc", type=int, 
                    choices=[0,1], help=f"insert ice conc NRT NSIDC, default={iconc}")
parser.add_argument("--ithkn", type=int, 
                    choices=[0,1], help=f"insert ice thickn climatology, default={ithkn}")
parser.add_argument("--hsnow", type=int, 
                    choices=[0,1], help=f"insert snow depth climatology, default={hsnow}")
parser.add_argument("--snitd", type=int, 
                    choices=[0,1], help=f"use snow redistribution over ice, default={snitd}")
parser.add_argument("--regn", help=f"region where restart is being updated", 
                    choices=['south','north','global'], required=True, type=str)
parser.add_argument("--fyaml", 
                    help=f"YAML with local directories, filenames, restart dates, default={fyaml}",
                    default=fyaml,
                    type=str)
args = parser.parse_args()

regn      = args.regn
iconc     = args.iconc  if args.iconc  is not None else iconc
ithkn     = args.ithkn  if args.ithkn  is not None else ithkn
hsnow     = args.hsnow  if args.hsnow  is not None else hsnow
snitd     = args.snitd  if args.snitd  is not None else snitd
fyaml     = args.fyaml  if args.fyaml  is not None else fyaml 

print(f"Reading YAML with restart info: {fyaml}\n")
with open(fyaml) as ff:
  config_rest = safe_load(ff)

# Output Restart dates if missing - same as input
# input dates are deduced from restart input file
rdate_out = config_rest["restart_time"]["rdate_out"] 
rhr_out   = config_rest["restart_time"]["rhr_out"]

# Input restart:
# Where original CICE6 restart file is located:
pthrst_in = config_rest["cice_paths"]["pth_in"]
flrst_in  = config_rest["rest_names"]["flrst_in"]

if flrst_in is None:
  raise RuntimeError("Input restart file name is missing in YAML")
else:
  flrst_in_start  = flrst_in

with xarray.open_dataset(os.path.join(pthrst_in, flrst_in)) as ds:
  year  = ds.attrs["myear"]
  month = ds.attrs["mmonth"]
  day   = ds.attrs["mday"]
  sec   = ds.attrs["msec"]

  rdate_in = int(year*10000 + month*100 + day)
  rhr_in   = sec // 3600

if rdate_out is None:
  rdate_out = rdate_in
if rhr_out is None:
  rhr_out = rhr_in

# Output restart:
pthrst_out = config_rest["cice_paths"]["pth_out"].format(regn=regn)
flrst_out  = config_rest["rest_names"]["flrst_out"]

# If restart out is not specified, check if template has been provided 
# to use in constructing file name:
#if flrst_out_tmp is None:
#  flrst_out_tmp = 'cice_restart.YYYYMMDD.HH'  # default
flrst_out_tmp = config_rest["rest_names"]["flrst_tmp"]

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


 
