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

rest_date = 20250103
rest_hr   = 0
rhr_out   = rest_hr
flrst_in  = None
flrst_out = None
iconc  = 1  # insert iconc from NRT NSIDC
ithkn  = 0  # insert ithkn clim
hsnow  = 0  # insert snow depth clim
snphys = 0  # snow physics on / off
snitd  = 0  # snow redistribution over ice 

parser = argparse.ArgumentParser()
parser.add_argument("--rdate_in", help=f"restart date input file, default={rest_date}", type=int)
parser.add_argument("--rhr_in", help=f"input restart hour = 0, ..., 23, default={rest_hr}", type=int)
parser.add_argument("--rdate_out", help=f"output restart date if different from {rest_date}", type=int)
parser.add_argument("--rhr_out", help=f"output restart hour if date is different from {rest_hr}", type=int)
parser.add_argument("--pth_in", help="input restart directory with original file", type=str)
parser.add_argument("--pth_out", help="output restart directory where new file be dumped", type=str)
parser.add_argument("--flrst_in", help=f"rest file in, otherwise name constructed from {rest_date}", type=str)
parser.add_argument("--flrst_out", help=f"new rest file, otherwise name constr. from rdate_out", type=str)
parser.add_argument("--iconc", type=int, 
                    choices=[0,1], help=f"insert ice conc NRT NSIDC, default={iconc}")
parser.add_argument("--ithkn", type=int, 
                    choices=[0,1], help=f"insert ice thickn climatology, default={ithkn}")
parser.add_argument("--hsnow", type=int, 
                    choices=[0,1], help=f"insert snow depth climatology, default={hsnow}")
#parser.add_argument("--snphys", type=int, 
#                    choices=[0,1], help=f"use snow physics/metamorphism, default={snphys}")
parser.add_argument("--snitd", type=int, 
                    choices=[0,1], help=f"use snow redistribution over ice, default={snitd}")
parser.add_argument("--regn", help=f"region where restart is being updated", 
                    choices=['south','north','global'], required=True, type=str)
args = parser.parse_args()

regn      = args.regn
rdate_in  = args.rdate_in  if args.rdate_in  is not None else rest_date
rhr_in    = args.rhr_in    if args.rhr_in    is not None else rest_hr
rdate_out = args.rdate_out if args.rdate_out is not None else rest_date
rhr_out   = args.rhr_out   if args.rhr_out   is not None else rhr_in
flrst_in  = args.flrst_in  if args.flrst_in  else None
flrst_out = args.flrst_out if args.flrst_out else None
iconc     = args.iconc     if args.iconc     is not None else iconc
ithkn     = args.ithkn     if args.ithkn     is not None else ithkn
hsnow     = args.hsnow     if args.hsnow     is not None else hsnow
snitd     = args.snitd     if args.snitd     is not None else snitd
pth_in    = args.pth_in    if args.pth_in    else None
pth_out   = args.pth_out   if args.pth_out   else None
#snphys    = args.snphys    if args.snphys    is not None else snphys

syst_info = os.uname()
machine = syst_info.nodename

if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
elif 'an' in machine:
  print("Running on PPAN node:", machine)
  node_nm = "ppan"
else:
  print("Unknown machine:", machine)

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

# Input restart:
# Where original CICE6 restart file is located:
if pth_in is None:
  pthrst_in = os.path.join(pths_ufs[node_nm]["MOM6"]["pthrest"],'new')
else:
  pthrst_in = pth_in

if flrst_in is None:
  dnmbR = mtime.rdate2datenum(rdate_in*100 + rhr_in)
  yrR, mmR, ddR, hrR = mtime.datevec(dnmbR, round_hrs=True)[:4]
  nsecR = hrR*3600
  flrst_in_start = f"cice_model.res.{yrR}{mmR:02d}{ddR:02d}.{nsecR:06d}.nc"
else:
  flrst_in_start  = flrst_in

# Output restart:
if pth_out is None:
  pthrst_out = os.path.join(pths_ufs[node_nm]["MOM6"]["pthrest"],f'cice6_{regn}')
else:
  pthrst_out = pth_out

os.makedirs(pthrst_out, exist_ok=True)

flrst_out_start = flrst_out

# Temporary file names passed along the processes
# construct from restart input if provided otherwise use template file name
if flrst_in is not None:
  fltmp_base = os.path.splitext(flrst_in_start)[0]
else: 
  fltmp_base = f'cice_restart.{rdate_out}.{rhr_out:02d}'

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
    "--ithkn", str(ithkn),
    "--regn", str(regn),
  ]

  if flrst_out is not None:
    cmd_iconc += ["--flrst_out", flrst_out]

  subprocess.run(cmd_iconc, check=True)

# Insert snow depth:
if hsnow == 1:
  # Update input/output restart files:
  # Output directory keep the same
  if ithkn == 1 or iconc == 1:
    # Continued:
    flrst_in = flrst_out
    pthrst_in = pthrst_out 
    fltmp_base = os.path.splitext(flrst_out)[0]
  else:
    flrst_in = flrst_in_start

  flrst_out = f'{fltmp_base}.hsnow.nc'
  if snitd == 0:
    # Final update:
    if flrst_out_start is not None:
      flrst_out = flrst_out_start
   
  print(f"Driver: running hsnow --> CICE6 restart")
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
    "--regn", str(regn),
  ]
    
  if flrst_out is not None:
    cmd_hsnow += ["--flrst_out", flrst_out]

  subprocess.run(cmd_hsnow, check=True)

if snitd == 1:
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
    "--pth_in", str(pthrst_in),
    "--flrst_in", str(flrst_in),
    "--pth_out", str(pthrst_out),
    "--flrst_out", str(flrst_out),
  ]

  if flrst_out is not None:
    cmd_snitd += ["--flrst_out", flrst_out]

  subprocess.run(cmd_snitd, check=True)


 
