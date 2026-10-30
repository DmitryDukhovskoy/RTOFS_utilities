"""
  Driver to prepare CICE6 restart for specified
  setup of initial conditions

  Updated version:
  allows only ithkn insertion without iconc
  parameters specified in yaml 

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

# SPecify all paths, files, parameters in yaml
#fyaml = 'cice6rest_files_SFS_ML.yaml'
parser = argparse.ArgumentParser()
parser.add_argument(
  "--regn",
  help=f"region where restart is being updated, default=global", 
  choices=['south','north','global'],
  default="global", 
)
parser.add_argument(
  "--fyaml", 
  help=f"YAML with local directories, filenames, restart dates",
  required=True,
)
args = parser.parse_args()

print(f"Reading YAML with restart info: {args.fyaml}\n")
with open(args.fyaml) as ff:
  config_rest = safe_load(ff)

PARAMS = {
  "regn"  : args.regn,
  "fyaml" : args.fyaml,
  "iconc" : config_rest["params"]["iconc"],
  "ithkn" : config_rest["params"]["ithkn"],
  "hsnow" : config_rest["params"]["hsnow"],
  "snitd" : config_rest["params"]["snitd"],
  "enmb"  : config_rest["params"]["enmb"], 
  "sstmom": config_rest["params"]["sstmom"],
}

print("Input parameters")
for name, value in PARAMS.items():
  if isinstance(value, (int, float)):
    print(f"{name}  : {value}")

if all(value == 0 for value in PARAMS.values()
       if isinstance(value, (int, float))):
  print("All input numeric parameters are zero, no insertion specified, stopping ...") 
  sys.exit(0) 

regn   = PARAMS["regn"]
iconc  = PARAMS["iconc"]
ithkn  = PARAMS["ithkn"]
hsnow  = PARAMS["hsnow"]
snitd  = PARAMS["snitd"]
fyaml  = PARAMS["fyaml"]
enmb   = PARAMS["enmb"]
sstmom = PARAMS["sstmom"]

ins_iconc = iconc == 1
ins_ithkn = ithkn == 1
ins_hsnow = hsnow == 1
ins_snitd = snitd == 1
ins_sstmom = sstmom == 1

# Output Restart dates if missing - same as input
# input dates are deduced from restart input file
# or provided in YAML
def to_int_or_none(val):
  return None if val is None else int(val)

rdate_out = to_int_or_none(config_rest["restart_time"]["rdate_out"])
rhr_out   = to_int_or_none(config_rest["restart_time"]["rhr_out"])
rdate_in  = to_int_or_none(config_rest["restart_time"]["rdate_in"])
rhr_in    = to_int_or_none(config_rest["restart_time"]["rhr_in"])

# Original input restart
# Where original CICE6 restart file is located:
if enmb is None:
  pthrst_in = config_rest["cice_paths"]["pth_in"]
  flrst_in  = config_rest["rest_names"]["flrst_in"]
else:
  pthrst_in = config_rest["cice_mem_paths"]["pth_in"].format(enmb=f"{enmb:03d}")
  flrst_in  = config_rest["rest_names"]["mem"]["flrst_in"]

if flrst_in is None:
  raise RuntimeError("Input restart file name is missing in YAML")

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
      raise RuntimeError(f"Missing expected time attributes in restart file: {err}")

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
if enmb is None:
  flrst_out_tmp = config_rest["rest_names"]["flrst_tmp"]
else:
  flrst_out_tmp = config_rest["rest_names"]["mem"]["flrst_tmp"]

# Keep requested final output restart name if specified in yaml:
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
  flrst_time_new = mc6util.flname_replace_date(flrst_in_start, rdate_out, rhr_out)
  fltmp_base = os.path.splitext(flrst_time_new)[0]


os.makedirs(pthrst_out, exist_ok=True)

# Restart tracking:
restart_in_dir  = pthrst_in
restart_in_file = flrst_in_start

print(f"Restart input  directory:\n  {restart_in_dir}")
print(f"Restart input  file:\n  {restart_in_file}")
print(f"Restart output directory:\n  {pthrst_out}")
print(
    "Configured final output file "
    f"(None means constructed from {fltmp_base}):\n"
    f"  {flrst_out_start}"
)
print(" ==== START INSERTION ====\n\n")


# Insert ice concentration and ice thickness if requested:
if ins_ithkn or ins_iconc:
  # Create temporary restart if other changes are following:
  if ins_iconc and not ins_ithkn:
    flrst_out = f'{fltmp_base}.iconc.nc'
  elif ins_ithkn and not ins_iconc:
    flrst_out = f'{fltmp_base}.ithkn.nc'
  else:
    flrst_out = f'{fltmp_base}.iconc_ithkn.nc'

  # If no further CICE restart updates are needed,
  # use the configured final output filename if provided.
  if not ins_hsnow and not ins_snitd:
    if flrst_out_start is not None:
      flrst_out = flrst_out_start 
      
  print(f"Driver: running ithkn iconc --> CICE6 restart")
  print(f"{restart_in_file} --> {flrst_out}")

  cmd_iconc = [
    "python", "insert_iconc_ithkn_cice6rest_global.py",
    "--rdate", str(rdate_in),
    "--rhr", str(rhr_in),
    "--rdate_out", str(rdate_out),
    "--rhr_out", str(rhr_out),
    "--pth_in", str(restart_in_dir),
    "--flrst_in", str(restart_in_file),
    "--pth_out", str(pthrst_out),
    "--flrst_out", str(flrst_out),
    "--iconc", str(iconc),
    "--ithkn", str(ithkn),
    "--regn", str(regn),
    "--fyaml", str(fyaml),
  ]

  subprocess.run(cmd_iconc, check=True)

  # Update current restart after successful insertio
  restart_in_dir  = pthrst_out
  restart_in_file = flrst_out
  fltmp_base = os.path.splitext(flrst_out)[0]


# Insert snow depth:
if ins_hsnow:
  # Construct the next output from the current base name.
  flrst_out = f"{fltmp_base}.hsnow.nc"

  # If snow physics is not the next operation, this is
  # the final CICE update.
  if not ins_snitd:
    if flrst_out_start is not None:
        flrst_out = flrst_out_start

  print("Driver: inserting hsnow --> CICE6 restart")
  print(f"{restart_in_file} --> {flrst_out}")

  cmd_hsnow = [
      "python",
      "insert_hsnow_cice6rest_global.py",
      "--rdate", str(rdate_in),
      "--rhr", str(rhr_in),
      "--rdate_out", str(rdate_out),
      "--rhr_out", str(rhr_out),
      "--pth_in", str(restart_in_dir),
      "--flrst_in", str(restart_in_file),
      "--pth_out", str(pthrst_out),
      "--flrst_out", str(flrst_out),
      "--regn", str(regn),
      "--fyaml", str(fyaml),
  ]

  subprocess.run(cmd_hsnow, check=True)

  # Update current restart after successful insertion.
  restart_in_dir  = pthrst_out
  restart_in_file = flrst_out
  fltmp_base = os.path.splitext(flrst_out)[0]

# Add snow ITD and snow physics variables
# Need to turn on flags in ice_in to use ITD and snow physics in CICE6
if ins_snitd:
  flrst_out = f"{fltmp_base}.snphys.nc"

  # This is the final CICE update.
  if flrst_out_start is not None:
    flrst_out = flrst_out_start

  print(
      "\nDriver: adding snow ITD and physics variables "
      "--> CICE6 restart"
  )
  print(f"{restart_in_file} --> {flrst_out}")

  cmd_snitd = [
    "python",
    "add_snowvar_snowphys_cice6rest_global.py",
    "--rdate", str(rdate_in),
    "--rhr", str(rhr_in),
    "--rdate_out", str(rdate_out),
    "--rhr_out", str(rhr_out),
    "--pth_in", str(restart_in_dir),
    "--flrst_in", str(restart_in_file),
    "--pth_out", str(pthrst_out),
    "--flrst_out", str(flrst_out),
  ]

  subprocess.run(cmd_snitd, check=True)

  # Update current restart after successful insertion.
  restart_in_dir  = pthrst_out
  restart_in_file = flrst_out
  fltmp_base = os.path.splitext(flrst_out)[0]


if ins_sstmom: 
  # Adjust MOM6 SST  under sea ice

  flmom_in = config_rest["mom6_names"]["flrst_in"]
  if flmom_in is None:
    rhr_mom_in = config_rest["mom6_names"]["rhr_in"]
    nsec = int(rhr_mom_in * 3600)
    flmom_in = f"{rdate_in}.{nsec:06d}.MOM.res.nc"

  print(f"Running MOM6 SST update, YAML file = {fyaml}, MOM orig restart = {flmom_in}")

  cmd_sst = [
    "python", "correct_surfT_mom6_restart.py",
    "--flmom_in", str(flmom_in),
    "--flice", str(restart_in_file),
    "--rdate_out", str(rdate_out),
    "--fyaml", str(fyaml)
  ]

  subprocess.run(cmd_sst, check=True)

print("\n ==== ALL REQUESTED UPDATES COMPLETED ====\n")
print(f"Latest CICE restart directory: {restart_in_dir}")
print(f"Latest CICE restart file:      {restart_in_file}")

