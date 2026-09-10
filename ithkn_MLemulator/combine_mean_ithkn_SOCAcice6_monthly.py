"""
  from GDAS / SOCA sea ice fields
  Combine temporary mean ice volume and ice thickness / over sea ice 
  daily fields into monthly

  called from derive_monthly_meanithkn.sh

  after deriving daily values:
  Daily: calc_mean_ithkn_SOCAcice6.py
  The script is called from 
  derive_monthly_meanithkn.sh

"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from yaml import safe_load
import argparse
from scipy.interpolate import interp1d

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
import mod_mom6 as mmom6
from mod_mom6 import dx_dy

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--yr", help="Year of ice averaging", type=int, required=True)
parser.add_argument("--mm", 
  help="Month of ice averaging", 
  type=int, 
  nargs="+",
  required=True)
args  = parser.parse_args()

regn  = args.regn
YR    = args.yr
MMl   = args.mm

fyaml = 'paths_ML.yaml'
with open(fyaml) as ff:
  pths_ml = safe_load(ff)

# Load daily mean ice vol and thkn
pthprd  = os.path.join(pths_ml["PRED"]["pthprd"],'tmp')

for MM in MMl:
  prefix = f"cice6_mean_ivol_ithkn_{YR}{MM:02d}"
  suffix = f"_{regn}.npz"

  FLS = sorted(
      os.path.join(pthprd, f)
      for f in os.listdir(pthprd)
      if f.startswith(prefix) and f.endswith(suffix)
  )

  assert len(FLS) > 3, (
      f"Not enough daily ithkn files for {YR}/{MM:02d} "
      f"in {pthprd}: found {len(FLS)}"
  )

  IVOL = []
  ITHKN = []
  for dfls in FLS:
    print(f"Loading tmp file: {dfls}")
    tmp = np.load(dfls)
    ivol = tmp["IVOL"]
    ithkn = tmp["ITHKN"]

    IVOL.append(ivol)
    ITHKN.append(ithkn)

  IVOL = np.asarray(IVOL)
  ITHKN = np.asarray(ITHKN)

  ivol_mn = np.mean(IVOL)
  ithkn_mn = np.mean(ITHKN)
    
  print(f"{YR}/{MM:02d}: Monthly ice vol = {ivol_mn*1e-9:.3f} km3, ice thkn = {ithkn_mn:.2f}m")

  # dump monthly values:
  pthimn = os.path.join(pths_ml["PRED"]["pthprd"],'mean_ithkn')
  os.makedirs(pthimn, exist_ok=True)
  fltmp   = f"cice6_mean_ithkn_{YR}{MM:02d}_{regn}.npz"
  dfliceout = os.path.join(pthimn, fltmp)

  print(f"Saving monthly mean ice vol & thickness  --> {dfliceout}")
  np.savez(dfliceout, IVOL=ivol_mn, ITHKN=ithkn_mn)



