"""
  Calc mean ice volume and ice thickness / over sea ice 
  from GDAS / SOCA sea ice fields
  Process 1 day

  The script is called from 
  derive_monthly_meanithkn.sh

  Fields from HPSS, use scripts to fetch fields
  get_soca_ice_restart.sh

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
parser.add_argument("--rdate", help="Date of ithkn prediction YYYMMDD", type=int, required=True)
args  = parser.parse_args()

regn       = args.regn
rdate      = args.rdate

fyaml = 'paths_ML.yaml'
with open(fyaml) as ff:
  pths_ml = safe_load(ff)

dnmbR = mtime.rdate2datenum(rdate)
YRr, MMr, DDr = mtime.datevec(dnmbR)[:3]

pthdata   = pths_ml["GDAS"]["pthdata"]  # root dir for processed data
pthice    = pths_ml["GDAS"]["pthsoca_ice"].format(YR=YRr, MM=MMr, DD=DDr)
regn_name = pths_ml[regn]["name"]
regn_lat0 = pths_ml[regn]["lat0"]  # bounding lat for region

# Get MOM6 grid:
pthgrid    = pths_ml["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
     
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
DX, DY = dx_dy(hlon, hlat)
Acell = DX*DY

  
with xr.open_dataset(dftopo_mom) as dstopo:
  depth = dstopo['depth'].values.squeeze()
  
# Convert all positive values -> land (100) and ocean (<0):
HH = np.where(depth < 1.e-20, 100., -depth)

jdm, idm = HH.shape
LMsk = HH < 0

# Set mask for region and ML lat bounds
# Consider only interior basins away from marginal ice zone
assert HH.shape == hlat.shape, f"Shape mismatch: HH={HH.shape}, hlat={hlat.shape}"

if regn == 'north':
  lat0 = 65.
  DOMAIN = (hlat > lat0) & (LMsk)
elif regn == 'south':
  lat0 = -60.
  DOMAIN = (hlat < lat0) & (LMsk)


# Read SOCA ice fields:
flice = pths_ml["GDAS"]["flsoca_ice"].format(YR=YRr, MM=MMr, DD=DDr)
dflice = os.path.join(pthice, flice)

assert os.path.isfile(dflice), f"File not found: {dflice}"
with xr.open_dataset(dflice) as dsice:
  aicen = dsice["aicen"].values  # ice conc by cats
  vicen = dsice["vicen"].values   # ice vol by cats m3/m2_ice

# Check shapes / dimensions:
assert aicen.ndim == 3, "aicen n dimensions not 3"
assert vicen.shape == aicen.shape, "Check vicen and aicen shapes"
assert aicen.shape[1:] == HH.shape, "Check aicen grid dim, should match HH"

ncat = aicen.shape[0]

# Derive aggregated ice conc and ice volume in m3/m2_ice
aice = np.nansum(aicen, axis=0)
aice[aice > 1] = 1.
vice_cell = np.nansum(aicen * vicen, axis=0)  # m3/m2_cell

# Derive mean ice thickness over ice covered area:
#hice = np.divide(vice_cell, aice, out=np.zeros_like(vice_cell), where=aice > 0)

# Apply regional mask:
aice = np.nan_to_num(aice, nan=0.0)
aice[~DOMAIN] = np.nan
vice_cell = np.nan_to_num(vice_cell, nan=0.0)
vice_cell[~DOMAIN] = np.nan 

# Ice volume
IVOL = np.nansum(vice_cell * Acell) # m3 of ice

# Total ice-covered area [m2]
iarea = np.nansum(Acell * aice)

# Mean ice thkn over ice-covered area
ITHKN = IVOL / iarea if iarea > 0 else np.nan

print(f"{rdate}: Mean ice volume = {IVOL*1e-9:.3f} km3, ice thickness = {ITHKN:.2f}m")

# Save numpy binary to be combined later
pthprd  = os.path.join(pths_ml["PRED"]["pthprd"],'tmp')
os.makedirs(pthprd, exist_ok=True)
fltmp   = f"cice6_mean_ivol_ithkn_{rdate}_{regn}.npz"
dfliceout = os.path.join(pthprd, fltmp)

print(f"Saving ice vol & thickness to a temporary file --> {dfliceout}")
np.savez(dfliceout, IVOL=IVOL, ITHKN=ITHKN)


