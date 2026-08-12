"""
  Derive dynamic predictor: sea ice divergence averaged over some area
  GLORYS

  Two approaches are possible (using divergence theorem)
  - average spatial integral of div(U) * dA
  - average integra of U flux across the boundary: U*n*dl, n - is normal comp. to the segment
  boundary flux is prefereable (as it conserved divergence)
  1st approach is straight forward but not exact

  Use time stamps and J,I grid points from the response
  variable (ithkn), that need to be run first


  For faster processing, run unstaging script before this:
  /home/Dmitry.Dukhovskoy/scripts/GLORYS_anls/unstage_glorys.sh

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray as xr
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, cm
from yaml import safe_load
import argparse
from pathlib import Path

#ROOT = Path(__file__).resolve().parent

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
import mod_time as mtime
import mod_glorys as mglr 
from mod_mom6 import dx_dy

#from MyPython.mod_cice6_utils import change_base_template, flname_replace_date

parser = argparse.ArgumentParser()
parser.add_argument("--dxy", 
               help=f"Min dist (km) between data points (~corr.scale), to skip close i,j points", 
               default=50,
               type=int)
parser.add_argument("--ys", help="Year start, default=1993", default=1993, type=int)
parser.add_argument("--ye", help="Year end, default=2025", default=2025, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'], 
                    required=True, type=str)
parser.add_argument("--load", help="Load saved sst tmp file, continue from last record (1), start from time 0 (0)", 
                  choices=[0,1], required=True, type=int)
args = parser.parse_args()

dxy   = args.dxy    
YS    = args.ys
YE    = args.ye
regn  = args.regn
load_saved = args.load == 1

fld_name = 'divU'

regions = {
    "north": ("Arctic", 65.0),
    "south": ("Antarctic", -60.0),
}
regn_name, lat0 = regions[regn]
varnmu = 'usi'
varnmv = 'vsi'

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

DIRS = {
  "pthithkn" : config_predictor["linregr"]["pthithkn"],
  "pthiconc" : config_predictor["linregr"]["pthiconc"],
  "pthsst"   : config_predictor["linregr"]["pthsst"],
  "pthssh"   : config_predictor["linregr"]["pthssh"],
  "pthui"    : config_predictor["linregr"]["pthui"],
  "pthvi"    : config_predictor["linregr"]["pthvi"],
  "ptht2m"   : config_predictor["linregr"]["ptht2m"].format(regn_name=regn_name),
  "pthout"   : config_predictor["linregr"]["pthout"],
  "ithkntmp" : config_predictor["linregr"]["ithkntmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "iconctmp" : config_predictor["linregr"]["iconctmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "ssttmp"   : config_predictor["linregr"]["ssttmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "divutmp"  : config_predictor["linregr"]["divutmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  }

# Read time array ad J,I sample grid points:
pthout = DIRS["pthout"]
fltmp = DIRS["ithkntmp"]
dflithkn = os.path.join(pthout, fltmp)
  
print(f"Loading saved {dflithkn}, will start from last saved record")
if not os.path.isfile(dflithkn):
  print(f"Missing tmp file {dflithkn}\n  start from time = 0")
else:
  data = np.load(dflithkn)
  JG = data["JG"]
  IG = data["IG"]
  DNMB = data["DNMB"]
  nrecs = len(DNMB)

  #Find last saved record, no nans in the column:
  #processed = np.all(np.isfinite(YY), axis=0)
  #Nproc = np.count_nonzero(processed)
  print(f"Number of time records = {nrecs}")


# Read GLORYS grid:
pthice = os.path.join(DIRS["pthsst"],f"{YS}")
rdate = f"{YS*10000+100+1}"

# Find file:
dflglr = mglr.find_file(rdate, pthice)

with xr.open_dataset(dflglr) as dsice:
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values

# 0 <= lon < 360
LON = (LON + 360) % 360
hlon, hlat = np.meshgrid(LON, LAT)
DX, DY = dx_dy(hlon, hlat)
Acell = DX*DY

def calc_divU(dltI, dltJ, ii, jj, U2d, V2d, Acell, DX, DY):
  """
    Compute average div u ice over specified region
    Assuming output fieds are at the cell centers
  """
  divU = 0
  jdm, idm = U2d.shape

  Intgr = 0.0
  # Define the box around the grid point:
  iS = int(ii - dltI)
  iE = int(ii + dltI)
  iS = np.max([iS, 0])
  iE = np.min([iE, idm-1])

  jS = int(jj - dltJ)
  jE = int(jj + dltJ)
  # Boundaries, better - for global grid
  # use grid points at the opposite side
  jS = np.max([0, jS])
  jE = np.min([jE, jdm-1])

  # Integrate along the boundary:
  # Note different sign of outward norm vectors 
  # along the box sides
  Intgr = (
    np.sum(U2d[jS:jE+1, iE] * DY[jS:jE+1, iE])
    - np.sum(U2d[jS:jE+1, iS] * DY[jS:jE+1, iS])
    + np.sum(V2d[jE, iS:iE+1] * DX[jE, iS:iE+1])
    - np.sum(V2d[jS, iS:iE+1] * DX[jS, iS:iE+1])
  )

  # subtract half of the four corner contributions
  Intgr -= 0.5 * (
    U2d[jS, iE] * DY[jS, iE]
    + U2d[jE, iE] * DY[jE, iE]
    - U2d[jS, iS] * DY[jS, iS]
    - U2d[jE, iS] * DY[jE, iS]
    + V2d[jE, iS] * DX[jE, iS]
    + V2d[jE, iE] * DX[jE, iE]
    - V2d[jS, iS] * DX[jS, iS]
    - V2d[jS, iE] * DX[jS, iE]
  )

  # Space-Average divergence:
  Area = np.sum(Acell[jS:jE+1, iS:iE+1])
  assert Area > 0, f"ii={ii}, jj={jj}, dltI={dltI}, dltJ={dltJ}, Area = {Area}"
  div_uice = Intgr / Area

  return div_uice


# Construct predictor sst time series for all locations, 
# Or load previously saved
pthout = DIRS["pthout"]
fltmp = DIRS["divutmp"]
dfltmp = os.path.join(pthout, fltmp)

irec_start = 0
YY = None
if load_saved:
  print(f"Loading saved {dfltmp}, will start from last saved record")
  if not os.path.isfile(dfltmp):
    print(f"Missing tmp file {dfltmp}\n  start from time = 0")
  else:
    # Do not load saved JG, IG from this file: 
    # Keep using those saved in ithkn - should be identical
    data = np.load(dfltmp)
    YY = data["YY"]
    DNMB_check = data["DNMB"]  

    # Check that this is the right time series:
    dtmp = np.floor(np.abs(DNMB - DNMB_check))
    assert np.max(dtmp) == 0, "Check DNMB - dates do not match with saved time series"

    #Find last saved record, no nans in the column:
    processed = np.all(np.isfinite(YY), axis=0)
    irec_start = np.count_nonzero(processed)
    print(f"Next record to start {irec_start}")

#YY, JG, IG = construct_ifld(DNMB, DIRS, JG, IG, dfltmp, irec_start, YY, Acell, DX, DY, dxy, dump_tstp=10)

"""
  Construct time series of response variable (ithkn)
  2D: locations x time
"""
npnts = len(JG)
nrecs = len(DNMB)
dump_tstp = 10
if YY is None:
  YY = np.zeros((npnts, nrecs), dtype=float)*np.nan
for irec, dnmb0 in enumerate(DNMB):
  YR, MM, DD = mtime.datevec(dnmb0)[:3]

  if irec < irec_start:
    print(f"Skipping, already processed: {YR}/{MM:02d}/{DD:02d}")
    continue
  
  print(f"Reading {fld_name} {YR}/{MM:02d}/{DD:02d}")

  rdate = int(YR*1e4 + MM*100 + DD)

  pthu = os.path.join(DIRS["pthui"],f"{YR}")
  dflice = mglr.find_file(rdate, pthu)
  with xr.open_dataset(dflice) as dsice:
    A2d = dsice['usi'].isel(time=0).values.squeeze()

  # Treat nans as no ice grid cells
  U2d = np.nan_to_num(A2d, nan=0.0)

  pthv = os.path.join(DIRS["pthvi"],f"{YR}")
  dflice = mglr.find_file(rdate, pthv)
  with xr.open_dataset(dflice) as dsice:
    A2d = dsice['vsi'].isel(time=0).values.squeeze()

  # Treat nans as no ice grid cells
  V2d = np.nan_to_num(A2d, nan=0.0)

  fld_pnts = []
  for jj, ii in zip(JG, IG):
    # Estimate box size based on min distance criterion
    dltX = DX[jj,ii]*1e-3  # km
    dltY = DY[jj,ii]*1e-3  # km
    dltI = int(np.ceil(dxy / dltX))
    dltJ = int(np.ceil(dxy / dltY))
    div_uice = calc_divU(dltI, dltJ, ii, jj, U2d, V2d, Acell, DX, DY)
    fld_pnts.append(div_uice)

  YY[:,irec] = np.asarray(fld_pnts)

  assert YY.shape[0] == len(IG), f"Check YY shape does not match IG {YY.shape}"

  if (irec + 1) % dump_tstp == 0: 
    print(f"TMP step: Saving {fld_name} time series  --> {dfltmp}")
    np.savez(dfltmp,
           YY=YY,
           JG=JG,
           IG=IG,
           DNMB=DNMB)

# Save:
print(f"Final Saving ithkn time series and IG, JG --> {dfltmp}")
np.savez(dfltmp,
       YY=YY,
       JG=JG,
       IG=IG,
       DNMB=DNMB)


f_check = False
if f_check:
  # Land mask:
  LMsk = None
  pthssh = os.path.join(DIRS["pthssh"],f"{YS}")
  dflssh = mglr.find_file(rdate, pthssh)
  with xr.open_dataset(dflssh) as dszos:
    SSH = dszos['zos'].isel(time=0).values.squeeze()

  LMsk = np.where(np.isfinite(SSH),1,0)

  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
  ax1.contour(LMsk, [0.99], linestyles='solid', colors=[(0.5,0.8,1)])
  ax1.scatter(IG, JG, s=5, color=(0.8,0.4,0), marker='.')


