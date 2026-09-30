"""
Version 2

  Derive dynamic predictor: ocean surf temp (model layer 1, ~ 1m thick)
  GLORYS

  Use time stamps and J,I grid points from the response
  variable (ithkn), that need to be run first

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
import mod_icepredict as micpr

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="Region to process", choices=['north','south'], 
                    required=True, type=str)
parser.add_argument(
  "--load",
  help="Load saved sst tmp file, continue from last record (1), start from time 0 (0 default)",
  choices=[0,1],
  default=0,
  type=int
)

args = parser.parse_args()

regn       = args.regn
load_saved = args.load == 1
fld_name   = 'sst'
varnm      = "thetao"

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

# Load parameters:
regn_name = config_predictor["regn"][regn]["name"]
lat0      = config_predictor["regn"][regn]["lat_bnd"]
tstep     = config_predictor["params"]["tstep"]
dxy       = config_predictor["params"]["dxy"]
YS        = config_predictor["params"]["ys"]
YS        = config_predictor["params"]["ys"]
YE        = config_predictor["params"]["ye"]
tstep_era = config_predictor["params"]["tstep"]

intgr_time = 90  # Time for freeze degree days accumulation, back from current time
dump_tstp = 50

DIRS = {
  "pthithkn" : config_predictor["linregr"]["pthithkn"],
  "pthiconc" : config_predictor["linregr"]["pthiconc"],
  "pthsst"   : config_predictor["linregr"]["pthsst"],
  "pthssh"   : config_predictor["linregr"]["pthssh"],
  "ptht2m"   : config_predictor["linregr"]["ptht2m"].format(regn_name=regn_name),
  "pthout"   : config_predictor["linregr"]["pthout"],
  "ithkntmp" : config_predictor["linregr"]["ithkntmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "iconctmp" : config_predictor["linregr"]["iconctmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "ssttmp"   : config_predictor["linregr"]["ssttmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  }

def ocean_Tfreeze(s0):
  Tfrz = -0.0575 * s0 + 1.71e-3 * (s0)**(3/2) / 2 - 2.15e-4*(s0)**2
  return Tfrz

# Bound by the lowest possible Tfrz for expected max surface S in high lats
sst_min = round(ocean_Tfreeze(38), 3)

# Get time steps, grid points:
pthinfo = config_predictor["params"]["pthinfo"]
fltime = config_predictor["params"]["fltime"].format(regn=regn, tstep=tstep)
flij   = config_predictor["params"]["flij"].format(regn=regn)
flgrid = config_predictor["params"]["flgrid"]
dfltime = os.path.join(pthinfo, fltime)
dflij   = os.path.join(pthinfo, flij)
dflgrid = os.path.join(pthinfo, flgrid)

assert os.path.isfile(dfltime), f"Time steps file is missing: {dfltime}"
assert os.path.isfile(dflij), f"Subsample grid points file is missing: {dflij}"

# Read saved date numbers and 
# Read time array ad J,I sample grid points
DNMB = np.load(dfltime)
A = np.load(dflij)
JG = A["JG"]
IG = A["IG"]
npnts = IG.shape
nrecs = len(DNMB)

print(f"N of grid points: {npnts}, N time steps: {nrecs}")

# GLORYS grid:
A = np.load(dflgrid)
hlon = A["LON"]
hlat = A["LAT"]
LMsk = A["LMsk"]

# 0 <= lon < 360
hlon = (hlon + 360) % 360


def construct_ifld(DNMB, DIRS, JG, IG, dfltmp, irec_start, YY, dump_tstp=20):
  """
    Construct time series of response variable (ithkn)
    2D: locations x time
  """
  npnts = len(JG)
  nrecs = len(DNMB)
  if YY is None:
    YY = np.zeros((npnts, nrecs), dtype=float)*np.nan
  for irec, dnmb0 in enumerate(DNMB):
    YR, MM, DD = mtime.datevec(dnmb0)[:3]

    if irec < irec_start:
      print(f"Skipping, already processed: {YR}/{MM:02d}/{DD:02d}")
      continue
    
    print(f"irec={irec}, Reading {fld_name} {YR}/{MM:02d}/{DD:02d}")

    pthice = os.path.join(DIRS["pthsst"],f"{YR}")
    rdate = int(YR*1e4 + MM*100 + DD)
    dflice = mglr.find_file(rdate, pthice)

    with xr.open_dataset(dflice) as dsice:
      A2d = dsice[varnm].isel(time=0, depth=0).values.squeeze()
 
    fld_pnts = A2d[JG,IG]
    
    # Treat nans as no ice grid cells
    #A2d = np.nan_to_num(A2d, nan=0.0)
    # Cup min T to min possible freezing T in high lat
    #assert np.min(fld_pnts) > -2., f"Min SST is too cold: {np.min(fld_pnts)}"
    fld_pnts[fld_pnts < sst_min] = sst_min

    assert np.max(fld_pnts < 30.), f"Max SST is too warm: {np.max(fld_pnts)}"
 
    YY[:,irec] = fld_pnts
  
    if (irec + 1) % dump_tstp == 0: 
      print(f"TMP step: Saving {fld_name} time series  --> {dfltmp}")
      np.savez(dfltmp,
             YY=YY,
             JG=JG,
             IG=IG,
             DNMB=DNMB)

  if irec_start < len(DNMB):
    # No need to save if already everything processed
    print(f"END TMP step: Saving {fld_name} time series and IG, JG --> {dfltmp}")
    np.savez(dfltmp,
           YY=YY,
           JG=JG,
           IG=IG,
           DNMB=DNMB)

  return YY, JG, IG

# Construct predictor sst time series for all locations, 
# Or load previously saved
pthout = DIRS["pthout"]
fltmp = DIRS["ssttmp"]
dfltmp = os.path.join(pthout, fltmp)

irec_start = 0
YY = None
if load_saved:
  print(f"Loading saved {dfltmp}, will start from last saved record")
  if not os.path.isfile(dfltmp):
    print(f"Missing tmp file {dfltmp}\n  start from time = 0")
  else:
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

YY, JG, IG = construct_ifld(DNMB, DIRS, JG, IG, dfltmp, irec_start, YY, dump_tstp=20)
  
# Save:
print(f"Final Saving ithkn time series and IG, JG --> {dfltmp}")
np.savez(dfltmp,
       YY=YY,
       JG=JG,
       IG=IG,
       DNMB=DNMB)


f_check = False
if f_check:
  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
  ax1.contour(LMsk, [0.99], linestyles='solid', colors=[(0.5,0.8,1)])
  ax1.scatter(IG, JG, s=5, color=(0.8,0.4,0), marker='.')


