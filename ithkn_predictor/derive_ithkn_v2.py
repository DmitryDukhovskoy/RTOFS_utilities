"""
Version 2: for both polar regions
  Grid points and time steps are prepared in define_time_IJpnts.py
  Grid points with no ice cover during multiple years are eliminated


  Derive response variable: ice thickness
  extract all time steps at specified grid points


  Use every N-days of daily data
  and subsample high-res. GLORYS fields: there is lots of spatial
  correlation, no need to do every grid point

  GLORYS ice thicknesses are ice_mean (i.e. sum(hice(k)*aice(k))/sum(aice(k)))
  There are negative values, e.g. 
  sithick_mercatorglorys12v1_gl12_mean_19940621_R19940622.nc
  j0=1771, i0=1852 
  A2d[j0,i0]: -24.217963172122836
  average over previous Nav time steps
  Large values - capped at ithkn_max
  

The daily data is on uda:
/uda/Global_Ocean_Physics_Reanalysis/global/daily/siconc/

/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick

ERA5 T2m fields:
use 1-hr fields on PPAN --> derive daily mean
/archive/uda/ERA5/Hourly_Data_On_Single_Levels/reanalysis/global/1hr-timestep/annual_file-range/Temperature_and_Pressure/T_2m
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
import mod_icepredict as micepr

#from MyPython.mod_cice6_utils import change_base_template, flname_replace_date

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="Region to process", choices=['north','south'], 
                    required=True)
parser.add_argument(
    "--load", 
    help="Load saved ithkntmp, continue from last record (1), start from time 0 (0)", 
    choices=[0,1], 
    required=True, 
    type=int
)
args = parser.parse_args()

regn  = args.regn
load_saved = args.load == 1


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

ithkn_max = 4.  # cap max ice thickness
Ntime_avrg = 3 # for unrealistic (<0) values, use average values from previous records

DIRS = {
  "pthithkn" : config_predictor["linregr"]["pthithkn"],
  "pthiconc" : config_predictor["linregr"]["pthiconc"],
  "pthsst"   : config_predictor["linregr"]["pthsst"],
  "pthssh"   : config_predictor["linregr"]["pthssh"],
  "ptht2m"   : config_predictor["linregr"]["ptht2m"].format(regn_name=regn_name),
  "pthout"   : config_predictor["linregr"]["pthout"],
  "ithkntmp" : config_predictor["linregr"]["ithkntmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  }

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

DNMB = np.load(dfltime)
A = np.load(dflij)
JG = A["JG"]
IG = A["IG"]
npnts = IG.shape
nrecs = len(DNMB)

# GLORYS grid:
A = np.load(dflgrid)
hlon = A["LON"]
hlat = A["LAT"]
LMsk = A["LMsk"]

# 0 <= lon < 360
hlon = (hlon + 360) % 360


# check DNMB array:
check_dnmb = True
if check_dnmb:
  micepr.check_dnmb_array(DNMB, print_months=False)


Npnts = len(IG)
print(f"For dxy={dxy} km and lat0={lat0:.2f}, Selected N pnts={Npnts}")

IG = np.asarray(IG, dtype=int)
JG = np.asarray(JG, dtype=int)

def find_index(j0, i0, JG, IG):
  """
    Find closest index from the sample points
    to j0,i0
  """
  D = np.sqrt((JG-j0)**2 + (IG-i0)**2)
  indx = np.argmin(D)

  return indx

def construct_ithkn(DNMB, DIRS, JG, IG, dfltmp, irec_start, YY, dump_tstp=20):
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
    
    print(f"irec={irec} Reading ithkn {YR}/{MM:02d}/{DD:02d}")

    pthice = os.path.join(DIRS["pthithkn"],f"{YR}")
    rdate = int(YR*1e4 + MM*100 + DD)
    dflice = mglr.find_file(rdate, pthice)

    with xr.open_dataset(dflice) as dsice:
      A2d = dsice['sithick'].isel(time=0).values.squeeze()
 
    # Treat nans as no ice grid cells
    A2d = np.nan_to_num(A2d, nan=0.0)
 
    fld_pnts = A2d[JG,IG]
    #assert np.min(fld_pnts) >= 0, f"Found negative thicknesses {np.min(fld_pnts)}"
    if np.min(fld_pnts) < 0:
      Ineg = np.where(fld_pnts < 0)[0]
      print(f"Found {len(Ineg)} negative thicknesses {np.min(fld_pnts)}")
      for ineg in Ineg:
        j0 = JG[ineg]
        i0 = IG[ineg]
        havrg = np.nanmean(YY[ineg,irec-Ntime_avrg:irec])
        assert np.isfinite(havrg), "havrg is not finite" 
        print(f"i={i0} j={j0}, Replacing negative ithkn {fld_pnts[ineg]:.6f} with mean {havrg:.6f}")
        fld_pnts[ineg] = havrg 
             
    YY[:,irec] = fld_pnts
  
    if (irec + 1) % dump_tstp == 0: 
      print(f"TMP step: Saving ithkn time series and IG, JG --> {dfltmp}")
      np.savez(dfltmp,
             YY=YY,
             JG=JG,
             IG=IG,
             DNMB=DNMB)

  if irec_start < len(DNMB):
    # No need to save if already everything processed
    print(f"END TMP step: Saving ithkn time series and IG, JG --> {dfltmp}")
    np.savez(dfltmp,
           YY=YY,
           JG=JG,
           IG=IG,
           DNMB=DNMB)

  return YY, JG, IG

# Construct response (ithkn) time series concatenating all locations, 
# locations with 0 ithkn will be eliminated from IG, JG
# Or load previously saved
pthout = DIRS["pthout"]
fltmp = DIRS["ithkntmp"]
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
    JG = data["JG"]
    IG = data["IG"]
    DNMB_check = data["DNMB"]  

    # Check that this is the right time series:
    dtmp = np.floor(np.abs(DNMB - DNMB_check))
    assert np.max(dtmp) == 0, "Check DNMB - dates do not match with saved time series"

    #Find last saved record, no nans in the column:
    processed = np.all(np.isfinite(YY), axis=0)
    irec_start = np.count_nonzero(processed)
    print(f"Next record to start {irec_start}")

YY, JG, IG = construct_ithkn(DNMB, DIRS, JG, IG, dfltmp, irec_start, YY, dump_tstp=10)
  

# Save:
if irec_start < len(DNMB):
  print(f"Final Saving ithkn time series and IG, JG --> {dfltmp}")
  np.savez(dfltmp,
         YY=YY,
         JG=JG,
         IG=IG,
         DNMB=DNMB)



check_recs = False
if check_recs:
  nrecs = len(DNMB)
  years = np.zeros((nrecs))
  mms   = np.zeros((nrecs))
  dds   = np.zeros((nrecs))

  for irec, d0 in enumerate(DNMB):
    yr, mm, dd = mtime.datevec(d0)[:3]
    years[irec] = yr
    mms[irec] = mm
    dds[irec] = dd 

  Nyrs = []
  for yr0 in range(int(years[0]), int(years[-1])+1):
    iy = len(np.where(years==yr0)[0])
    Nyrs.append(iy)

  Nyrs = np.asarray(Nyrs, dtype=int)



f_check = False
if f_check:
  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
  ax1.contour(LMsk, [0.99], linestyles='solid', colors=[(0.5,0.8,1)])
  ax1.scatter(IG, JG, s=5, color=(0.8,0.4,0), marker='.')


