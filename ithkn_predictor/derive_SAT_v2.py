"""
  version 2

  Derive dynamic predictor: ERA5 surface air temperature (T2m)
  Antarctic or Arctic 

  Use 1hr ERA5 fields on PPAN / uda 
  Need to derive daily then map on GLORYS analysis indices J,I

  Need mapping indices gmapi to map ERA5 --> GLORYS grid
  find_remap_indx_era5_to_GLORYS.py

  Use time stamps and J,I grid points defined in define_time_IJpnts.py 

  Data:
/uda/ERA5/Hourly_Data_On_Single_Levels/reanalysis/global/1hr/annual_file-range/Temperature_and_Pressure/2m-temperature/0.25x0.25

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
from mod_misc1 import dist_sphcrd
from mod_mom6 import dx_dy

#from MyPython.mod_cice6_utils import change_base_template, flname_replace_date

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="Region to process", choices=['north','south'], 
                    default='south')
parser.add_argument(
    "--load", 
    help="Load saved sst tmp file, continue from last record (1), start from time 0 (0)", 
    choices=[0,1], 
    required=True, 
    type=int
)
args = parser.parse_args()

regn  = args.regn
load_saved = args.load == 1

dump_tstp = 10
fld_name = 'SAT'


fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

regn_name = config_predictor["regn"][regn]["name"]
lat0      = config_predictor["regn"][regn]["lat_bnd"]
tstep     = config_predictor["params"]["tstep"]
dxy       = config_predictor["params"]["dxy"]
YS        = config_predictor["params"]["ys"]
YS        = config_predictor["params"]["ys"]
YE        = config_predictor["params"]["ye"]

DIRS = {
  "pthithkn" : config_predictor["linregr"]["pthithkn"],
  "pthiconc" : config_predictor["linregr"]["pthiconc"],
  "pthsst"   : config_predictor["linregr"]["pthsst"],
  "pthssh"   : config_predictor["linregr"]["pthssh"],
  "pthui"    : config_predictor["linregr"]["pthui"],
  "pthvi"    : config_predictor["linregr"]["pthvi"],
  "ptht2m"   : config_predictor["linregr"]["ptht2m_1hr"],
  "pthgmapi" : config_predictor["linregr"]["pthgmapi"],
  "pthout"   : config_predictor["linregr"]["pthout"],
  "ithkntmp" : config_predictor["linregr"]["ithkntmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "iconctmp" : config_predictor["linregr"]["iconctmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "ssttmp"   : config_predictor["linregr"]["ssttmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "divutmp"  : config_predictor["linregr"]["divutmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "sattmp"   : config_predictor["linregr"]["sattmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
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

print(f"N of grid points: {npnts}, N time steps: {nrecs}")

# GLORYS grid:
A = np.load(dflgrid)
hlon = A["LON"]
hlat = A["LAT"]
LMsk = A["LMsk"]

# 0 <= lon < 360
hlon = (hlon + 360) % 360


# Load gmapi:
pthindx = DIRS["pthgmapi"]
flout = f"gmapi_ERA5_to_GLORYS_{regn}.nc"
dflout = os.path.join(pthindx, flout)
with xr.open_dataset(dflout) as ds:
  LONE = ds["era_longit"].values
  LATE = ds["era_latit"].values
  IGLR = ds["glorys_indx"].values
  JGLR = ds["glorys_jndx"].values
  IERA = ds["era_indx"].values
  JERA = ds["era_jndx"].values

LONE = (LONE + 360) % 360

# Construct predictor T2m time series for all locations, 
# Or load previously saved
pthout = DIRS["pthout"]
fltmp = DIRS["sattmp"]
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


def find_era_indx(jj, ii, JERA, IERA, JGLR, IGLR):
  """
    Given glorys grid pnt (ii,jj) 
    Find corresponding ERA grd pnt
    using gmapi 
  """
  DD = (JGLR-jj)**2 + (IGLR-ii)**2
  idx = np.argmin(DD)
  assert np.floor(DD[idx]) == 0, f"Could not match GLORYS index jj={jj} ii={ii}"
  
  return JERA[idx], IERA[idx]

def check_era_glorys_coord(LATE, LONE, IE, JE, hlon, hlat, IG, JG, dmax=27e3):
  """
   Check that gmapi indices are correct
   dmax ~ ERA5 res. (0.25 degree)
  """
  print("Checking gmapi indices ...")
  DD_max = 0
  for k, (jj, ii) in enumerate(zip(JG, IG)):
    jje = JE[k]
    iie = IE[k]
    lon_era = LONE[iie]
    lat_era = LATE[jje]
    lon_glr = hlon[jj,ii]
    lat_glr = hlat[jj,ii]

    DD = dist_sphcrd(lat_era, lon_era, lat_glr, lon_glr)
    if k > 0 and k % 500 == 0:
      prc = k/len(JG)*100.
      print(f"  {prc:.2f}% processed, k={k} dist = {DD:.4f} m")

    DD_max = np.max([DD_max, DD])

    if DD > dmax:
      print(f"ERA point is {DD}m apart from GLORYS, k={k}")
      print(f"ERA lon={lon_era} lat={lat_era}")
      print(f"GLORYS lon={lon_glr} lat={lat_glr}")
      raise Exception("ERR: check not passed")

  print(f"Checked gmapi: OK, overall max dist = {DD_max} m")

# Find GLORYS - ERA5 pairs if not saved:
flg2e = f"glorys2era_pairs_{regn}.npz"
dflg2e = os.path.join(pthinfo, flg2e)
if os.path.isfile(dflg2e):
  print(f"Reading ERA5 indices corresponding GLORYS grid points")
  A = np.load(dflg2e)
  JE = A["JE"]
  IE = A["IE"]
  
else:
  print("Finding JERA, IERA to match JGLR, IGLR")
  JE = np.empty(len(JG), dtype=int)
  IE = np.empty(len(IG), dtype=int)

  for k, (jj, ii) in enumerate(zip(JG, IG)):
    if k > 0 and k % 500 == 0:
      prc = k/len(JG)*100.
      print(f"  {prc:.2f}% processed")
    JE[k], IE[k] = find_era_indx(jj, ii, JERA, IERA, JGLR, IGLR)

  print(f"Saving GLORYS --> ERA indices: {dflg2e}")
  np.savez(dflg2e, JE=JE, IE=IE)  

check_gmapi = True
if check_gmapi:
  # If failed - check IE, JE, may need to redefine those again
  check_era_glorys_coord(LATE, LONE, IE, JE, hlon, hlat, IG, JG)

"""
  Construct time series of response variable SAT
  2D: locations x time
"""
npnts = len(JG)
nrecs = len(DNMB)
YRold = 1900
if YY is None:
  YY = np.zeros((npnts, nrecs), dtype=float)*np.nan
for irec, dnmb0 in enumerate(DNMB):
  YR, MM, DD = mtime.datevec(dnmb0)[:3]

  if irec < irec_start:
    print(f"Skipping, already processed: {YR}/{MM:02d}/{DD:02d}")
    continue
  
  print(f"Reading {fld_name} {YR}/{MM:02d}/{DD:02d}")

  rdate = int(YR*1e4 + MM*100 + DD)

  #flt2m = f"era5_2mTemp_daily7day_Arctic_{YR}.nc"
  flt2m = f"ERA5_reanalysis_sLevels_1hr_0.25x0.25_2m-temperature_{YR}.nc"
  ptht2m = DIRS['ptht2m']
  dflt2m = os.path.join(ptht2m, flt2m)

  # Read time coord for a new year
  if YR != YRold:
    YRold = YR
    with xr.open_dataset(dflt2m, decode_times=False) as ds:
      Time_hrs = ds['time'].values
    dnmbS = mtime.datenum([1900,1,1])
    TM = (dnmbS + Time_hrs / 24.0)

  # hourly ---> Daily
  # Find hourly records in this day:
  idx = np.flatnonzero((TM >= dnmb0) & (TM < dnmb0 + 1)) # indices where condition is True
  if len(idx) == 0:
    print(f"WARNING: No hourly data found for {YR}/{MM:02d}/{DD:02d}")
    continue

  idx1 = idx[0]
  idx2 = idx[-1]

  # Average
  icc = 0
  T2d = None
  with xr.open_dataset(dflt2m) as dsice:
    for ill in range(idx1, idx2+1): 
      A2d = dsice['t2m'].isel(time=ill).values.squeeze()

      if T2d is None:
        T2d = np.zeros_like(A2d, dtype=float)

      icc += 1
      T2d += A2d - 273.15  # K --> C

  T2d /= icc
  YY[:,irec] = T2d[JE, IE]

  # Temporary save:
  if (irec + 1) % dump_tstp == 0: 
    print(f"TMP step: Saving {fld_name} time series  --> {dfltmp}")
    np.savez(dfltmp,
           YY=YY,
           JG=JG,
           IG=IG,
           DNMB=DNMB)

# Final save
if irec_start < len(DNMB):
  # No need to save if already everything processed
  print(f"END TMP step: Saving {fld_name} time series and IG, JG --> {dfltmp}")
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
  #ax1.scatter(IG, JG, s=5, color=(0.8,0.4,0), marker='.')

  # Plot ERA5 SAT at the sample locations:
  ir0 = 0
  AA = YY[:,ir0]
  sc = ax1.scatter(
    IG, JG,
    c=AA,
    cmap='jet',
    s=20,          # marker size
    vmin=-30,        # optional color scale limits
    vmax=0
 )

  plt.colorbar(sc, ax=ax1, label='SAT, degC')

