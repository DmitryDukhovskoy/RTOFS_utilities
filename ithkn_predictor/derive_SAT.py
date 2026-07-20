"""
  Derive dynamic predictor: ERA5 surface air temperature

  Need mapping indices gmapi to map ERA5 --> GLORYS grid
  find_remap_indx_era5_to_GLORYS.py

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
from mod_misc1 import dist_sphcrd
from mod_mom6 import dx_dy

#from MyPython.mod_cice6_utils import change_base_template, flname_replace_date

parser = argparse.ArgumentParser()
parser.add_argument("--dxy", help=f"Min dist (km) between data points (~corr.scale), to skip close i,j points", 
                    type=int, required=True)
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

dump_tstp = 10
fld_name = 'SAT'

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
  "pthgmapi" : config_predictor["linregr"]["pthgmapi"],
  "pthout"   : config_predictor["linregr"]["pthout"],
  "ithkntmp" : config_predictor["linregr"]["ithkntmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "iconctmp" : config_predictor["linregr"]["iconctmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "ssttmp"   : config_predictor["linregr"]["iconctmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "divutmp"  : config_predictor["linregr"]["divutmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "sattmp"   : config_predictor["linregr"]["sattmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  }


def derive_time(YS, YE, DIRS, regn_name):
  DNMB = None
  time_stmp = []
  ndays_era = 7   # freq. of saved era5 fields

  print(f"Deriving time array from ERA5 fields")
  # Derive Time from saved atm. fields:
  ptht2m = DIRS['ptht2m']
  for YR in range(YS,YE+1):
    flnm = f"era5_2mTemp_daily{ndays_era}day_{regn_name}_{YR}.nc"
    dflnm = os.path.join(ptht2m, flnm)
    assert os.path.isfile(dflnm), f"Missing ERA5: {dflnm}, check ndays flag"

    dnmb0 = mtime.datenum([YR,1,1])
    with xr.open_dataset(dflnm, decode_times=False) as ds:
      Time = ds["valid_time"].values
      DYR = dnmb0 + Time
    time_stmp.append(DYR)

  DNMB = np.concatenate(time_stmp)
  return DNMB

DNMB = derive_time(YS, YE, DIRS, regn_name)


# Read time array ad J,I sample grid points:
# Time array should match ERA5 extracted fields
pthout = DIRS["pthout"]
flithkn = DIRS["ithkntmp"]
dflithkn = os.path.join(pthout, flithkn)
  
print(f"Loading saved {dflithkn}, will start from last saved record")
if not os.path.isfile(dflithkn):
  print(f"Missing tmp file {dflithkn}\n  start from time = 0")
else:
  data = np.load(dflithkn)
  JG = data["JG"]
  IG = data["IG"]
  DNMB_check = data["DNMB"]
  nrecs = len(DNMB)

  # Check that this is the right time series:
  dtmp = np.floor(np.abs(DNMB - DNMB_check))
  assert np.max(dtmp) == 0, "Check DNMB - dates do not match with saved time series"


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
#DX, DY = dx_dy(hlon, hlat)
#Acell = DX*DY

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

# Construct predictor sst time series for all locations, 
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

# Find GLORYS - ERA5 pairs:
print("Finding JERA, IERA to match JGLR, IGLR")
JE = np.empty(len(JG), dtype=int)
IE = np.empty(len(IG), dtype=int)
for k, (jj, ii) in enumerate(zip(JG, IG)):
  if k > 0 and k % 500 == 0:
    prc = k/len(JG)*100.
    print(f"  {prc:.2f}% processed")
  JE[k], IE[k] = find_era_indx(jj, ii, JERA, IERA, JGLR, IGLR)

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

check_gmapi = True
if check_gmapi:
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

  flt2m = f"era5_2mTemp_daily7day_Arctic_{YR}.nc"
  ptht2m = DIRS['ptht2m']
  dflt2m = os.path.join(ptht2m, flt2m)
  if YR != YRold:
    YRold = YR
    with xr.open_dataset(dflt2m, decode_times=False) as ds:
      Time = ds['valid_time'].values
    dnmbS = mtime.datenum([YR,1,1])
    TM = (Time + dnmbS).astype(int)

  #idx = np.where(np.isclose(TM, dnmb0))[0]
  idx = np.where(TM == int(dnmb0))[0]
  if len(idx) == 0:
      raise ValueError(f"No matching day for {dnmb0} {YR}/{MM}/{DD}")
  iday = idx[0]

  with xr.open_dataset(dflt2m) as dsice:
    A2d = dsice['t2m'].isel(valid_time=iday).values.squeeze()

  T2d = A2d - 273.15  # K --> C

  fld_pnts = []
  for jje, iie in zip(JE, IE):
    fld_pnts.append(T2d[jje,iie])

  YY[:,irec] = np.asarray(fld_pnts)

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

