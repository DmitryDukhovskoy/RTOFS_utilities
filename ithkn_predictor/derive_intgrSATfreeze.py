"""
  Derive dynamic predictor: sqrt of the number of freeze degree days
  Following Zubov's relation: h2 + 50h = 8 IFDD, 
  IFDD = sum of (Tfrz - Tair), when Tair < Tfrz

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
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d

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

intgr_time = 90  # Time for freeze degree days accumulation, back from current time
dump_tstp = 20
fld_name = 'intgrSATfreeze'
Tfrz = -1.85    # ocea freezing T
ndays_era = 7   # freq. of saved era5 fields

regions = {
    "north": ("Arctic", 65.0),
    "south": ("Antarctic", -60.0),
}
regn_name, lat0 = regions[regn]

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
  "ssttmp"   : config_predictor["linregr"]["ssttmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "divutmp"  : config_predictor["linregr"]["divutmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "sattmp"   : config_predictor["linregr"]["sattmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "dfrztmp"  : config_predictor["linregr"]["dfrztmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  }


def derive_time(YS, YE, DIRS, regn_name, ndays_era):
  DNMB = None
  time_stmp = []

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

DNMB = derive_time(YS, YE, DIRS, regn_name, ndays_era)
dnmbStart = DNMB[0]
# Add previous days for integrating SAT
# integrating freezing deegre days
# Previous (to start) year should exist !
dnmbS = DNMB[0]   # actual start day
dnmbP = dnmbS - intgr_time - 1  # previous intgr time preiod, start day
Ypr, Mpr, Dpr = mtime.datevec(dnmbP)[:3]
DNMBprv = derive_time(Ypr, Ypr, DIRS, regn_name, ndays_era)

# Find closest time:
idx0 = max(np.argmin(abs(DNMBprv - dnmbP)) - 2, 0) # add extra index
nrec_prev = len(DNMBprv) - idx0     # how many records to keep for intgr Tfrz
# prepand first days for integrating Tfrz before the start:
DNMB_run = DNMB.copy()
DNMB = np.concatenate((DNMBprv[idx0:], DNMB))

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

  # Check that this is the right time series:
  dtmp = np.floor(np.abs(DNMB_run - DNMB_check))
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
fltmp = DIRS["dfrztmp"]
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
    DNMBprv = data["DNMBprv"]
    DNMB_check = data["DNMB"]  

    # Check that this is the right time series:
    dtmp = np.floor(np.abs(DNMB_run - DNMB_check))
    assert np.max(dtmp) == 0, "Check DNMB - dates do not match with saved time series"

    DNMB = np.concatenate((DNMBprv, DNMB_check))
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
assert len(JG)==len(IG)
JE = np.empty(len(JG), dtype=int)
IE = np.empty(len(IG), dtype=int)
# Alternative to distance-approach, build a lookup table:
gmapi = {(jg,ig):(je,ie)
         for jg,ig,je,ie in zip(JGLR,IGLR,JERA,IERA)}

for k, (jj, ii) in enumerate(zip(JG, IG)):
  if k > 0 and k % 500 == 0:
    prc = k/len(JG)*100.
    print(f"  {prc:.2f}% processed")
  #JE[k], IE[k] = find_era_indx(jj, ii, JERA, IERA, JGLR, IGLR)
  key = (int(jj), int(ii))
  if key not in gmapi:
    raise ValueError(f"Missing gmapi entry for GLORYS index {key}")
  JE[k], IE[k] = gmapi[key]


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

def intgr_Tfrz(SATprv, intgr_time, Tfrz, days_frz):
  """
    Use Zubov definition of accumulated freeze degree days
    sum(Tfrz-T), when SAT < Tfrz
    Linear interpolation is safer
  """
  Tintrp = np.arange(days_frz[-1] - intgr_time, days_frz[-1]+1)
  #cs = CubicSpline(days_frz, SATprv, axis=1)
  # SATs during the requested previous Ndays
  #SATi = cs(Tintrp)

  interp = interp1d(days_frz, SATprv,
                  axis=1,
                  kind='linear')
  SATi = interp(Tintrp)

  T2frz = np.where(SATi < Tfrz, SATi, np.nan) 
  intgrFDD = np.nansum(Tfrz - T2frz, axis=1)  # integrated Freeze degree days
 
  f_chck = False
  if f_chck:
    iC=15
    Ti = SATi[iC,:]
    T0 = SATprv[iC,:]
    ax1.cla()
    ax1.plot(Tintrp, Ti)
    ax1.plot(days_frz, T0,'.-')

  return intgrFDD 


check_gmapi = True
if check_gmapi:
  check_era_glorys_coord(LATE, LONE, IE, JE, hlon, hlat, IG, JG)

SATprv = np.zeros((len(IE), nrec_prev)) # should include intgr. period
days_frz = np.zeros(nrec_prev)    # time stamps of saved SAT
"""
  Construct time series of response variable SAT
  2D: locations x time

  First N records will be skipped before actual start date
  To populate SAT array with temp for integrating Freez. days
"""
iStart = np.where(DNMB == dnmbStart)[0][0]
npnts = len(JG)
nrecs = len(DNMB_run)
irec = 0
YRold = 1900
if YY is None:
  YY = np.zeros((npnts, nrecs), dtype=float)*np.nan
for irec0, dnmb0 in enumerate(DNMB):
  YR, MM, DD = mtime.datevec(dnmb0)[:3]

  print(f"Reading ERA5 SAT {YR}/{MM:02d}/{DD:02d}")

  rdate = int(YR*1e4 + MM*100 + DD)

  #flt2m = f"era5_2mTemp_daily7day_Arctic_{YR}.nc"
  flt2m = f"era5_2mTemp_daily{ndays_era}day_{regn_name}_{YR}.nc"
  ptht2m = DIRS['ptht2m']
  dflt2m = os.path.join(ptht2m, flt2m)
  if YR != YRold:
    if YRold != 1900:
        ds_t2m.close()
    YRold = YR
    ds_t2m = xr.open_dataset(dflt2m, decode_times=False)
    Time = ds_t2m['valid_time'].values
    dnmb_day1 = mtime.datenum([YR,1,1])
    TM = (Time + dnmb_day1).astype(int)

  #idx = np.where(np.isclose(TM, dnmb0))[0]
  idx = np.where(TM == int(dnmb0))[0]
  if len(idx) == 0:
    raise ValueError(f"No matching day for {dnmb0} {YR}/{MM}/{DD}")
  iday = idx[0]
  A2d = ds_t2m['t2m'].isel(valid_time=iday).values.squeeze()
  T2d = A2d - 273.15  # K --> C

  fld_pnts = []
  Tsurf = []
  for jje, iie in zip(JE, IE):
    Tsurf.append(T2d[jje,iie])

  # Update SAT with previous SATs:
  # add current day at the end and delete 1st column 
  SATprv[:, :-1] = SATprv[:, 1:]
  SATprv[:,-1] = Tsurf
  days_frz[:-1] = days_frz[1:]
  days_frz[-1] = dnmb0

  if dnmb0 < dnmbS:
    continue

  irec = irec0 - iStart
  assert irec >= 0
  assert irec < YY.shape[1]
  if irec < irec_start:
    continue

  assert np.all(days_frz > 0), "days_frz not populated, there are 0s"
  assert np.all(np.diff(days_frz)>0), "days_frz not increasing" 
  fld_pnts = intgr_Tfrz(SATprv, intgr_time, Tfrz, days_frz)

  YY[:,irec] = np.asarray(fld_pnts)

  if (irec + 1) % dump_tstp == 0: 
    print(f"TMP step: Saving {fld_name} time series  --> {dfltmp}")
    np.savez(dfltmp,
           YY=YY,
           JG=JG,
           IG=IG,
           JE=JE,
           IE=IE,
           DNMBprv=DNMB[:iStart],
           DNMB=DNMB[iStart:])

if irec_start < len(DNMB):
  # No need to save if already everything processed
  print(f"END TMP step: Saving {fld_name} time series and IG, JG --> {dfltmp}")
  np.savez(dfltmp,
         YY=YY,
         JG=JG,
         IG=IG,
         JE=JE,
         IE=IE,
         DNMBprv=DNMB[:iStart],
         DNMB=DNMB[iStart:])

ds_t2m.close()


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

