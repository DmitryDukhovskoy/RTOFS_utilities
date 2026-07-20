"""
  As a first step in developing ice thickness predictor
  train a linear reg. model to test different predictors
  and check if the taks is feasable, e.g. if it can be treated
  as linear problem

  Use every N-days of daily data
  and subsample high-res. GLORYS fields: there is lots of spatial
  correlation, no need to do every grid point

The daily data is on uda:
/uda/Global_Ocean_Physics_Reanalysis/global/daily/siconc/

/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick

ERA5 T2m fields:
  Downloaded every 7-day daily mean 2m SAT for specified region from ERA5 website
  https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics?tab=download

OR use 1-hr fields on PPAN --> derive daily mean
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

#from MyPython.mod_cice6_utils import change_base_template, flname_replace_date

parser = argparse.ArgumentParser()
parser.add_argument("--ndays", help="Use N-th day, 0, ndays, 2*ndays, ..., =0 - derive days from T2m data", 
                    default=0, type=int)
parser.add_argument("--dxy", help=f"Min dist (km) between data points (~corr.scale), to skip close i,j points", 
                    type=int, required=True)
parser.add_argument("--ys", help="Year start, default=1993", default=1993, type=int)
parser.add_argument("--ye", help="Year end, default=2025", default=2025, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'], 
                    required=True, type=str)
args = parser.parse_args()

ndays = args.ndays
dxy   = args.dxy    
YS    = args.ys
YE    = args.ye
regn  = args.regn

ndays_era = 7   # freq. of saved era5 fields

if regn == 'north':
  regn_name = 'Arctic'
  lat0 = 65.
elif regn == 'south':
  regn_name = 'Antarctic'
  lat0 = -60.

DIRS = {
  "pthithkn" : "/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick/",
  "pthiconc" : "/uda/Global_Ocean_Physics_Reanalysis/global/daily/siconc/",
  "pthsst"   : "/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick/",
  "pthssh"   : "/uda/Global_Ocean_Physics_Reanalysis/global/daily/zos/",
  "ptht2m"   : f"/archive/Dmitry.Dukhovskoy/data/ERA5/{regn_name}",
  "pthout"   : "/work/Dmitry.Dukhovskoy/anls_output/GLORYS_anls/ice_linregr",
  }

def derive_time(YS, YE, ndays, DIRS, regn_name):
  DNMB = None
  time_stmp = []
  if ndays > 0:
    print(f"Deriving time array for {ndays} skip days")
    for YR in range(YS,YE+1):
      dnmb0 = mtime.datenum([YR,1,1])
      dnmbE = mtime.datenum([YR,12,31])
      DYR = np.arange(dnmb0, dnmbEi+1, ndays)
      time_stmp.append(DYR)

    DNMB = np.concatenate(time_stmp)  
    return DNMB

  else:
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

DNMB = derive_time(YS, YE, ndays, DIRS, regn_name)

# Read GLORYS grid:
pthice = os.path.join(DIRS["pthithkn"],f"{YS}")
rdate = f"{YS*10000+100+1}"

# Find file:
dflglr = mglr.find_file(rdate, pthice)

with xr.open_dataset(dflglr) as dsice:
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values

# 0 <= lon < 360
LON = (LON + 360) % 360
hlon, hlat = np.meshgrid(LON, LAT)

# Land mask:
LMsk = None
pthssh = os.path.join(DIRS["pthssh"],f"{YS}")
dflssh = mglr.find_file(rdate, pthssh)
with xr.open_dataset(dflssh) as dszos:
  SSH = dszos['zos'].isel(time=0).values.squeeze()

LMsk = np.where(np.isfinite(SSH),1,0)


# Define domain:
if regn == 'north':
  DOMAIN = (hlat > lat0) & (LMsk == 1)
elif regn == 'south':
  DOMAIN = (hlat < lat0) & (LMsk == 1)

# Exclude N. Atlantic and Barents Sea:
DOMAIN[:1893,2154:2795] = False
DOMAIN[:1764,:258] = False   # N. Bering 
DOMAIN[:1807,1991:2243] = False # Iceland Sea

# Subset points based on minimum distance criterion:
# Estimate dltJ, assuming 1dgr ~= 111 km:
Rearth = 6371
dy1dgr = Rearth*np.pi/180

dlty_dgr = np.diff(LAT)[0]  # regular grid
dlty_km = dlty_dgr * dy1dgr
dltJ = int(dxy / dlty_km)
if dltJ == 0:
  dltJ = 1 

dltx_dgr = np.diff(LON)[0]


JG = []
IG = []
nj, ni = DOMAIN.shape
jj0 = np.argmax(DOMAIN.any(axis=1)) # gives the index of the 1st row containing any valid grid pnt
for jj in range(jj0, nj, dltJ):
  phi = LAT[jj]
  dx1dgr = np.cos(np.deg2rad(phi)) * Rearth * np.pi / 180
  dltx_km = dltx_dgr * dx1dgr
  dltI = max(1, int(np.ceil(dxy / dltx_km))) 
 
  iold = -np.inf
  for ii in np.flatnonzero(DOMAIN[jj]):
    if ii - iold >= dltI:
      IG.append(ii)
      JG.append(jj)
      iold = ii

Npnts = len(IG)
print(f"For dxy={dxy} km and lat0={lat0:.2f}, Selected N pnts={Npnts}")

IG = np.asarray(IG, dtype=int)
JG = np.asarray(JG, dtype=int)

def update_icepnts(YY, JG, IG):
  """
    Eliminate points that never have ice
  """
  indx_zeros = np.all(YY == 0, axis=1)
  ikeep = ~indx_zeros
  YY = YY[ikeep, :]
  JG = JG[ikeep]
  IG = IG[ikeep]

  return YY, JG, IG

def construct_ithkn(DNMB, DIRS, JG, IG):
  """
    Construct time series of response variable (ithkn)
    2D: locations x time
  """
  npnts = len(JG)
  nrecs = len(DNMB)
  YY = np.empty((npnts, nrecs), dtype=float)
  for irec, dnmb0 in enumerate(DNMB):
    YR, MM, DD = mtime.datevec(dnmb0)[:3]
    print(f"Reading ithkn {YR}/{MM:02d}/{DD:02d}")

    pthice = os.path.join(DIRS["pthithkn"],f"{YR}")
    rdate = int(YR*1e4 + MM*100 + DD)
    dflice = mglr.find_file(rdate, pthice)

    with xr.open_dataset(dflice) as dsice:
      A2d = dsice['sithick'].isel(time=0).values.squeeze()
 
    # 0 values are nans?
    A2d = np.nan_to_num(A2d, nan=0.0)
 
    fld_pnts = A2d[JG,IG]
    YY[:,irec] = fld_pnts
  
    if irec%tmp_tstp == 0: 
      print(f"TMP step: Saving ithkn time series and IG, JG --> {dfltmp}")
      np.savez(dfltmp,
             YY=YY,
             JG=JG,
             IG=IG,
             DNMB=DNMB)

      YY, JG, IG = update_icepnts(YY, JG, IG)

  YY, JG, IG = update_icepnts(YY, JG, IG)

  return YY, IG, JG

# If needed: Construct response (ithkn) time series concatenating all locations, 
# locations with 0 ithkn will be eliminated from IG, JG
# Or load previously saved
pthout = DIRS["pthout"]
fltmp = f"ithkn_IJpnts_tser_{YS}-{YE}_{regn}.npz"
dfltmp = os.path.join(pthout, fltmp)
if derive_ithkn:
  YY, IG, JG = construct_ithkn(DNMB, dfltmp, tmp_tstp=10)
  
  # Save:
  print(f"Saving ithkn time series and IG, JG --> {dfltmp}")
  np.savez(dfltmp,
         YY=YY,
         JG=JG,
         IG=IG,
         DNMB=DNMB)

else:
  data = np.load(dfltmp)
  YY = data["YY"]
  JG = data["JG"]
  IG = data["IG"]
  DNMB_check = data["DNMB"]  

  # Check that this is the right time series:
  dtmp = np.floor(np.abs(DNMB - DNMB_check))
  assert np.max(dtmp) == 0, "Check DNMB - dates do nnot match with saved time series"



plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
ax1.contour(LMsk, [0.99], linestyles='solid', colors=[(0.5,0.8,1)])
ax1.scatter(IG, JG, s=5, color=(0.8,0.4,0), marker='.')


