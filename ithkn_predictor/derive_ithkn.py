"""
  Derive response variable: ice thickness
  extract all time steps at specified grid points
  Update grid points: eliminate 0-thickness grid points

  Save tmp file to be used by other scripts
  deriving predictor variables

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
parser.add_argument("--ndays", help="Use N-th day, 0, ndays, 2*ndays, ..., =0 - derive days from T2m ERA5", 
                    default=0, type=int)
parser.add_argument("--dxy", help=f"Min dist (km) between data points (~corr.scale), to skip close i,j points", 
                    type=int, required=True)
parser.add_argument("--ys", help="Year start, default=1993", default=1993, type=int)
parser.add_argument("--ye", help="Year end, default=2025", default=2025, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'], 
                    required=True, type=str)
parser.add_argument("--load", help="Load saved ithkntmp, continue from last record (1), start from time 0 (0)", 
                  choices=[0,1], required=True, type=int)
args = parser.parse_args()

ndays = args.ndays
dxy   = args.dxy    
YS    = args.ys
YE    = args.ye
regn  = args.regn
load_saved = args.load == 1

ndays_era = 7   # freq. of saved era5 fields
ithkn_max = 4.  # cap max ice thickness
Ntime_avrg = 3 # for unrealistic (<0) values, use average values from previous records

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
  "ptht2m"   : config_predictor["linregr"]["ptht2m"].format(regn_name=regn_name),
  "pthout"   : config_predictor["linregr"]["pthout"],
  "ithkntmp" : config_predictor["linregr"]["ithkntmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  }

def derive_time(YS, YE, ndays, DIRS, regn_name):
  DNMB = None
  time_stmp = []
  if ndays > 0:
    print(f"Deriving time array for {ndays} skip days")
    for YR in range(YS,YE+1):
      dnmb0 = mtime.datenum([YR,1,1])
      dnmbE = mtime.datenum([YR,12,31])
      DYR = np.arange(dnmb0, dnmbE+1, ndays)
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

# check DNMB array:
check_dnmb = True
if check_dnmb:
  micepr.check_dnmb_array(DNMB, print_months=False)

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

def update_icepnts(YY, JG, IG, hice_min=0.3, nprc=0.1):
  """
    Eliminate points that have ice > hice_min
    less than nprc of time
  """
  #indx_zeros = np.all(YY < hice_min, axis=1)
  #ikeep = ~indx_zeros
  ikeep = np.count_nonzero(YY >= hice_min, axis=1) >= nprc * YY.shape[1]
  nkeep = ikeep.sum()  # how many points retained
  print(f"Original count {YY.shape[0]}, keep points = {nkeep}")
  YY = YY[ikeep, :]
  JG = JG[ikeep]
  IG = IG[ikeep]

  return YY, JG, IG


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
  
# Eliminate unneeded grid points with no ice
YY, JG, IG = update_icepnts(YY, JG, IG)

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


