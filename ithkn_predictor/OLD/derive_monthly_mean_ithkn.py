"""
  Derive response variable: ice state
  that captures the long-term interannual trend in
  ice thicknesses - ice getting thineer in the Arctic in 2020s 
  compared to 1990s

  Possible predictors: last-year monthly ice volume in the Arctic or 
  area-mean ice thickness

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
from mod_mom6 import dx_dy

#from MyPython.mod_cice6_utils import change_base_template, flname_replace_date

parser = argparse.ArgumentParser()
parser.add_argument("--ys", help="Year start, default=1993", default=1993, type=int)
parser.add_argument("--ye", help="Year end, default=2025", default=2025, type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'], 
                    required=True, type=str)
args = parser.parse_args()

YS    = args.ys
YE    = args.ye
regn  = args.regn

dlt_days = 7   # days to skip
MDAYS = np.array([1,7,14,21,28])

# Only central Arctic 
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
  }

def derive_time_adhock(YS, YE, MDAYS, DIRS, regn_name):
  DNMB = []
  print(f"Deriving time array for {dlt_days} skip days")
  for YR in range(YS,YE+1):
    for MM in range(1,13):
      for DD in MDAYS:
        dnmb0 = mtime.datenum([YR,MM,DD])
        DNMB.append(dnmb0)

  return np.asarray(DNMB)

ptht2m = DIRS["ptht2m"]
DNMB = micepr.derive_time(YS, YE, ptht2m, regn_name, dlt_days)

# Read GLORYS grid:
pthice = os.path.join(DIRS["pthithkn"],f"{YS}")
#rdate = f"{YS*10000+100+1}"
rdate = YS*10000+100+1

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


nrec = len(DNMB)
nmonths = 12 * (YE-YS+1)
IVOL = np.zeros(nmonths)
ITHKM = np.zeros(nmonths)*np.nan
DNMBR = np.zeros(nmonths)*np.nan
mold = -1
imm = -1
iday = 0
for dnmb0 in DNMB:
  YR, MM, DD = mtime.datevec(dnmb0)[:3]
 
  if mold < 0:
    mold = MM
    yrold = YR
    ivol_mean = 0.
    ithkn_mean = 0.
  elif mold != MM:
    imm += 1
    # Monthly means:
    # Update records at the end of the month
    if iday > 0:
      ivol_mean /= iday
      ithkn_mean /= iday
      print(f"{yrold}/{mold:02d}: ivol_mean = {ivol_mean:.2f}km3, mean ithkn = {ithkn_mean:.2f}m")
      IVOL[imm] = ivol_mean
      ITHKM[imm] = ithkn_mean
    DNMBR[imm] = mtime.datenum([yrold, mold, 15])

    iday = 0
    ivol_mean  = 0.
    ithkn_mean = 0.
    mold = MM
    yrold = YR

  #print(f"Reading ithkn {YR}/{MM:02d}/{DD:02d}")
  pthice = os.path.join(DIRS["pthithkn"],f"{YR}")
  rdate = int(YR*1e4 + MM*100 + DD)
  dflice = mglr.find_file(rdate, pthice)
  if dflice is None:
    #raise FileNotFoundError(f"No ice thickness file for {rdate}")
    print(f"No ice thickness file for {rdate}")
    print("Skipping ...")
    continue

  with xr.open_dataset(dflice) as dsice:
    A2d = dsice['sithick'].isel(time=0).values.squeeze()
 
  # Treat nans as no ice grid cells and mask out outside domain
  ITHKN = np.nan_to_num(A2d, nan=0.0)
  ITHKN[~DOMAIN] = np.nan

  # Ice conc:
  pthiconc = os.path.join(DIRS["pthiconc"],f"{YR}")
  rdate = int(YR*1e4 + MM*100 + DD)
  dfliconc = mglr.find_file(rdate, pthiconc)
  if dfliconc is None:
    raise FileNotFoundError(f"No ice thickness file for {rdate}")

  with xr.open_dataset(dfliconc) as dsice:
    A2d = dsice['siconc'].isel(time=0).values.squeeze()

  # Treat nans as no ice grid cells and mask out outside domain
  ICONC = np.nan_to_num(A2d, nan=0.0)
  ICONC[~DOMAIN] = np.nan
 
  ivol = np.nansum(ITHKN*ICONC*Acell)  # m3
  iarea = np.nansum(Acell * ICONC)
  ivol_mean += ivol * 1e-9             # km3

  # Mean ice thickness - over ice area!
  ithkn_day = ivol / iarea
  ithkn_mean += ithkn_day

  iday += 1

# Save last month:
if iday > 0:
  ivol_mean /= iday
  ithkn_mean /= iday
  imm += 1
  IVOL[imm] = ivol_mean
  ITHKM[imm] = ithkn_mean
  DNMBR[imm] = mtime.datenum([yrold, mold, 15])

pthout = DIRS["pthout"]
fltmp = f"GLORYS_monthly_icevol_ithknmn_{regn}_{YS}_{YE}.npz"
dflout = os.path.join(pthout, fltmp)
print(f"Saving ice vol and mean ice thickness --> {dflout}")
np.savez(dflout,
       IVOL=IVOL,
       ITHKM=ITHKM,
       DNMB=DNMBR)




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


