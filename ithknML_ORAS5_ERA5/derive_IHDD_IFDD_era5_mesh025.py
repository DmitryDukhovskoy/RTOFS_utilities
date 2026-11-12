"""
  Derive dynamic predictor: the number of freeze degree days
  Following Zubov's relation: h2 + 50h = 8 IFDD, 
  IFDD = sum of (Tfrz - Tair), when Tair < Tfrz

  Heat degree days
  in analogy to Zubov's IFDD, but opposite
  which may be important for melt season

  Check whether sqrt(IFDD) is used as a predictor
  Use ERA5 T2m hourly fields
  Interpolated IFDD onto mesh025

  Need mapping indices gmapi to map ERA5 --> UFS mesh025 grid
  Here, nearest neighbour is used, fast and good enough
  Save as numpy binary npy / npz file or netcdf

  For this script to work, 6-hr GDAS atm output should be downloaded 
  and averaged by days, saved as numpy binary

  For deriving several daily means, with fetching GDAS from HPSS, prior to interpolation:
  use sbatch derive_dailyT2m_Ndays.sh --sdate 20250701 --edate 20250731 --dt 3
  This will save daily felds at 3 day intervals for 2025/07

  ML models developed / trained on PPAN

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray as xr
from yaml import safe_load
from mpl_toolkits.basemap import Basemap, cm
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
from mod_misc1 import dist_sphcrd
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_regmom as mrmom

def check_era_mesh025_coord(LATE, LONE, IE, JE, hlon, hlat, IG, JG, dmax=27e3, vocal=False):
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
    if vocal and (k > 0 and k % 20000) == 0:
      prc = k/len(JG)*100.
      print(f"  {prc:.2f}% processed, k={k} dist = {DD:.4f} m")

    DD_max = np.max([DD_max, DD])

    if DD > dmax:
      print(f"ERA point is {DD}m apart from mesh025, k={k}")
      print(f"ERA lon={lon_era} lat={lat_era}")
      print(f"mesh025 lon={lon_glr} lat={lat_glr}")
      raise Exception("ERR: check not passed")

  print(f"Checked gmapi: OK, overall max dist = {DD_max} m")

def read_daily_era5(dfl, dnmb0):
  """
  Derive daily fields from ERA5 hourly
  """
  yr, mm, dd = mtime.datevec(dnmb0)[:3]
  rdate0 = yr*10000 + mm*100 + dd

  with xr.open_dataset(dfl) as ds:
    # Check that all 24 hourly records are present
    if ds.sizes["Time"] != 24:
      raise ValueError(f"Expected 24 hourly records, found {ds.sizes['Time']}")

    # Check that all records belong to the requested day rday0=YYYYMMDD
    dates = ds.time.dt.strftime("%Y%m%d").astype(int).values

    if not np.all(dates == rdate0):
      raise ValueError(
          f"Time mismatch: requested {rdate0}, "
          f"file contains {np.unique(dates)}"
      )

    # Daily averages
    T2m_mean = ds["temperature"].mean(dim="Time", skipna=True).values
  
    # Convert K --> C if needed:
    if np.max(T2m_mean) > 270 and np.min(T2m_mean) > 150:
      T2m_mean -= 273.15

  return T2m_mean


parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--rdate", help="Date of ithkn prediction YYYMMDD", type=int, required=True)
parser.add_argument(
  "--field",
  help="Output field: ifdd (Freeze Degr. Days) or ihdd (Heat Degr. Days), or both",
  choices=['ifdd','ihdd','all'],
  required=True
)
args  = parser.parse_args()

regn       = args.regn
rdate      = args.rdate
fld_name   = args.field.lower()

if fld_name == 'all':
  save_ifdd = True
  save_ihdd = True
else:
  save_ifdd = fld_name == 'ifdd'
  save_ihdd = fld_name == 'ihdd'

dnmbR = mtime.rdate2datenum(rdate)
YRr, MMr, DDr = mtime.datevec(dnmbR)[:3]

fyaml = "config_MLithkn_oras5.yaml"
with open(fyaml) as ff:
  config_ml = safe_load(ff)

pthdata    = config_ml["MOM6"]["pthdata"]  # root dir for processed data
pthgmapi   = config_ml["ERA5"]["pthgmapi"]
regn_lat0  = config_ml["regn"][regn]["lat_bnd"]
dxy        = config_ml["params"]["dxy"]             # ice length scale
Tfrz       = config_ml["params"]["Tfrz"]
intgr_frz  = config_ml["params"]["intgr_frz"]
#tstep_era  = config_ml["params"]["tstep"]     # ERA5 daily fields, skipping days
sqrt_frzdays = config_ml["params"]["sqrt_frz"]   # use sqrt(IFDD) or not 

intgr_frz = int(intgr_frz)

# Construct time array for integrating IFDD and IHDD
dnmbR = mtime.rdate2datenum(rdate)
YRr, MMr, DDr = mtime.datevec(dnmbR)[:3]
dnmbS = int(dnmbR - intgr_frz + 1)
dnmbE = int(dnmbR)
DNMB = np.arange(dnmbS, dnmbE+1)

# Include last day:
#if DNMB[-1] != dnmbE:
#  DNMB = np.append(DNMB, dnmbE)

# Mapping index pairs mesh025 - ERA5:
dflgmapi = os.path.join(pthgmapi, 'gmapi_ERA5_to_UFSmesh025_global.nc')
with xr.open_dataset(dflgmapi) as ds:
  LONE = ds["era_longit"].values
  LATE = ds["era_latit"].values
  IERA = ds["era_indx"].values.astype(int)
  JERA = ds["era_jndx"].values.astype(int)
  I025 = ds["mesh025_indx"].values.astype(int)
  J025 = ds["mesh025_jndx"].values.astype(int)

# mesh025 grid and mask
# Read mesh025 grid:
pthgrid = config_ml["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xr.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, 1., -HH)
jdm, idm = HH.shape
#LMsk = HH < 0
LMsk = np.where(HH<0, 1, 0)

# Mask out not needed latitudes:
if regn == 'south':
  LMsk = np.where(hlat > regn_lat0, 0, LMsk)
elif regn == 'north':
  LMsk = np.where(hlat < regn_lat0, 0, LMsk)
else:
  raise RuntimeError(f"regn not matched: {regn}")

IFDD = None
IHDD = None
print("Deriving IFDD and IHDD")
for irec0, dnmb0 in enumerate(DNMB):
  YR, MM, DD = mtime.datevec(dnmb0)[:3]

  print(f"{irec0}:   Processing ERA5 T2m {YR}/{MM:02d}/{DD:02d}")
  pthera = config_ml["ERA5"]["pthera"].format(YR=YR)
  flera = f"ERA5_forcing_{YR}-{MM:02d}-{DD:02d}.nc"
  dflera = os.path.join(pthera, flera)
  if not os.path.isfile(dflera):
    raise FileNotFoundError(f"Not found: {dflera}")
 
  # Note N pole is fliped, lat goes from 90 to -90 
  T2m = read_daily_era5(dflera, dnmb0) 

  #  Use Zubov definition of accumulated freeze degree days
  #  sum(Tfrz-T), when SAT < Tfrz
  # = 0 otherwise
  # Similarly, integr. heat degree days
  freeze_degree = np.maximum(Tfrz - T2m, 0.0)
  heat_degree = np.maximum(T2m - Tfrz, 0.0)

  if IFDD is None:
    IFDD = freeze_degree.copy()
  else:
    IFDD += freeze_degree

  if IHDD is None:
    IHDD = heat_degree.copy()
  else:
    IHDD += heat_degree

# Process only active points, subsample 
Jsub, Isub = np.where(LMsk == 1)
JEsub = np.empty(len(Jsub), dtype=int)
IEsub = np.empty(len(Isub), dtype=int)

# Build a lookup table for fast finding corresponding ERA indices
print("Building remapping lookup table ...")
gmapi = {(jg,ig):(je,ie)
         for jg,ig,je,ie in zip(J025, I025, JERA, IERA)}

for k, (jj, ii) in enumerate(zip(Jsub, Isub)):
  if k > 0 and k % 50000 == 0:
    prc = k/len(Jsub)*100.
    print(f"  {prc:.2f}% processed")
  key = (int(jj), int(ii))
  if key not in gmapi:
    raise ValueError(f"Missing gmapi entry for UFS mesh025 index {key}")
  JEsub[k], IEsub[k] = gmapi[key]


check_gmapi = True
if check_gmapi:
  check_era_mesh025_coord(LATE, LONE, IEsub, JEsub, hlon, hlat, Isub, Jsub)


if save_ifdd: 
  print("Interpolating IFDD --> mesh025")

  IFDD_mesh025 = IFDD[JEsub, IEsub]

  assert np.all(np.isfinite(IFDD_mesh025)), "Non-finite IFDD_mesh025 values at active grid points"

  pthprd = config_ml["PRED"]["pthprd"]
  flfrz0 = config_ml["PRED"]["flifdd"].format(Ndays=intgr_frz, dxy=dxy, rdate=rdate, regn=regn)

  # Save numpy binary, only processed grid points
  flfrz = f"{flfrz0}.npz"
  dfliceout = os.path.join(pthprd, flfrz)
  print(f"Saving interpoalted IFDD --> {dfliceout}")
  np.savez(dfliceout, IFDD=IFDD_mesh025, JG=Jsub, IG=Isub, jdim=jdm, idim=idm)

if save_ihdd:
  # Interpolate to mesh025
  print("Interpolating IHDD --> mesh025")

  IHDD_mesh025 = IHDD[JEsub, IEsub]

  assert np.all(np.isfinite(IHDD_mesh025)), "Non-finite IHDD_mesh025 values at active grid points"

  pthprd = config_ml["PRED"]["pthprd"]
  flheat0 = config_ml["PRED"]["flihdd"].format(Ndays=intgr_frz, dxy=dxy, rdate=rdate, regn=regn)

  # Save numpy binary
  flheat = f"{flheat0}.npz"
  dfliceout = os.path.join(pthprd, flheat)
  print(f"Saving interpoalted IHDD  --> {dfliceout}")
  np.savez(dfliceout, IHDD=IHDD_mesh025, JG=Jsub, IG=Isub, jdim=jdm, idim=idm)


f_check = False
if f_check:
  DNMB_saved = DNMB[iStart:]

  plt.ion()

  # Plot ERA5 SAT at the sample locations:

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
  ax1.contour(HH,[0], colors=[(0.7,0.7,0.7)])
  #ax1.contour(LMsk, [0.99], linestyles='solid', colors=[(0.5,0.8,1)])
  #ax1.scatter(IG, JG, s=5, color=(0.8,0.4,0), marker='.')

  plt_ifdd = False

  if plt_ifdd:
    AA = IFDD_mesh025
    rmax = 1300
    stl = f"IFDD, {rdate}"
  else:
    AA = IHDD_mesh025
    rmax = 500
    stl = f"IFDD, {rdate}"

  sc = ax1.scatter(
    Isub, Jsub,
    c=AA,
    cmap='jet',
    s=20,          # marker size
    vmin=0,        # optional color scale limits
    vmax=rmax
  )

  ax1.set_title(stl)

  plt.colorbar(sc, ax=ax1)






