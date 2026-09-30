"""
Define time steps and IJ GLORYS grid points
to construct predictors / input fileds and response field
for ML emulators
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
from mod_misc1 import dist_sphcrd
from mod_mom6 import dx_dy

#from MyPython.mod_cice6_utils import change_base_template, flname_replace_date

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True)
parser.add_argument(
  "--adjice",
  help="Remove grid points from ML learning that has no ice during year",
  choices=[0,1],
  default=1
)
args = parser.parse_args()

regn  = args.regn
f_adjice = args.adjice == 1 

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

regn_name = config_predictor["regn"][regn]["name"]
lat0      = config_predictor["regn"][regn]["lat_bnd"]
dxy       = config_predictor["params"]["dxy"]
tstep     = config_predictor["params"]["tstep"]
YS        = config_predictor["params"]["ys"]
YE        = config_predictor["params"]["ye"]

def update_icepnts(YY, JG, IG, hice_min=0.1, nprc=0.1):
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

DNMB = micpr.derive_time(YS, YE, tstep)

# Derive I,J grid points from GLORYS grid
# Read GLORYS grid:
pthice = os.path.join(config_predictor["linregr"]["pthiconc"],f"{YS}")
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
pthssh = os.path.join(config_predictor["linregr"]["pthssh"],f"{YS}")
dflssh = mglr.find_file(rdate, pthssh)
with xr.open_dataset(dflssh) as dszos:
  SSH = dszos['zos'].isel(time=0).values.squeeze()

LMsk = np.where(np.isfinite(SSH),1,0)

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

# Define domain:
if regn == 'north':
  DOMAIN = (hlat >= lat0) & (LMsk == 1)
  # Exclude N. Atlantic and Barents Sea:
  DOMAIN[:1893,2154:2795] = False
  DOMAIN[:1764,:258] = False   # N. Bering 
  DOMAIN[:1807,1991:2243] = False # Iceland Sea

elif regn == 'south':
  DOMAIN = (hlat <= lat0) & (LMsk == 1)

JG = []
IG = []
nj, ni = DOMAIN.shape
jj0 = np.argmax(DOMAIN.any(axis=1)) # index of the 1st row containing any valid grid pnt
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


# Use N last years to eliminate grid point 
# that are ice-free
if f_adjice:
  DNMBT = derive_time(YE-4, YE, tstep)
  npnts = len(JG)
  nrecs = len(DNMB)
  YY = np.zeros((npnts, nrecs), dtype=float)*np.nan

  for irec, dnmb0 in enumerate(DNMBT):
    YR, MM, DD = mtime.datevec(dnmb0)[:3]
    print(f"irec={irec} Reading ithkn {YR}/{MM:02d}/{DD:02d}")

    pthice = os.path.join(config_predictor["linregr"]["pthithkn"],f"{YR}")
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

  # Eliminate unneeded grid points with no ice
  JG0 = JG.copy()
  IG0 = IG.copy()
  if regn == 'south':
    YY, JG, IG = update_icepnts(YY, JG0, IG0, hice_min=0.05, nprc=0.05)
  elif regn == 'north':
    YY, JG, IG = update_icepnts(YY, JG0, IG0, hice_min=0.2, nprc=0.05)

  print(f"Original N grid points: {len(JG0)}, after ice-free eliminated: {len(JG)}\n")

# Save:
pthout = config_predictor["params"]["pthinfo"]
fltime = config_predictor["params"]["fltime"].format(regn=regn, tstep=tstep)
flij   = config_predictor["params"]["flij"].format(regn=regn)
flgrid = config_predictor["params"]["flgrid"]
dfltime = os.path.join(pthout, fltime)
dflij   = os.path.join(pthout, flij)
dflgrid = os.path.join(pthout, flgrid)

print(f"Saving TIME --> {dfltime}")
np.save(dfltime, DNMB)

print(f"Saving I,J grid points --> {dflij}")
np.savez(dflij, JG=JG, IG=IG)

print(f"Saving GLORYS grid --> {dflgrid}")
np.savez(dflgrid, LON=hlon, LAT=hlat, LMsk=LMsk)


f_check = False
if f_check:

  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
  ax1.contour(LMsk, [0.99], linestyles='solid', colors=[(0.5,0.8,1)])

  sc = ax1.scatter(
    IG, JG,
    c=SSH[JG,IG],
    cmap='jet',
    s=20,          # marker size
    vmin=-1,        # optional color scale limits
    vmax=1
 )

  plt.colorbar(sc, ax=ax1, label='SSH, m')




