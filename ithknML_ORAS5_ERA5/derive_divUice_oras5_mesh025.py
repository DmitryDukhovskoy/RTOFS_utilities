"""
  Derive dynamic predictor: sea ice divergence averaged over some area
  from ORAS5 fields

  ORAS5 fields derived by Phil

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
import glob

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
from mod_misc1 import dist_sphcrd
import mod_time as mtime
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_regmom as mrmom
from mod_mom6 import dx_dy
import mod_ml_emulator as mml


def calc_divU(dltI, dltJ, ii, jj, U2d, V2d, Acell, DX, DY):
  """
    Compute average div u ice over specified region
    Assuming output fieds are at the cell centers
  """
  divU = 0
  jdm, idm = U2d.shape

  Intgr = 0.0
  # Define the box around the grid point:
  iS = int(ii - dltI)
  iE = int(ii + dltI)
  iS = np.max([iS, 0])
  iE = np.min([iE, idm-1])

  jS = int(jj - dltJ)
  jE = int(jj + dltJ)
  # Boundaries, better - for global grid
  # use grid points at the opposite side
  jS = np.max([0, jS])
  jE = np.min([jE, jdm-1])

  # Integrate along the boundary:
  # Note different sign of outward norm vectors 
  # along the box sides
  Intgr = (
    np.sum(U2d[jS:jE+1, iE] * DY[jS:jE+1, iE])
    - np.sum(U2d[jS:jE+1, iS] * DY[jS:jE+1, iS])
    + np.sum(V2d[jE, iS:iE+1] * DX[jE, iS:iE+1])
    - np.sum(V2d[jS, iS:iE+1] * DX[jS, iS:iE+1])
  )

  # subtract half of the four corner contributions
  Intgr -= 0.5 * (
    U2d[jS, iE] * DY[jS, iE]
    + U2d[jE, iE] * DY[jE, iE]
    - U2d[jS, iS] * DY[jS, iS]
    - U2d[jE, iS] * DY[jE, iS]
    + V2d[jE, iS] * DX[jE, iS]
    + V2d[jE, iE] * DX[jE, iE]
    - V2d[jS, iS] * DX[jS, iS]
    - V2d[jS, iE] * DX[jS, iE]
  )

  # Space-Average divergence:
  Area = np.sum(Acell[jS:jE+1, iS:iE+1])
  assert Area > 0, f"ii={ii}, jj={jj}, dltI={dltI}, dltJ={dltJ}, Area = {Area}"
  div_uice = Intgr / Area

  return div_uice

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--rdate", help="Date of ithkn prediction YYYMMDD", type=int, required=True)
args  = parser.parse_args()

regn       = args.regn
rdate      = args.rdate

dnmbR = mtime.rdate2datenum(rdate)
YRr, MMr, DDr = mtime.datevec(dnmbR)[:3]

fyaml = "config_MLithkn_oras5.yaml"
with open(fyaml) as ff:
  config_ml = safe_load(ff)

pthdata    = config_ml["MOM6"]["pthdata"]  # root dir for processed data
regn_lat0  = config_ml["regn"][regn]["lat_bnd"]
dxy        = config_ml["params"]["dxy"]             # ice length scale

# mesh025 grid and mask
pthgrid = config_ml["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
hlon = (hlon + 360) % 360

with xr.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, 100., -HH)
jdm, idm = HH.shape
LMsk = HH < 0

# Set mask for region and ML lat bounds:
assert HH.shape == hlat.shape, f"Shape mismatch: HH={HH.shape}, hlat={hlat.shape}"

if regn == 'north':
  LMsk &= hlat >= regn_lat0
elif regn == 'south':
  LMsk &= hlat <= regn_lat0

# Mapping index pairs mesh025 - ORAS5:
pthgmapi = config_ml["ORAS5"]["pthgmapi"]
flgmapi  = config_ml["ORAS5"]["flgmapi"]
dflgmapi = os.path.join(pthgmapi, flgmapi)
with xr.open_dataset(dflgmapi) as ds:
  LONOR = ds["oras_longit"].values
  LATOR = ds["oras_latit"].values
  IORAS = ds["oras_indx"].values.astype(int)
  JORAS = ds["oras_jndx"].values.astype(int)
  I025  = ds["mesh025_indx"].values.astype(int)
  J025  = ds["mesh025_jndx"].values.astype(int)

LONOR = (LONOR + 360) % 360
DX, DY = dx_dy(LONOR, LATOR)
Acell_oras = DX*DY

# Create a look-up table for matching mesh025 -> ORAS5 index pairs
# Subsample mesh025 grid points to be processed
Jsub, Isub = np.where(LMsk == 1)
JORsub = np.empty(len(Jsub), dtype=int)
IORsub = np.empty(len(Isub), dtype=int)

# Build a lookup table for fast finding corresponding ORAS5 indices
print("Building remapping lookup table ...")
gmapi = {(jg,ig):(je,ie)
         for jg,ig,je,ie in zip(J025, I025, JORAS, IORAS)}

for k, (jj, ii) in enumerate(zip(Jsub, Isub)):
  if k > 0 and k % 50000 == 0:
    prc = k/len(Jsub)*100.
    print(f"  {prc:.2f}% processed")
  key = (int(jj), int(ii))
  if key not in gmapi:
    raise ValueError(f"Missing gmapi entry for UFS mesh025 index {key}")
  JORsub[k], IORsub[k] = gmapi[key]

check_gmapi = True
if check_gmapi:
  mml.check_oras_mesh025_coord(LATOR, LONOR, IORsub, JORsub, hlon, hlat, Isub, Jsub)


# Read ORAS5 ice fields:
# Find file
pthoras = config_ml["ORAS5"]["pthoras"]
dfloras = mml.find_oras_ice_file(pthoras, rdate)

with xr.open_dataset(dfloras) as dsice:
  # Read time in YYYYMMDD format
  Time = dsice["time_counter"].dt.strftime("%Y%m%d").astype(int).values 
  if not np.any(Time == rdate):
    raise ValueError(
        f"Time mismatch: requested {rdate}, "
        f"file contains {Time}"
    )

  itime = np.where(Time == rdate)[0][0]
  U2d = dsice["iicevelu"].isel(time_counter=itime).values
  V2d = dsice["iicevelv"].isel(time_counter=itime).values


# Treat nans as no ice grid cells
U2d = np.nan_to_num(U2d, nan=0.0)
V2d = np.nan_to_num(V2d, nan=0.0)

assert U2d.shape == V2d.shape, (
  f"U/V shape mismatch: U2d={U2d.shape}, V2d={V2d.shape}"
)

assert U2d.shape == Acell_oras.shape == DX.shape == DY.shape, (
  f"Grid shape mismatch: U2d={U2d.shape}, "
  f"Acell={Acell_oras.shape}, DX={DX.shape}, DY={DY.shape}"
)

# Calc divUice on ORAS5 grid
fld_pnts = []
icc = 0
npnts = len(JORsub)
ichck = int(npnts * 0.05)
print("Calculating divU ice:")
for jj, ii in zip(JORsub, IORsub):
  icc += 1
  nprc = icc / npnts * 100.
  if icc % ichck == 0 or icc == npnts:
    print(f"   processed {nprc:.2f}% ...")
  # Estimate box size based on min distance criterion
  dltX = DX[jj,ii]*1e-3  # km
  dltY = DY[jj,ii]*1e-3  # km
  dltI = int(np.ceil(dxy / dltX))
  dltJ = int(np.ceil(dxy / dltY))
  div_uice = calc_divU(dltI, dltJ, ii, jj, U2d, V2d, Acell_oras, DX, DY)
  fld_pnts.append(div_uice)
print(" ----- DONE ----")

# Selected ORAS points match mesh025 analysis points
DIVUI = np.asarray(fld_pnts)
assert DIVUI.shape[0] == npnts, f"Check DIVUI shape {DIVUI.shape}"

# Save numpy binary
pthprd   = config_ml["PRED"]["pthprd"]
flnm0    = config_ml["PRED"]["fldivui"].format(dxy=dxy, rdate=rdate, regn=regn)

fldivui = f"{flnm0}.npz"
dfliceout = os.path.join(pthprd, fldivui)

print(f"Final Saving divUice and Isub, Jsub --> {dfliceout}")
np.savez(dfliceout, DIVUI=DIVUI, J025=Jsub, I025=Isub, jdim=jdm, idim=idm)


f_check = False
if f_check:
  cf = 1e6

  # Insert analysis points into mesh025
  DU = np.zeros_like(LMsk, dtype=float)
  DU[Jsub, Isub] = DIVUI * cf
    
  # Land mask:
  DU[HH>=0] = np.nan

  clrmp = mclrmps.colormap_temperature_coldwarm()
  rmin = -3. 
  rmax = 3.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()         
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])

  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
    parallels = np.arange(-80,-10,10.)
    meridians = np.arange(-360,359.,45.)
  elif regn == 'north':
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
    parallels = np.arange(50, 90, 5)
    meridians = np.arange(-360, 359., 45.)
        
  xh, yh = m(hlon, hlat)

  m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
  m.drawcoastlines()

  img = m.pcolormesh(xh,yh, DU, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"divUice * {cf:.0e} s-1, ORAS5,  mesh025, {YR}/{MM:02d}/{DD:02d}")

  ax3 = fig1.add_axes([0.2, 0.08, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='both')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = 'derive_divUice_oras5_mesh025.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)

