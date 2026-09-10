"""
  Derive dynamic predictor: sea ice divergence averaged over some area
  from GDAS / SOCA sea ice fields

  Fields from HPSS, use scripts to fetch fields
  get_soca_ice_restart.sh

  Two approaches are possible (using divergence theorem)
  - average spatial integral of div(U) * dA
  - average integral of U flux across the boundary: U*n*dl, n - is normal comp. to the segment
  boundary flux is prefereable (as it conserved divergence)
  1st approach is straight forward but not exact
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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_sis2_relax as msisrlx
import mod_regmom as mrmom
from mod_mom6 import dx_dy

parser = argparse.ArgumentParser()
parser.add_argument("--dxy", 
  help=f"Min dist (km) between data points (~corr.scale), to skip close i,j points, default=50",
  type=int, 
  default=50)
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--rdate", help="Date of ithkn prediction YYYMMDD", type=int, required=True)
args  = parser.parse_args()

dxy        = args.dxy
regn       = args.regn
rdate      = args.rdate

fyaml = 'paths_ML.yaml'
with open(fyaml) as ff:
  pths_ml = safe_load(ff)

dnmbR = mtime.rdate2datenum(rdate)
YRr, MMr, DDr = mtime.datevec(dnmbR)[:3]

pthdata   = pths_ml["GDAS"]["pthdata"]  # root dir for processed data
pthice    = pths_ml["GDAS"]["pthsoca_ice"].format(YR=YRr, MM=MMr, DD=DDr)
regn_name = pths_ml[regn]["name"]
regn_lat0 = pths_ml[regn]["lat0"]  # bounding lat for ML emulator

# Get MOM6 grid:
pthgrid    = pths_ml["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
     
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
DX, DY = dx_dy(hlon, hlat)
Acell = DX*DY

with xr.open_dataset(dftopo_mom) as dstopo:
  depth = dstopo['depth'].data.squeeze()

# Convert all positive values -> land (100) and ocean (<0):
HH = np.where(depth < 1.e-20, 100., -depth)
 
jdm, idm = HH.shape
LMsk = HH < 0

# Set mask for region and ML lat bounds:
assert HH.shape == hlat.shape, f"Shape mismatch: HH={HH.shape}, hlat={hlat.shape}"

if regn == 'north':
  LMsk &= hlat >= regn_lat0
elif regn == 'south':
  LMsk &= hlat <= regn_lat0

# Grid points to process:
JG, IG = np.where(LMsk)

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


# Read SOCA ice fields:
flice = pths_ml["GDAS"]["flsoca_ice"].format(YR=YRr, MM=MMr, DD=DDr)
dflice = os.path.join(pthice, flice)

assert os.path.isfile(dflice), "File not found: {dflice}"
with xr.open_dataset(dflice) as dsice:
  U2d = dsice["uvel"].values
  V2d = dsice["vvel"].values

# Treat nans as no ice grid cells
U2d = np.nan_to_num(U2d, nan=0.0)
V2d = np.nan_to_num(V2d, nan=0.0)

assert U2d.shape == V2d.shape, (
    f"U/V shape mismatch: U2d={U2d.shape}, V2d={V2d.shape}"
)

assert U2d.shape == Acell.shape == DX.shape == DY.shape, (
    f"Grid shape mismatch: U2d={U2d.shape}, "
    f"Acell={Acell.shape}, DX={DX.shape}, DY={DY.shape}"
)


fld_pnts = []
icc = 0
npnts = len(JG)
ichck = int(npnts * 0.05)
print("Calculating divU ice:")
for jj, ii in zip(JG, IG):
  icc += 1
  nprc = icc / npnts * 100.
  if icc % ichck == 0 or icc == npnts:
    print(f"   processed {nprc:.2f}% ...")
  # Estimate box size based on min distance criterion
  dltX = DX[jj,ii]*1e-3  # km
  dltY = DY[jj,ii]*1e-3  # km
  dltI = int(np.ceil(dxy / dltX))
  dltJ = int(np.ceil(dxy / dltY))
  div_uice = calc_divU(dltI, dltJ, ii, jj, U2d, V2d, Acell, DX, DY)
  fld_pnts.append(div_uice)
print(" ----- DONE ----")

DIVUI = np.asarray(fld_pnts)
assert DIVUI.shape[0] == npnts, f"Check DIVUI shape does not match IG {DIVUI.shape}"

# Save numpy binary
pthprd   = pths_ml["PRED"]["pthprd"]
flnm0    = pths_ml["PRED"]["fldivui"].format(dxy=dxy, rdate=rdate, regn=regn)

fldivui = f"{flnm0}.npz"
dfliceout = os.path.join(pthprd, fldivui)

print(f"Final Saving divUice and IG, JG --> {dfliceout}")
np.savez(dfliceout, DIVUI=DIVUI, JG=JG, IG=IG, jdim=jdm, idim=idm)


f_check = False
if f_check:
  cf = 1e6
  DU = LMsk * 0.
  DU[JG,IG] = DIVUI * cf

  # Land mask:
  DU[HH>=0] = np.nan
  

  clrmp = mclrmps.colormap_temperature_coldwarm()
  rmin = -5. 
  rmax = 5.
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
  ax1.set_title(f"divUice * {cf:.0e} s-1, SOCA CICE6,  mesh025, {YR}/{MM:02d}/{DD:02d}")

  ax3 = fig1.add_axes([0.2, 0.08, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='both')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = 'derive_divUice_cice.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)


