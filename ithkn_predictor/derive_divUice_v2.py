"""
  Derive dynamic predictor: sea ice divergence averaged over some area
  GLORYS

  Two approaches are possible (using divergence theorem)
  - average spatial integral of div(U) * dA
  - average integral of U flux across the boundary: U*n*dl, n - is normal comp. to the segment
  boundary flux is prefereable (as it conserved divergence)
  1st approach is straight forward but not exact

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
from mod_mom6 import dx_dy

#from MyPython.mod_cice6_utils import change_base_template, flname_replace_date

def read_uvice(DIRS, rdate):
  pthu = os.path.join(DIRS["pthui"],f"{YR}")
  dflice = mglr.find_file(rdate, pthu)
  with xr.open_dataset(dflice) as dsice:
    A2d = dsice['usi'].isel(time=0).values.squeeze()

  # Treat nans as no ice grid cells
  U2d = np.nan_to_num(A2d, nan=0.0)

  pthv = os.path.join(DIRS["pthvi"],f"{YR}")
  dflice = mglr.find_file(rdate, pthv)
  with xr.open_dataset(dflice) as dsice:
    A2d = dsice['vsi'].isel(time=0).values.squeeze()

  # Treat nans as no ice grid cells
  V2d = np.nan_to_num(A2d, nan=0.0)

  return U2d, V2d

def calc_divU(dltI, dltJ, ii, jj, U2d, V2d, Acell, DX, DY):
  """
    Compute average div u ice over specified region
    Assuming output fieds are at the cell centers
    Units:
    Acell = m2,  DX=m, DY=m
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
parser.add_argument(
  "--regn",
  help="Region to process",
  choices=['north','south'],
  required=True
)
parser.add_argument(
  "--load",
  help="Load saved divU tmp file, continue from last record (1), start from time 0 (0 default)",
  choices=[0,1],
  default=0,
  type=int
)
parser.add_argument("--davrg", help="Averaging time (days) for divU", required=True, type=int)
args = parser.parse_args()

regn       = args.regn
davrg      = args.davrg
load_saved = args.load == 1

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

# Load parameters:
lat0      = config_predictor["regn"][regn]["lat_bnd"]
tstep     = config_predictor["params"]["tstep"]
dxy       = config_predictor["params"]["dxy"]
YS        = config_predictor["params"]["ys"]
YS        = config_predictor["params"]["ys"]
YE        = config_predictor["params"]["ye"]

#davrg     = 1       # use N daily average ice fields to calc divU
dump_tstp = 50


varnmu = 'usi'
varnmv = 'vsi'

DIRS = {
  "pthithkn" : config_predictor["linregr"]["pthithkn"],
  "pthiconc" : config_predictor["linregr"]["pthiconc"],
  "pthsst"   : config_predictor["linregr"]["pthsst"],
  "pthssh"   : config_predictor["linregr"]["pthssh"],
  "pthui"    : config_predictor["linregr"]["pthui"],
  "pthvi"    : config_predictor["linregr"]["pthvi"],
  "ptht2m"   : config_predictor["linregr"]["ptht2m_1hr"],
  "pthout"   : config_predictor["linregr"]["pthout"],
  "ithkntmp" : config_predictor["linregr"]["ithkntmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "iconctmp" : config_predictor["linregr"]["iconctmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "ssttmp"   : config_predictor["linregr"]["ssttmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "divutmp"  : config_predictor["linregr"]["divutmp"].format(
        YS=YS, YE=YE, dxy=dxy, davrg=davrg, regn=regn),
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
DX, DY = dx_dy(hlon, hlat)
Acell = DX*DY


# Construct predictor sst time series for all locations, 
# Or load previously saved
pthout = DIRS["pthout"]
fltmp = DIRS["divutmp"]
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

#YY, JG, IG = construct_ifld(DNMB, DIRS, JG, IG, dfltmp, irec_start, YY, Acell, DX, DY, dxy, dump_tstp=10)

#Construct time series of response variable (ithkn)
#2D: locations x time
print(f"Start derivation of divU ice, time averaging = {davrg} days")

npnts = len(JG)
nrecs = len(DNMB)
dump_tstp = 10
if YY is None:
  YY = np.zeros((npnts, nrecs), dtype=float)*np.nan
for irec, dnmb0 in enumerate(DNMB):
  YR, MM, DD = mtime.datevec(dnmb0)[:3]

  if irec < irec_start:
    print(f"Skipping, already processed: {YR}/{MM:02d}/{DD:02d}")
    continue
  
  print(f"Deriving divUice for {YR}/{MM:02d}/{DD:02d}")

  rdate = int(YR*1e4 + MM*100 + DD)

  # If davrg > 1: average ice velocity fields
  if davrg > 1:
    Usum = np.zeros(DX.shape)
    Vsum = np.zeros(DX.shape)
    icc = 0

    # No date before 1993:
    d93 = mtime.datenum([1993,1,1])
    day1 = max(d93, dnmb0 - davrg + 1)
    day2 = dnmb0 + 1
    for dnmbA in range(day1, day2):
      YR, MM, DD = mtime.datevec(dnmbA)[:3]
      rdate = int(YR*1e4 + MM*100 + DD)  
      U, V = read_uvice(DIRS, rdate)
      icc += 1
      Usum += U
      Vsum += V

    U2d = Usum / icc
    V2d = Vsum / icc
  else:
    U2d, V2d = read_uvice(DIRS, rdate)

  fld_pnts = []
  for jj, ii in zip(JG, IG):
    # Estimate box size based on min distance criterion
    dltX = DX[jj,ii]*1e-3  # km
    dltY = DY[jj,ii]*1e-3  # km
    dltI = int(np.ceil(dxy / dltX))
    dltJ = int(np.ceil(dxy / dltY))
    div_uice = calc_divU(dltI, dltJ, ii, jj, U2d, V2d, Acell, DX, DY)
    fld_pnts.append(div_uice)

  YY[:,irec] = np.asarray(fld_pnts)

  assert YY.shape[0] == len(IG), f"Check YY shape does not match IG {YY.shape}"

  if (irec + 1) % dump_tstp == 0: 
    print(f"TMP step: Saving divUice time series  --> {dfltmp}")
    np.savez(dfltmp,
           YY=YY,
           JG=JG,
           IG=IG,
           DNMB=DNMB)

# Save:
print(f"Final Saving divUice time series and IG, JG --> {dfltmp}")
np.savez(dfltmp,
       YY=YY,
       JG=JG,
       IG=IG,
       DNMB=DNMB)


f_check = False
if f_check:
  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
  ax1.contour(LMsk, [0.99], linestyles='solid', colors=[(0.5,0.8,1)])
  ax1.scatter(IG, JG, s=5, color=(0.8,0.4,0), marker='.')

  # More points in S. Ocean:
  JB,IB = np.where((np.abs(U2d) > 0) & (hlat < -55))

  # Plot stereogr porj
  # Plotting
  import mod_colormaps as mclrmps

  cff = 1.e5
  DIVU = np.zeros_like(LMsk, dtype=float)
  DIVU[JG,IG] = np.asarray(fld_pnts) * cff
  DIVU[JB,IB] = np.asarray(fld_pnts) * cff
  #DIVU = U2d.copy()
  DIVU[LMsk == 0] = np.nan
  
  clrmp = mclrmps.colormap_uv()
  rmin = -0.1
  rmax = 0.1
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  clrmp.set_under(color=[1,1,1])
  
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])

    
  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l', ax=ax1)
    parallels = np.arange(-80,-10,10.)
    meridians = np.arange(-360,359.,45.)

    # Subset region
    JJ = np.where(hlat[:, 0] <= -50)[0]

  elif regn == 'north':
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l', ax=ax1)
    parallels = np.arange(50, 90, 5)
    meridians = np.arange(-360, 359., 45.)

    # Subset region
    JJ = np.where(hlat[:, 0] >= 50)[0]

  hlat_s = hlat[JJ, :]
  hlon_s = hlon[JJ, :]
  AP_s   = DIVU[JJ, :]


  xh, yh = m(hlon_s, hlat_s)

  m.drawparallels(parallels, labels=[0,0,0,0])
  m.drawmeridians(meridians, labels=[0,0,0,0])
  m.drawcoastlines()

  img = m.pcolormesh(xh,yh, AP_s, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"divU {YR}/{MM:02d}/{DD:02d}")

  ax3 = fig1.add_axes([0.2, 0.06, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = f' @{machine}: derive_divUice_v2'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)

  
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

  xx = np.arange(iS,iE+1)
  yy = np.arange(jS,jE+1)
  uu = U2d[jS:jE+1,iS:iE+1]
  vv = V2d[jS:jE+1,iS:iE+1]

  # Contour of integration
  ax1.plot([iS,iE],[jE,jE],'c-')
  ax1.plot([iS,iE],[jS,jS],'c-')
  ax1.plot([iE,iE],[jS,jE],'c-')
  ax1.plot([iS,iS],[jS,jE],'c-')
  ax1.quiver(xx, yy, uu, vv, scale=10)



