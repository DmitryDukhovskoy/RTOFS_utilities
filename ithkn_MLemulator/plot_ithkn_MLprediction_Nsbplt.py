"""
  Plot ML preidcted ithkn fields
  for N specified dates
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray as xr
import matplotlib.colors as colors
from yaml import safe_load
from mpl_toolkits.basemap import Basemap, cm
import argparse
from mpl_toolkits.basemap import Basemap, cm

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
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])

from mod_utils_fig import bottom_text, colorbar_horiz
import mod_time as mtime
import mod_colormaps as mclrmps
import mod_mom6 as mmom6

regn = 'north'

parser = argparse.ArgumentParser()
parser.add_argument("--rdate", help="Prediction date YYYYMMDD",
                    required=True,
                    nargs="+",
                    type=int)
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str, default='north')
parser.add_argument("--ncol", help=f"N of columns for subplots, default=N dates", type=int)
parser.add_argument("--nrow", help=f"N of rows for subplots, default=1", type=int)
parser.add_argument("--model", help="ML model to use",
                   choices=['rf1','rf2','gbr2','gbr3'],
                   required=True)
args = parser.parse_args()

RDATES  = args.rdate
regn  = args.regn
ncol    = args.ncol if args.ncol is not None else len(args.rdate)
nrow    = args.nrow if args.nrow is not None else 1
ml_model  = args.model

fyaml = "paths_ML.yaml"
with open(fyaml) as ff:
  config_ml = safe_load(ff)

# Prediction grid points:
# Get MOM6 grid:
pthgrid    = config_ml["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xr.open_dataset(dftopo_mom) as dstopo:
  depth = dstopo['depth'].data.squeeze()

# Convert all positive values -> land (100) and ocean (<0):
HH = np.where(depth < 1.e-20, 100., -depth)

jdm, idm = HH.shape
LMsk = HH < 0


def plot_field(ax1, m, xh, yh, A2d, clrmp, rmin, rmax, sttl, regn):
  if regn == 'north':
    m.drawparallels(np.arange(60, 90, 10), labels=[0,0,0,0])
  elif regn == 'south':
    m.drawparallels(np.arange(-80, -50, 10), labels=[0,0,0,0])
  
  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,0])
  m.drawcoastlines()
  
  img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)

  ax1.set_title(sttl, fontsize=12)

  return img

clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.

clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north': 
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)

xh, yh = m(hlon,hlat) # GFS coords


print("Plotting ...")
plt.ion()
fig1 = plt.figure(1,figsize=(12, 9))
fig1.clf()  # Clear the figure
axes = fig1.subplots(nrows=nrow, ncols=ncol, squeeze=False)

fig1.subplots_adjust(
    left=0.02,
    right=0.99,
    top=0.95,
    bottom=0.1,
    wspace=0.05,
    hspace=0.1
)

for iplt, rdate in enumerate(RDATES):
  dnmbR = mtime.rdate2datenum(rdate*100)  # restart day nmb
  YR, MM, DD = mtime.datevec(dnmbR, round_hrs=True)[:3]

  # Load prediction:
  pthdump = config_ml["PRED"]["pthdump"]
  flfcst = f"{ml_model}_ithkn_fcast_GFSv17anls_{rdate}_{regn}.npz"
  dflfcst = os.path.join(pthdump, flfcst)
  print(f"Loading fcst: {dflfcst}")

  FCST = np.load(dflfcst)
  Yfcst = FCST["Yfcst"]
  JG    = FCST["JG"]
  IG    = FCST["IG"]
  idim_s = FCST["idim"]
  jdim_s = FCST["jdim"]

  assert HH.shape == (jdim_s,idim_s), (
    f"ML field for different grid? expected j x i: {jdim_s} {idim_s}"
  )

  A2d = np.zeros_like(HH) * np.nan
  A2d[LMsk] = 0
  A2d[JG,IG] = Yfcst


  sttl = f"ML {ml_model} ithkn using SOCA CICE6\n{YR}/{MM:02d}/{DD:02d}"

  irow = iplt // ncol
  icol = iplt % ncol
  ax1 = axes[irow, icol]

  img = plot_field(ax1, m, xh, yh, A2d, clrmp, rmin, rmax, sttl, regn)


#Plot colrbar
# Force Matplotlib to finalize axes positions
fig1.canvas.draw()

# Get FINAL positions of the subplot axes
positions = [
    ax.get_position()
    for ax in axes.flat
    if ax.get_visible()
]

left   = min(p.x0 for p in positions)
right  = max(p.x1 for p in positions)
bottom = min(p.y0 for p in positions)

# Colorbar dimensions
cb_height = 0.016
cb_gap = 0.016

pos_clrb = [
    left,
    bottom - cb_gap - cb_height,
    right - left,
    cb_height
]

# Create colorbar
ax3 = fig1.add_axes(pos_clrb)

clb = fig1.colorbar(
    img,
    cax=ax3,
    orientation='horizontal',
    extend='max'
)

#clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ticks = np.linspace(rmin, rmax, 11)
clb.set_ticks(ticks)
clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
clb.ax.tick_params(direction='in', length=10, labelsize=12)

fig1.canvas.draw()

pos_clb = ax3.get_position()
bot_clb = pos_clb.y0
pbtm = bot_clb - 0.05

btx = f'plot_ithkn_MLprediction_Nsbplt.py'
bottom_text(btx, pos=[0.02, pbtm], fsz=8)

