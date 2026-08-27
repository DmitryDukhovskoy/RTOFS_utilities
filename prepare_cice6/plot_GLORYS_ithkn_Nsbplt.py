"""
  GLORYS daily ice thickness fields 
  interpolated onto mesh025 
  see PPAN: ../ithkn_predictor/interp_GLORYSithkn_to_mesh025_month.py
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

regn = 'north'

parser = argparse.ArgumentParser()
parser.add_argument("--rdate", help="Prediction date YYYYMMDD",
                    required=True,
                    nargs="+",
                    type=int)
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str, default='north')
parser.add_argument("--ncol", help=f"N of columns for subplots, default=N dates", type=int)
parser.add_argument("--nrow", help=f"N of rows for subplots, default=1", type=int)
args = parser.parse_args()

RDATES  = args.rdate
regn  = args.regn
ncol    = args.ncol if args.ncol is not None else len(args.rdate)
nrow    = args.nrow if args.nrow is not None else 1


syst_info = os.uname()
machine = syst_info.nodename

if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
elif 'an' in machine:
  print("Running on PPAN node:", machine)
  node_nm = "ppan"
else:
  print("Unknown machine:", machine)

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xr.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

pthglr = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/GLORYS_ithkn_interp_UFSmesh025'

def read_icefld(rdate, dfglr, HH):
  dnmbR = mtime.rdate2datenum(rdate*100)  # restart day nmb
  YR, MM, DD = mtime.datevec(dnmbR, round_hrs=True)[:3]
  with xr.open_dataset(dfglr) as dcice:
    for varnm in ['hi_h', 'hi_d', 'ice_thkn']:
      if varnm in dcice:
        A2d = dcice[varnm].isel(time=DD-1).values.squeeze()
        break
    else:
      print(list(dcice.data_vars))
      raise KeyError(f"No aice_h or aice_d or ice_thkn in {dfglr}")

  A2d[HH >= 0] = np.nan

  return A2d

def plot_field(ax1, m, xh, yh, AP_s, clrmp, rmin, rmax, regn, sttl):
  plt.sca(ax1)
  if regn == 'north':
    m.drawparallels(np.arange(60, 90, 10), labels=[0,0,0,0])
  elif regn == 'south':
    m.drawparallels(np.arange(-80, -50, 10), labels=[0,0,0,0])

  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,0])
  m.drawcoastlines()

  img = ax1.pcolormesh(xh, yh, AP_s, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(sttl, fontsize=8)

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
  flglr = f"GLORYS_ithkn_interp_mesh025_{YR}{MM:02d}_{regn}.nc"
  dfglr = os.path.join(pthglr, flglr)
  if not os.path.isfile(dfglr):
    raise RuntimeError(f"File not found: {dfglr}")

  sttl = f"ithkn GLORYS mesh025\n{YR}/{MM:02d}/{DD:02d}"

  irow = iplt // ncol
  icol = iplt % ncol
  ax1 = axes[irow, icol]

  print(f"Reading {YR}/{MM}/{DD}: {dfglr}")

  A2d = read_icefld(rdate, dfglr, HH)
  img = plot_field(ax1, m, xh, yh, A2d, clrmp, rmin, rmax, regn, sttl)


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
clb.ax.tick_params(direction='in', length=12, labelsize=12)

fig1.canvas.draw()

pos_clb = ax3.get_position()
bot_clb = pos_clb.y0
pbtm = bot_clb - 0.05

btx = f' @{machine}: plot_GLORYSmesh025_ithkn_Nsbplt.py'
#bottom_text(btx, pos=[0.02,0.02], fsz=8)
bottom_text(btx, pos=[0.02, pbtm], fsz=8)



