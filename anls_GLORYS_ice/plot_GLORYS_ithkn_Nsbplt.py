"""
  GLORYS daily - both h/shperes

  Plot Interpolated CryoSatice thickn. monthly fileds 
  summer months only

I have a subset of monthly unfilled siconc and sithick for the Arctic on analysis:
/work1/tjc/datasets/glorys/GLOBAL_MULTIYEAR_PHY_001_030/monthly/not_filled/GLORYS_arctic.199301-202412.siconc.nc

The daily data is on uda, you can find it here:
/uda/Global_Ocean_Physics_Reanalysis/global/daily/siconc/

/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray
import matplotlib.colors as colors
from yaml import safe_load
from mpl_toolkits.basemap import Basemap, cm
import argparse

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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_time as mtime

regn = 'north'
field_name = 'ithkn'

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



# Find file:
def find_file(rdate, pthice):
  from pathlib import Path
  try:
    dflice = next(
        Path(pthice).glob(
            f"*_mean_{rdate}_R*.nc"
        )
    )
    print(f"Found file: {dflice}")
    return dflice
  except StopIteration:
    print(f"No file found for {rdate} in {pthice}")

def read_glorys(rdate, LMsk, hlat_s, hlon_s, JJ):
  dnmbR = mtime.rdate2datenum(rdate*100)  # restart day nmb
  YR, MM, DD = mtime.datevec(dnmbR, round_hrs=True)[:3]
  pthice = f"/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick/{YR}"
  #dflice = os.path.join(pthice, flice)
  dflice = find_file(rdate, pthice)
  if not os.path.isfile(dflice):
    raise RuntimeError(f"File not found: {dflice}")

  with xarray.open_dataset(dflice) as dsice:
    A2d = dsice['sithick'].isel(time=0).data.squeeze()
    LON = dsice['longitude'].values
    LAT = dsice['latitude'].values 


  if LMsk is None or hlat_s is None or hlon_s is None:
    pthssh = f"/uda/Global_Ocean_Physics_Reanalysis/global/daily/zos/{YR}"
    dflssh = find_file(rdate, pthssh)
    with xarray.open_dataset(dflssh) as dszos:
      SSH = dszos['zos'].isel(time=0).values.squeeze()

    LMsk = np.where(np.isfinite(SSH),1,0)

    # Subset region
    hlon, hlat = np.meshgrid(LON, LAT)
    if regn == 'south':
      JJ = np.where(hlat[:, 0] <= -50)[0]
    elif regn == 'north':
      JJ = np.where(hlat[:, 0] >= 50)[0]

    hlat_s = hlat[JJ, :]
    hlon_s = hlon[JJ, :]

  AP = A2d.copy()
  if LMsk is not None:
    Jocn = (LMsk == 1) & (~np.isfinite(AP))   # open ocean
    AP[Jocn] = 0.

  AP_s   = AP[JJ, :]

  return AP_s, hlat_s, hlon_s, LMsk, JJ

clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.2, 0.2, 0.2])
clrmp.set_under(color=[1,1,1])


#m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l', ax=ax1)
if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)



def plot_field(ax1, fig1, m, xh, yh, AP_s, clrmp, rmin, rmax, regn, dnmb):
  YR, MM, DD = mtime.datevec(dnmb)[:3]
  plt.sca(ax1)
  if regn == 'north':
    m.drawparallels(np.arange(60, 90, 10), labels=[0,0,0,0])
  elif regn == 'south':
    m.drawparallels(np.arange(-80, -50, 10), labels=[0,0,0,0])

  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,0])
  m.drawcoastlines()

  img = ax1.pcolormesh(xh, yh, AP_s, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"ithkn GLORYS\n{YR}/{MM:02d}/{DD:02d}",fontsize=8)

  return img

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

xh = yh = None
LMsk = hlat_s = hlon_s = None
JJ = None
for iplt, rdate in enumerate(RDATES):
  dnmb0 = mtime.rdate2datenum(rdate)
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  print(f"Reading {YR0}/{MM0:02d}/{DD0:02d}")

  irow = iplt // ncol
  icol = iplt % ncol
  ax1 = axes[irow, icol]

  AP_s, hlat_s, hlon_s, LMsk, JJ = read_glorys(rdate, LMsk, hlat_s, hlon_s, JJ)
  if xh is None or yh is None:
    xh, yh = m(hlon_s, hlat_s)

  img = plot_field(ax1, fig1, m, xh, yh, AP_s, clrmp, rmin, rmax, regn, dnmb0)


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

btx = f' @{machine}: plot_GLORYS_ithkn_Nsbplt.py'
#bottom_text(btx, pos=[0.02,0.02], fsz=8)
bottom_text(btx, pos=[0.02, pbtm], fsz=8)


