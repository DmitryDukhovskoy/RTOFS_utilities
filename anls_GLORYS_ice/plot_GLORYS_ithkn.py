"""
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
parser.add_argument("--rdate", help="Plot date: YYYYMMDD, years=1993, ...", required=True, type=int)
args = parser.parse_args()

rdate = args.rdate

dnmbR = mtime.rdate2datenum(rdate*100)  # restart day nmb
YR, MM, DD = mtime.datevec(dnmbR, round_hrs=True)[:3]


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


pthice = f"/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick/{YR}"

# Find file:
from pathlib import Path
try:
  dflice = next(
      Path(pthice).glob(
          f"*_mean_{rdate}_R*.nc"
      )
  )
  print(f"Found file: {dflice}")
except StopIteration:
  print(f"No file found for {rdate} in ${pthice}")


#dflice = os.path.join(pthice, flice)
if not os.path.isfile(dflice):
  raise RuntimeError(f"File not found: {dfliceout}")

with xarray.open_dataset(dflice) as dsice:
  A2d = dsice['sithick'].isel(time=0).data.squeeze()
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values 

hlon, hlat = np.meshgrid(LON, LAT)

clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.2, 0.2, 0.2])
clrmp.set_under(color=[1,1,1])

AP = A2d.squeeze()
#AP[HH >= 0] = np.nan   # land
#AP[np.isnan(AP) & (HH < 0)] = -1.  # ocean

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.15, 0.8, 0.8])

m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l', ax=ax1)
xh, yh = m(hlon, hlat)

m.drawparallels(np.arange(60, 90, 5), labels=[0,0,0,0])
m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,0])
m.drawcoastlines()

img = m.pcolormesh(xh,yh, AP, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(f"ithkn GLORYS {YR}/{MM:02d}/{DD:02d}")

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ticks = np.linspace(rmin, rmax, 11)
clb.set_ticks(ticks)
clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
clb.ax.tick_params(direction='in', length=12)

btx = f' @{machine}: plot_GLORYS_ithkn.py'
bottom_text(btx, pos=[0.02,0.02], fsz=8)


