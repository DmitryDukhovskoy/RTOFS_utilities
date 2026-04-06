"""
  Plot snow depth - daily clim deriver from
  NASA SSM/I daily fields 
  interpolated to mesh025 

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
from copy import copy
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
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob

init_date = 20250103
regn = 'south'

# hs_h - grid cell mean (!) snow thickness, m
# snow_ai - snowfall rate cm/day (in liquid water equivalent !)
# dsnow_h - snow formation (cm/day) - can be > or < 0
# snoice_h - snow-ice formation (cm/day)
# melts_h  - top snow melt (cm/day)
parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--mm", help="Month to plot", type=int, required=True)
parser.add_argument("--dd", help="day month to plot", type=int, required=True)
parser.add_argument("--punit", help="plot units default=m", choices=['cm','m'], type=str)
args = parser.parse_args()

regn = args.regn if args.regn else regn
MM   = args.mm if args.mm else None
DD   = args.dd if args.dd else None
plot_units = args.punit if args.punit else 'm'

syst_info = os.uname()
machine = syst_info.nodename

if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
else:
  print("Unknown machine:", machine)

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

# Read interpolated snow depths:
# Snow depth climatology, Interpolated fields mesh025:
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthsnow = os.path.join(pthdata,'snow_nasa','daily_clim')
flhsn = f"SSMI_hsnow_mesh025_1440x1080_dailyclim_{MM:02d}_south.nc" 
dflhsn = os.path.join(pthsnow,flhsn)
print(f"Reading interpolated hsnow {dflhsn}")
with xarray.open_dataset(dflhsn) as ds_snow:
  units = ds_snow['snow_depth'].attrs.get('units')
  if units == 'm':
    if plot_units == 'm':
      cff = 1.
    else:
      cff = 0.01
  elif units == 'cm':
    if plot_units == 'm':
      cff = 100.
    else:
      cff = 1.
  else:
    raise Exception(f"input units are not recognized: {units}")

  HSi = cff * ds_snow['snow_depth'].isel(time=DD-1).data.squeeze()
  LON = ds_snow['lon'].data
  LAT = ds_snow['lat'].data


units = plot_units
#clrmp = mclrmps.colormap_uv()
clrmp = mclrmps.colormap_temp()
if plot_units == 'cm':
  rmin = 0
  rmax = 40
else:
  rmin = 0
  rmax = 0.4

clrmp.set_bad(color=[0.2, 0.2, 0.2])

sttl = f'hsnow SSM/I daily clim {MM:02d}/{DD:02d}'
sinfo = f"Daily climatology snow depth over ice {regn}, from SSM/I daily fields"

plt.ion()

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
  xl1 = -8.5e6
  xl2 = -1.e6
  yl1 = xl1
  yl2 = xl2

xh, yh = m(LON, LAT) # CICE6 coordinates

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])
m.drawcoastlines()

# draw parallels.
if regn == 'south':
  parallels = np.arange(-80,-10,10.)
else:
  parallels = np.arange(40, 89, 10)

m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,359.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

img = ax1.pcolormesh(xh, yh, HSi, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

if regn == 'south':
  ax1.set_xlim([xl1, xl2])
  ax1.set_ylim([yl1, yl2])
  ax1.invert_yaxis()
  ax1.invert_xaxis()

ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
if rmin < 0:
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
else:
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')

ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.03, 0.8, 0.06])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')

btx = 'map_hsnow_dailyclim_SSMI.py'
bottom_text(btx, pos=[0.2, 0.01])


