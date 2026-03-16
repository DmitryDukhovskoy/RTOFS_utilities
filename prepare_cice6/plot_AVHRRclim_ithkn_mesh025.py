"""
  Plot monthly clim AVHRR fields 
  see: derive_mnthclim_AVHRR_ithkn_mesh025


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

fld_name = 'ithkn' 
YRS = 2016
YRE = 2025

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"Region to process",
                    choices=['north','south'], required=True, type=str)
parser.add_argument("--mm", help="Month to plot", required=True, type=int) 
args = parser.parse_args()

regn = args.regn if args.regn is not None else None
MM = args.mm if args.mm else None

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

pthdata   = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthice    = os.path.join(pthdata, 'AVHRR_albedo_ithkn','tmp')
fliceout  = f'AVHRR_{fld_name}_mnthclim_{YRS}-{YRE}_1440x1080_{regn}.nc'
dfliceout = os.path.join(pthice,fliceout)

with xarray.open_dataset(dfliceout) as dsice:
  A2d = dsice['ice_thkn'].isel(time=MM-1).data.squeeze()
  hlon = dsice['lon'].data
  hlat = dsice['lat'].data


clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)

xh, yh = m(hlon,hlat) # GFS coords

lon_cntrs = [x for x in range(-180,180,45)]
lat_cntrs = [x for x in range(-80,90,10)]

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.82, 0.82])

# Plot original field on grid:
img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
m.drawparallels(parallels, labels=[0,0,0,0])
m.drawmeridians(meridians, labels=[0,0,0,0])

sttl = f'ithkn AVHRR clim MM={MM:02d}\n {dfliceout}'
ax1.set_title(sttl, fontsize=10)    


ax3 = fig1.add_axes([0.1, 0.06, 0.8, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_AVHRRclim_ithkn_mesh025.py'
bottom_text(btx, pos=[0.08,0.02], fsz=8) 



