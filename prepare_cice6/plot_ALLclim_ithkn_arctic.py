"""
  12 suplots of monthly ice thickness clim in Arctic

  Combined climatologies: 
    AVHRR monthly clim (2016-2025)
    combined CryoSat AWI mnthly clim (2015-2024 Jan-Apr, Oct-Dec) + NSIDC EASE 100 for summer months (May - Sept)
    ICESat-2: laser monthly summer (May-Aug 2019-2020/2021)

   derived:
   derive_mnthclimALL_ithkn_arctic_mesh025.py

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib  
import xarray
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
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

regn = 'north'
fld_name = 'ithkn'
min_hi = 0.6      # for plotting, mask low thickness

print(f"Min ice thickness threhold = {min_hi:.3f}")

parser = argparse.ArgumentParser()
parser.add_argument("--mm", help="Month to plot", required=True, type=int)
args = parser.parse_args()
  
field_name = 'ithkn'
MM = args.mm

syst_info = os.uname() 
machine = syst_info.nodename
  
if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea node:", machine)
  node_nm = "gaea"
elif 'an' in machine:
  print("Running on PPAN node:", machine)
  node_nm = "ppan"
else:
  print("Unknown machine:", machine)
    
fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)


pthdata   = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthice    = os.path.join(pthdata, 'ithkn_clim_combined')
fliceout  = 'ithkn_mnthclim_cryo_avhrr_ices_1440x1080_north.nc'
dfliceout = os.path.join(pthice,fliceout)

with xarray.open_dataset(dfliceout) as dsice:
  hlon = dsice['lon'].data
  hlat = dsice['lat'].data


clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.2, 0.2, 0.2])
clrmp.set_under(color=[1,1,1])

clrmp.set_bad(color=[0.2, 0.2, 0.2])
clrmp.set_under(color=[1,1,1])

from mpl_toolkits.basemap import Basemap, cm
m = Basemap(projection='npstere', boundinglat=60, lon_0=-10,resolution='l')
#m = Basemap(projection='npstere', boundinglat=53, lon_0=-10,resolution='l')
xh, yh = m(hlon,hlat) # GFS coords

def plot_field(ax1, fig1, m, xR, yR, A2d, clrmp, rmin, rmax, plt_clrb, sttl=[]):
  fig1.sca(ax1)
  m.drawcoastlines()
  parallels = np.arange(40,89,10.)
  meridians = np.arange(-360,359.,45.)

  img = ax1.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  m.drawparallels(parallels,labels=[0,0,0,0])
  m.drawmeridians(meridians,labels=[0,0,0,0])
  ax1.set_title(sttl)

  # extend: min, max, both
  if plt_clrb:
    ax2 = fig1.add_axes([0.92,0.1,0.013,0.8])
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')

    ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
    ax2.set_yticklabels(ax2.get_yticks())
    ticklabs = clb.ax.get_yticklabels()
    #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
    clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
    clb.ax.tick_params(direction='in', length=12)

  return ax1

plt.ion()
fig1 = plt.figure(1,figsize=(9, 9))
fig1.clf()  # Clear the figure
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

imo = MM - 1

print(f"Plotting {MM:02d}")

with xarray.open_dataset(dfliceout) as dsice:
  A2d = dsice['ice_thkn'].isel(time=imo).data.squeeze()

if min_hi > 0. and min_hi is not None:
  A2d[A2d <= min_hi] = 0.


sttl = f'ithkn fused clim {MM:02d}'

plt_clrb = True
ax1 = plot_field(ax1, fig1, m, xh, yh, A2d, clrmp,rmin,rmax,plt_clrb,sttl=sttl)

sinfo = 'combined clim (AVHRR+CryoSat AWI + NSIDC EASE100 + ICESat-2)'
bottom_text(sinfo, pos=[0.05,0.08], fsz=10, ipwd=0)

btx = 'plot_ALLclim_ithkn_arctic.py'
bottom_text(btx, pos=[0.05,0.05], fsz=10)


