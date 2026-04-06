"""
  Plot Interpolated NSIDC sea ice concentration 
  to MOM6/CICE6 mesh025 grid

  gmapi indices: get_gmapi_NSIDC_to_mesh025.py

  NSIDC fields from 
  https://noaadata.apps.nsidc.org/NOAA/G02202_V6/north/daily/2025/

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

YR = 2025
MM = 1 
DD = 1

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", type=str, required=True)
parser.add_argument("--yr", help=f"year of NSIDC data, default={YR}", type=int)
parser.add_argument("--mm", help=f"month of NSIDCS data to plot", required=True, type=int)
parser.add_argument("--dd", help="month day to plot", required=True, type=int)
args = parser.parse_args()
  
regn = args.regn if args.regn else None
YR   = args.yr if args.yr else YR
MM   = args.mm if args.mm else MM
DD   = args.dd if args.dd else DD
  
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
    
# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
    
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape


# Interpolated fields:
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")
fliceout = f'NSIDC_iconc_interp_mesh025_{jdm}x{idm}_{YR}{MM:02d}_{regn}.nc'
dfliceout = os.path.join(pthnsidc,fliceout)

print(f'Loading interpolated ice conc {dfliceout}')
with xarray.open_dataset(dfliceout) as dsint:
  AI = dsint['ice_conc'].isel(time=DD-1).data.squeeze()

clrmp = mclrmps.colormap_conc()
rmin = 0.
rmax = 1.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

# Interpolated iconc
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
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.13, 0.8, 0.8])

m.drawparallels(parallels, labels=[0,0,0,0])
m.drawmeridians(meridians, labels=[0,0,0,0])
img = ax1.pcolormesh(xh, yh, AI, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(f'interp NSIDC iconc {YR}/{MM:02d}/{DD:02d}')

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_interp_NSIDC_iconc_mesh025.py'
bottom_text(btx)

