"""
  Plot interpolated VIIRS ice thickness to mesh025

  ice thickness Arctic
  daily fields, 750-m resolution

  1 or several days
  using OpeNDaP or from downloaded files (not recommended - too big)

https://coastwatch.noaa.gov/cwn/products/viirs-sea-ice-concentration-ice-thickness-ice-surface-temperature.html

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
import mod_time as mtime
import mod_colormaps as mclrmps
import mod_mom6 as mmom6

regn  = 'north'

parser = argparse.ArgumentParser()
parser.add_argument("--idate", help="Date to plot YYYYMMDD, 2021-2025, MM=9", required=True, type=int)
args = parser.parse_args()

date_plot = args.idate
dnmbI = mtime.rdate2datenum(date_plot*100)  # init. day nmb
YR, MM, DD = mtime.datevec(dnmbI)[:3]

assert MM == 9, f"VIIRS daily fields are only for Sept. requested MM={MM}"
assert YR > 2020 and YR < 2026, f"VIIRS data for 2021-2025, requested YR={YR}"
  
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

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

# Convert to negative depths:
if np.nanmin(HH) > -1e-6:
  HH = np.where(HH < 1.e-6, np.nan, HH) # assuming land ~0
  HH = -HH
  HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape
LMsk = np.where(HH<0, 1, 0)

# Mask out not needed latitudes:
LMsk = np.where(hlat < 50, 0, LMsk)

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthviirs = os.path.join(pthdata,'VIIRS_ithkn_highres','interp_daily')
flnm = f"ithkn_VIIRS_arctic_{YR}{MM:02d}{DD:02d}_1080x1440.nc"
dflnm = os.path.join(pthviirs, flnm)

print(f"Reading {dflnm}")
with xarray.open_dataset(dflnm) as dsice:
  A2d = dsice['ice_thkn'].isel(time=0).data.squeeze()

clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.2, 0.2, 0.2])
clrmp.set_under(color=[1,1,1]) 

AP = A2d.squeeze()
AP[HH >= 0] = np.nan   # land
AP[np.isnan(AP) & (HH < 0)] = -1.  # ocean

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
ax1.set_title(f"ithkn VIIRS {YR}/{MM:02d}/{DD:02d}, interp mesh025")

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ticks = np.linspace(rmin, rmax, 11)
clb.set_ticks(ticks)
clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_VIIRS_arctic_ithkn_mesh025.py'
bottom_text(btx, pos=[0.02,0.02], fsz=8)





