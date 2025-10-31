"""
  CHeck monhtly hsnow fields
  interpolated on mesh025 grid

  see: interp_AMSR_hsnow_antarct_mesh025
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
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)
  
regn = 'south'  # only south region has been done so far

parser = argparse.ArgumentParser()
parser.add_argument("--mm", help="month to plot", type=int, required=True)
args = parser.parse_args()

MM = args.mm if args.mm else None
    
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

# monthly snow depth, Antarctica, original grid:
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthsnow = os.path.join(pthdata,'snow_nasa')
flsnow = 'AMSR_Antarctic_hsnow_month_clim_1998_2007.nc'
dflsnow = os.path.join(pthsnow, flsnow)

print(f"Loading {dflsnow}")
with xarray.open_dataset(dflsnow) as ds_snow:
  LON = ds_snow['lon'].data
  LAT = ds_snow['lat'].data
  HS = ds_snow['snow_depth'].isel(time=MM-1).data.squeeze()

# Interpolated fields:
fliceout = f'SSMI_hsnow_interp_mesh025_1080x1440_mnthly_clim_{regn}.nc'
dfliceout = os.path.join(pthsnow,fliceout)
print(f"Reading interpolated hsnow {dfliceout}")
with xarray.open_dataset(dfliceout) as ds_intrp:
  HSi = ds_intrp['snow_depth'].isel(time=MM-1).data.squeeze()



clrmp = mclrmps.colormap_temp()
rmin = 0.
rmax = 50.

clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]


if regn == 'south':
  xl1 = -8.e6
  xl2 = -1.2e6
  yl1 = xl1
  yl2 = xl2


plt.ion()

m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
xh, yh = m(hlon,hlat) # CICE6 coordinates

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()

# Plot intepolated hsnwo
ax1 = plt.axes([0.05, 0.4, 0.45, 0.45])
m.drawcoastlines()

# draw parallels.
parallels = np.arange(-80,-10,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,359.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

img = ax1.pcolormesh(xh, yh, HSi, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

ax1.set_xlim([xl1, xl2])
ax1.set_ylim([yl1, yl2])
ax1.invert_yaxis()
ax1.invert_xaxis()

ax1.set_title(f"SSM/I hsnow, cm, interp mesh025, M={MM}")


# Plot original field:
LON1 = LON.copy()
LON2 = LON.copy()
LON1 = np.where(LON1 < -175, np.nan, LON1)
LON2 = np.where(LON2 > 172, np.nan, LON2)
LON3 = np.where(LON < 0, LON+360., LON)
LON3 = np.where(LON3 > 350., np.nan, LON3)
lon_cntr1 = [x for x in range(-180,0,45)]  # grey -180:0
lon_cntr2 = [x for x in range(45,178,45)]  # blue: 0 to 180 E
lat_cntr = [x for x in range(-80,-20,10)]

# Plot data on original grid:
ax21 = plt.axes([0.52, 0.4, 0.45, 0.45])
ax21.pcolormesh(HS, cmap=clrmp, vmin=rmin, vmax=rmax)
ax21.axis('scaled')
ax21.set_xlim([23, 290])
ax21.set_ylim([40, 309])

# Check longitudes:
cs = ax21.contour(LON1, lon_cntr1, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
ax21.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
cs2 = ax21.contour(LON2, lon_cntr2, linestyles='solid', colors=[(0.,0.5,0.9)], linewidths=1)
ax21.clabel(cs2, inline=True, fontsize=10, fmt="%.1f")
cs3 = ax21.contour(LON1,[0], linestyles='solid', colors=[(1,0.,0.)], linewidths=2)
ax21.clabel(cs3, inline=True, fontsize=12, fmt="%.1f")
ax21.contour(LON3,[180], linestyles='solid', colors=[(0.6,0.,0.9)], linewidths=1)
ax21.contour(LAT, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
ax21.set_xticks([])
ax21.set_yticks([])

ax21.set_title(f"SSM/I hsnow NASA grid, M={MM}")

ax2 = plt.axes([0.2, 0.3, 0.6, 0.02])
if rmin < 0:
  clb = plt.colorbar(img, cax=ax2, orientation='horizontal', extend='both')
else:
  clb = plt.colorbar(img, cax=ax2, orientation='horizontal', extend='max')

ax2.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_xticklabels(ax2.get_xticks())
ticklabs = clb.ax.get_xticklabels()
#  clb.ax.set_xticklabels(ticklabs,fontsize=10)
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


btx = 'check_interp_hsnow_antarct.py'
bottom_text(btx, pos=[0.1, 0.2])



