"""
  Check Interpolated CryoSat ice thickness on mesh025

  gmapi indices:
  get_gmapi_CryoSat_to_mesh025.py

  interpolation:
  interp_CryoSat_ithkn_antarct_mesh025.py

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
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

yrS = 2011
yrE = 2020
MM  = 2
regn = 'south'

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--yrS", help=f"start year of CryoSat clim, default={yrS}", type=int)
parser.add_argument("--yrE", help=f"end year of CryoSat clim, default={yrE}", type=int)
parser.add_argument("--mm", help=f"month of CryoSat data to plot", type=int, required=True)
args = parser.parse_args()
  
regn = args.regn if args.regn else regn
yrS  = args.yrS if args.yrS else yrS
yrE  = args.yrE if args.yrE else yrE
MM   = args.mm if args.mm else None
  
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
    
#hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

jdim, idim = HH.shape

# Read CryoSat clim on original grid:
# Monthly ice thickness, Antarctica, original grid:
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pth0   = os.path.join(pthdata, 'CryoSat2_antarctic_ice_snow_thkn')
pthice = os.path.join(pth0,'clim')
fhice  = f'CryoSat_ithkn_mnthly_clim_316x332_{regn}.nc'
dfhice = os.path.join(pthice, fhice)

print(f'Reading ice thickn climatology {dfhice}')
with xarray.open_dataset(dfhice) as ds_hice:
  LONS = ds_hice['lon'].data
  LATS = ds_hice['lat'].data
  AA = ds_hice['ice_thickness'].isel(time=MM-1).squeeze()

fliceout = f'CryoSat_hice_mnthclim_{yrS}_{yrE}_mesh025_{idim}x{jdim}_{regn}.nc'
dfliceout = os.path.join(pthice,fliceout)
print(f'Reading interpolated ice thickness --> {dfliceout}')
with xarray.open_dataset(dfliceout) as ds_hi:
  LON = ds_hi['lon'].data
  LAT = ds_hi['lat'].data
  AI = ds_hi['ice_thkn'].isel(time=MM-1).squeeze()



clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

LON1 = LONS.copy()
LON2 = LONS.copy()
LON1 = np.where(LON1 < -175, np.nan, LON1)
LON2 = np.where(LON2 > 172, np.nan, LON2)
LON3 = np.where(LONS < 0, LONS+360., LONS)
LON3 = np.where(LON3 > 350., np.nan, LON3)
lon_cntr1 = [x for x in range(-180,0,45)]  # grey -180:0
lon_cntr2 = [x for x in range(45,178,45)]  # blue: 0 to 180 E
lat_cntr = [x for x in range(-80,-20,10)]

print("Plotting ...")

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.05, 0.3, 0.4, 0.4])

ax1.pcolormesh(AA, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.contour(LON1, lon_cntr1, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
ax1.contour(LON2, lon_cntr2, linestyles='solid', colors=[(0.,0.5,0.9)], linewidths=1)
ax1.contour(LON1,[0], linestyles='solid', colors=[(1,0.,0.)], linewidths=1)
ax1.contour(LON3,[180], linestyles='solid', colors=[(0.6,0.,0.9)], linewidths=1)
ax1.contour(LATS, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
ax1.invert_yaxis()
ax1.axis('scaled')
ax1.set_title(f'CryoSat ice thickn clim, MM={MM:02d}')


# Interpolated iconc
xl1 = -8.e6
xl2 = -1.2e6           
yl1 = xl1
yl2 = xl2 

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
#lons, lats = m.makegrid(idim, jdim) # get lat/lons of ny by nx evenly spaced grid.
xh, yh = m(LON,LAT) # GFS coords

ax2 = plt.axes([0.55, 0.3, 0.4, 0.4])
# draw parallels.
parallels = np.arange(-80,-10,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,359.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

img = ax2.pcolormesh(xh, yh, AI, cmap=clrmp, vmin=rmin, vmax=rmax)
ax2.set_xlim([xl1, xl2])
ax2.set_ylim([yl1, yl2])
ax2.invert_yaxis()
ax2.invert_xaxis()

ax2.set_title(f'CryoSate ithkn interp to mesh025, MM={MM:02d}')

ax3 = fig1.add_axes([0.2, 0.2, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'check_interp_CryoSate_ithkn_mesh025.py'
bottom_text(btx)








