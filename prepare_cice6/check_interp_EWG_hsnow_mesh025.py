"""
  Check
  Interpolated snow climatology from EWG snow data
  to MOM6/CICE6 mesh025 grid
  2 options: snow depth, liquid water equivalent

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
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

regn = 'north'

parser = argparse.ArgumentParser()
parser.add_argument("--snfld", help=f"snow field to interpolate",
                    choices=['sndpth','swe'], required=True, type=str)
parser.add_argument("--mm", help="month to plot", required=True, type=int)
args = parser.parse_args()
                             
snfld = args.snfld if args.snfld else None
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

# Get lon/lat for EWG Atlas data
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthnsidc = os.path.join(pthdata,'Warren_snow_clim_EWG_atlas/snow_clim_interp')
pthdump  = os.path.join(pthdata,"gmapi_NSIDC")
fgmapi  = f'EWG_Atlas_MOM6_gmapi_1440x1080_north.nc'
dfgmapi = os.path.join(pthdump, fgmapi)
with xarray.open_dataset(dfgmapi) as dgmapi:  
  LON = dgmapi['longit'].data
  LAT = dgmapi['latit'].data

# Interpolated fields:
if snfld == 'sndpth':
  flsnidc = f'hsnow_EWGatlas_interp_mesh025_{jdm}x{idm}_{regn}.nc'
  varnm = 'snow_depth'
  flsnow = f'snow_depth.{MM}.1954_1991.dat'
else:
  flsnidc = f'swe_EWGatlas_interp_mesh025_{jdm}x{idm}_{regn}.nc'
  varnm = 'swe'
  flsnow = f'swe.{MM}.1954_1991.dat'

print(f"Processing {MM} {varnm} ...")
with xarray.open_dataset(os.path.join(pthnsidc,flsnidc)) as ds_nsidc:
  AI = ds_nsidc[varnm].isel(time=MM-1).data.squeeze()

# Original fields:
pthsnow = os.path.join(pthdata,'Warren_snow_clim_EWG_atlas/DATA/GRIDDED_FIELDS/SNOW_DEPTH')
dflnm = os.path.join(pthsnow,flsnow)
data_hsnow = np.loadtxt(dflnm)
assert data_hsnow.size == 23 * 23, f"Check file length expected {23*23} lines"
AA = data_hsnow.reshape((23, 23))
AA = np.where(AA > 9999., np.nan, AA) * 0.01  # cm --> m  


clrmp = mclrmps.colormap_temp()
rmin = 0.
rmax = 0.4
clrmp.set_bad(color=[0.2, 0.2, 0.2])

LON1 = LON.copy()
LON2 = LON.copy()
LON1 = np.where(LON1 < -175, np.nan, LON1)
LON2 = np.where(LON2 > 172, np.nan, LON2)
LON3 = np.where(LON < 0, LON+360., LON)
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
ax1.contour(LAT, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
ax1.invert_yaxis()
ax1.axis('scaled')
ax1.set_title(f'EWG Atlas snow: {varnm} MM={MM:02d}, m')

m = Basemap(projection='npstere',boundinglat=50,lon_0=-45,resolution='l')
parallels = np.arange(40,89,10.)
meridians = np.arange(-360,359.,45.)

#lons, lats = m.makegrid(idim, jdim) # get lat/lons of ny by nx evenly spaced grid.
xh, yh = m(hlon,hlat) # GFS coords

ax2 = plt.axes([0.55, 0.3, 0.4, 0.4])
# draw parallels.
if regn == 'south':
  parallels = np.arange(-80,-10,10.)
else:
  parallels = np.arange(50,89,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,359.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

img = ax2.pcolormesh(xh, yh, AI, cmap=clrmp, vmin=rmin, vmax=rmax)
ax2.set_title(f'EWG atlas {varnm}  MM={MM:02d} interp to mesh025')

ax3 = fig1.add_axes([0.2, 0.2, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'check_interp_EWG_hsnow_mesh025.py'
bottom_text(btx)


f_ij = False
if f_ij:
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  img1 = ax1.pcolormesh(AI, cmap=clrmp, vmin=rmin, vmax=rmax)







