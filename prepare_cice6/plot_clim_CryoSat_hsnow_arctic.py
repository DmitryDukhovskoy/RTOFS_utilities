"""
  Plot 1 month from monthly hsnow clim Arctic

  Climatology derived from 
  Interpolated NSIDC CryoSat snow monthly fileds 2018-2021
  winter months only

  Warren (EWG Atlas) snow depth climatology - for summer months

  Both data sets have been interpolated onto 025 mesh

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
import pandas as pd
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
field_name = 'sndpth'

parser = argparse.ArgumentParser()
parser.add_argument("--mm", help="Month to plot", 
                    choices=[1,2,3,4,5,6,7,8,9,10,11,12], required=True, type=int)
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

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
    
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
LMsk = np.where(HH<0, 1, 0)

fliceout = 'CryoSat_EWG_hsnow_mnthclim_mesh025_1440x1080_north.nc'
varnm = 'snow_depth'

pthclim = os.path.join(pthdata,'CryoSat_arctic_ice_snow_thkn/clim')
dfliceout = os.path.join(pthclim,fliceout)
print(f'Processing climatology fields  --> {dfliceout}')
with xarray.open_dataset(dfliceout) as dsn:
  units = dsn['snow_depth'].attrs.get("units", None)
  if units == 'cm' or units == 'centimeter':
    cff = 0.01   # cm --> m
  elif units == 'm' or units == 'meter':
    cff = 1.
  A2d = cff * dsn[varnm].isel(time=MM-1).data.squeeze()

clrmp = mclrmps.colormap_temp()
rmin = 0.
rmax = 0.4
clrmp.set_bad(color=[0.2, 0.2, 0.2])
clrmp.set_under(color=[1,1,1])

from mpl_toolkits.basemap import Basemap, cm
m = Basemap(projection='npstere', boundinglat=60, lon_0=-10,resolution='l')
xh, yh = m(hlon,hlat) # GFS coords

plt.ion()
fig1 = plt.figure(1,figsize=(9, 9))
fig1.clf()  # Clear the figure

ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
parallels = np.arange(40,89,10.)
meridians = np.arange(-360,359.,45.)

img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
m.drawparallels(parallels,labels=[0,0,0,0])
m.drawmeridians(meridians,labels=[0,0,0,0])

sttl = f"CryoSat hsnow clim, m, interp mesh025, M={MM}\n"
sttl = sttl + f"{pthclim}\n"
sttl = sttl + f'{fliceout}'

ax1.set_title(sttl, fontsize=10)

# extend: min, max, both
ax2 = plt.axes([0.2, 0.06, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax2, orientation='horizontal', extend='both')
#clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')

ax2.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_xticklabels(ax2.get_xticks())
ticklabs = clb.ax.get_xticklabels()
#  clb.ax.set_xticklabels(ticklabs,fontsize=10)
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_clim_CryoSat_hsnow_arctic.py'
bottom_text(btx, pos=[0.05,0.02], fsz=10)


