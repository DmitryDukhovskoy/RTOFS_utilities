"""
  Plot original ICESat2 ice thickness
  on native grid
  monthly fields

  Summer ice and snow depth, Arctic

  https://zenodo.org/records/18004849

  Monthly gridded summer Arctic sea ice thickness from ICESat-2, v1
  Creators
  Petty, Alek Aaron (Producer)1, 2
  ORCID icon
  Description
  Monthly gridded summer Arctic sea ice thickness from ICESat-2. Produced by combining Release 006 ATL10 freeboards with SnowModel-LG snow loading, processed as in the IS2SITMOGR4 winter Arctic thickness dataset (https://nsidc.org/data/IS2SITMOGR4). Associated paper is currently under review in Journal of Glaciology.

ICESat-2 (Ice, Cloud, and land Elevation Satellite-2) uses a space-based laser altimeter to measure surface elevation changes—especially ice sheets and sea ice—down to centimeter accuracy.

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib  
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
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

regn = 'north'
parser = argparse.ArgumentParser()
parser.add_argument("--yr", help=f"year of ICES data", 
                    choices=[2019,2020,2021],
                    required=True, type=int)
parser.add_argument("--mm", help=f"month of ICES data to plot", 
                    choices=[5,6,7,8],
                    required=True,type=int)
args = parser.parse_args()
  
YR   = args.yr if args.yr else None
MM   = args.mm if args.mm else None
 
if YR == 2021 and MM > 7:
  raise Exception(f"In 2021, last month is July, MM={MM}")

 
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
    
#pthfld, flname = mc6util.pathfname_icesnow_mesh025(fyaml, node_nm, 'ICESat2_orig', YR=YR, MM=MM)
pthfld = os.path.join(pthdata, 'ICESat2_arctic_summer_ithkn_hsnow')
flname = f"IS2SIT_SUMMER_01_{YR}{MM:02d}_006_001.nc"
dflname = os.path.join(pthfld, flname)

# There are many versions of ice thickness estimates in the file
# ice_thickness_sm_e5_int: Monthly mean gridded and smoothed/interpolated sea ice thickness 
#                          calculated using redistributed SnowModel-LG snow loading 
#                          with ERA5 forcing (Liston et al., 2021, 10.5067/27A0P5M6LZBI, SM) 
#                          and fixed ice density (916 kg/m3)
varnm = 'ice_thickness_sm_e5_int'

print(f"Reading {dflname}")

with xarray.open_dataset(dflname) as ds_ices:
  A2d = ds_ices[varnm].values.squeeze()
  LON = ds_ices['longitude'].values
  LAT = ds_ices['latitude'].values
  units = ds_ices[varnm].attrs.get("units", None)

if units == 'cm' or units == 'centimeters':
  cff2m = 100.
elif units == 'm' or units == 'meters':
  cff2m = 1.

A2d = A2d * cff2m


clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 4.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
parallels = np.arange(50, 90, 5)
meridians = np.arange(-360, 359., 45.)

xh, yh = m(LON, LAT) 

print("Plotting ...")

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.13, 0.8, 0.8])

m.drawparallels(parallels, labels=[0,0,0,0])
m.drawmeridians(meridians, labels=[0,0,0,0])
m.drawcoastlines()
img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(f'ICESat2 summer {varnm} {YR}/{MM:02d}')

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_origICESat2_mnth_iconc.py'
bottom_text(btx)


