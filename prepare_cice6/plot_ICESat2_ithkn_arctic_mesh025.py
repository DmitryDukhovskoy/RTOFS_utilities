"""
  Plot Interpolated CryoSatice thickn. monthly fileds 
  summer months only

  When averaging, 0-thickness is ignored, i.e. only hice > 0 is considered for 
  computing the multi-year mean to track mean ice thickness in the grid cell 
  only when it presents there


  Summer ice and snow depth, densities, Arctic

  https://zenodo.org/records/18004849
  
  Monthly gridded summer Arctic sea ice thickness from ICESat-2, v1
  Creators
  Petty, Alek Aaron (Producer)1, 2
  ORCID icon
  Description
  Monthly gridded summer Arctic sea ice thickness from ICESat-2. 
  Produced by combining Release 006 ATL10 freeboards with SnowModel-LG snow loading, 
  processed as in the IS2SITMOGR4 winter Arctic thickness dataset (https://nsidc.org/data/IS2SITMOGR4). 
    
ICESat-2 (Ice, Cloud, and land Elevation Satellite-2) uses a space-based laser altimeter to measure surface elevation changes—especially ice sheets and sea ice—down to centimeter accuracy.

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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6

regn = 'north'
field_name = 'ithkn' 

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="Year to plot", choices=[2019,2020,2021], required=True, type=int)
parser.add_argument("--mm", help="Month to plot", choices=[5,6,7,8], required=True, type=int)
args = parser.parse_args()
 
YR = args.yr
MM = args.mm
 
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
pthintrp = os.path.join(pthdata,'ICESat2_arctic_summer_ithkn_hsnow','interp_mesh025')
fliceout = f"{field_name}_ICESat2_arctic_{YR}{MM:02d}_{jdm}x{idm}.nc"
dfliceout = os.path.join(pthintrp, fliceout)
if not os.path.isfile(dfliceout):
  raise RuntimeError(f"File not found: {dfliceout}")

with xarray.open_dataset(dfliceout) as dsice:
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
ax1.set_title(f"ithkn ICESat-2 laser altim. {YR}/{MM:02d}, interp mesh025")

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ticks = np.linspace(rmin, rmax, 11)
clb.set_ticks(ticks)
clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_ICESat2_ithkn_arctic_mesh025.py'
bottom_text(btx, pos=[0.02,0.02], fsz=8)


