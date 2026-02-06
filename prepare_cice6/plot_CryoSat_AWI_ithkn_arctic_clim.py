"""
  plot 1 month of  winter monthly snow depth clim in Arctic

  Interpolated CryoSat AWI  snow or ice thickn. monthly fileds 
  winter months only
  see:
  interp_CryoSat_arcticAWI_iceflds_mesh025.py

  Summer months:
  Warren (EWG Atlas) snow depth climatology
  interp_EWG_snow_mesh025.py

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
import mod_cice6_utils as mc6util

regn = 'north'
field = 'ithkn'

parser = argparse.ArgumentParser()
parser.add_argument("--field", help=f"Field to interpolate: snow depth or ice thickn",
                    choices=['sndpth','ithkn'], type=str)
parser.add_argument("--mm", help=f"Winter month to plot", choices=[1,2,3,4,10,11,12],
                   required=True, type=int)
args = parser.parse_args()
  
field_name = args.field if args.field else field
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

#pthclim = os.path.join(pthdata,'CryoSat_AWI_arctic_ithkn/clim')
if field_name == 'ithkn':
  varnm = 'ice_thkn'
  pthclim, fliceout = mc6util.pathfname_icesnow_mesh025(fyaml, node_nm, 'ithkn_AWI')
dfliceout = os.path.join(pthclim,fliceout)

print(f'Processing climatology fields  --> {dfliceout}')
with xarray.open_dataset(dfliceout) as dsn:
  A3d = dsn[varnm].data
  mnth_saved = dsn['time'].values
  irec = np.where(mnth_saved == MM)[0]

  A2d = A3d[irec[0],:].squeeze()
  A2d[HH >= 0] = np.nan   # land
  A2d[np.isnan(A2d) & (HH < 0)] = -1.  # ocean

clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.2, 0.2, 0.2])
clrmp.set_under(color=[1,1,1])

from mpl_toolkits.basemap import Basemap, cm
#m = Basemap(projection='npstere', boundinglat=60, lon_0=-45,resolution='l')
m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
xh, yh = m(hlon,hlat) # GFS coords

plt.ion()
fig1 = plt.figure(1,figsize=(9, 9))
fig1.clf()  # Clear the figure

ax1 = plt.axes([0.05, 0.1, 0.8, 0.8])

m.drawcoastlines()
parallels = np.arange(40,89,10.)
meridians = np.arange(-360,359.,45.)

img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
m.drawparallels(parallels,labels=[0,0,0,0])
m.drawmeridians(meridians,labels=[0,0,0,0])

sttl = f"{field_name} AWI monthly clim MM={MM:02d}"
ax1.set_title(sttl)

ax2 = fig1.add_axes([0.9,0.1,0.015,0.8])
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')

ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=10)


btx = 'plot_CryoSat_AWI_ithkn_arctic_clim.py'
bottom_text(btx, pos=[0.05,0.05], fsz=10)


