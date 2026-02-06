"""
  Plot Interpolated CryoSat ice thickness on mesh025

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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6

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

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthice  = os.path.join(pthdata, 'CryoSat2_antarctic_ice_snow_thkn','clim')
fliceout = f'CryoSat_hice_mnthclim_{yrS}_{yrE}_mesh025_{idim}x{jdim}_{regn}.nc'
dfliceout = os.path.join(pthice,fliceout)
print(f'Reading interpolated ice thickness --> {dfliceout}')
with xarray.open_dataset(dfliceout) as ds_hi:
  LON = ds_hi['lon'].data
  LAT = ds_hi['lat'].data
  AI = ds_hi['ice_thkn'].isel(time=MM-1).squeeze()


AI = np.where(np.isnan(AI), 0., AI)
AI[HH >=0] = np.nan

clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 4.
clrmp.set_bad(color=[0.2, 0.2, 0.2])


print("Plotting ...")


m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
parallels = np.arange(-80,-10,10.)
meridians = np.arange(-360,359.,45.)
xh, yh = m(LON,LAT) # GFS coords


plt.ion()
fig1 = plt.figure(1,figsize=(9,9))

plt.clf()
ax1 = plt.axes([0.1, 0.13, 0.8, 0.8])

m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
img1 = ax1.pcolormesh(xh, yh, AI, cmap=clrmp, vmin=rmin, vmax=rmax)
#ax1.contour(xh,yh,HH,[0], linestyles='solid', colors=[(0.,0.,0.)], linewidths=1)
ax1.set_title(f'{fliceout}, ithkn clim MM={MM:02d}\n {pthice}')

# Colorbars
ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img1, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


btx = 'plot_interp_CryoSate_ithkn_mesh025.py'
bottom_text(btx, pos = [0.02,0.02])


