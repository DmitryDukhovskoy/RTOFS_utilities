"""
  Check original fields of monthly ice thicknesses
  from CPOM ESA 

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

YR  = 2021
MM  = 10
regn = 'north'

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help=f"year of CryoSat data, default={YR}", type=int)
parser.add_argument("--mm", help=f"month of CryoSat data to plot, default={MM}", type=int)
args = parser.parse_args()
  
YR  = args.yr if args.yr else YR
MM  = args.mm if args.mm else MM
  
# Read CryoSat clim on original grid:
# Monthly ice thickness, Antarctica, original grid:
#pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdata = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/'
pthice  = os.path.join(pthdata, 'CryoSat_CPOM_ESA_arctic_ithkn')
fhice  = f'thk_{YR}_{MM}.map.nc'
dfhice = os.path.join(pthice, fhice)

print(f'Reading ice thickn climatology {dfhice}')
with xarray.open_dataset(dfhice) as ds_hice:
  LONS = ds_hice['longitude'].data
  LATS = ds_hice['latitude'].data
  AA = ds_hice['thickness'].data


clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

print("Plotting ...")

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.15, 0.8, 0.8])

m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l', ax=ax1)
xh, yh = m(LONS,LATS) 

m.drawparallels(np.arange(60, 90, 5), labels=[1,0,0,0])
m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,1])
m.drawcoastlines()

sc = ax1.scatter(xh, yh, c=AA, s=10, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(f"ithkn, CPOM ESA, {YR}/{MM:02d}")

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(sc, cax=ax3, orientation='horizontal', extend='max')
ticks = np.linspace(rmin, rmax, 11)
clb.set_ticks(ticks)
clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_ithk_CryoSat_CPOM_ESA.py'
bottom_text(btx, pos=[0.02,0.02], fsz=8)








