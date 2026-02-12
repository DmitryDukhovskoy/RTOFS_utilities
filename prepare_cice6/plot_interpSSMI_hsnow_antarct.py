"""
  Plot SSMI monhtly hsnow fields
  interpolated on mesh025 grid

  see: interp_SSMI_hsnow_monthly_antarct_mesh025.py
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
  
regn = 'south'  # only south region has been done so far
yrS = 1992
yrE = 2007

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

# Interpolated fields
# monthly snow depth, Antarctica
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthsnow = os.path.join(pthdata,'snow_nasa','monthly_clim')
fliceout = f'SSMI_hsnow_mnthclim_{yrS}_{yrE}_mesh025_1440x1080_{regn}.nc'
dfliceout = os.path.join(pthsnow,fliceout)

print(f"Reading interpolated hsnow {dfliceout}")
cff = 1.
with xarray.open_dataset(dfliceout) as ds_intrp:
  HSi = ds_intrp['snow_depth'].isel(time=MM-1).data.squeeze()
  units = ds_intrp['snow_depth'].attrs.get("units", None)

  if units == 'cm' or units == 'centimeter':
    HSi = 0.01 * HSi   # cm --> m


clrmp = mclrmps.colormap_temp()
rmin = 0.
rmax = 0.4

clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]


print('Plotting ...')
plt.ion()

m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
parallels = np.arange(-80,-10,10.)
meridians = np.arange(-360,359.,45.)

xh, yh = m(hlon,hlat) # CICE6 coordinates

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()

# Plot intepolated hsnwo
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()

m.drawparallels(parallels,labels=[0,0,0,0])
m.drawmeridians(meridians,labels=[0,0,0,0])

img = ax1.pcolormesh(xh, yh, HSi, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

sttl = f"SSM/I hsnow, m, interp mesh025, M={MM}\n" 
sttl = sttl + f"{pthsnow}\n"
sttl = sttl + f'{fliceout}'
ax1.set_title(sttl, fontsize=10)

ax2 = plt.axes([0.2, 0.06, 0.6, 0.02])
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


btx = 'plot_interpSSMI_hsnow_antarct.py'
bottom_text(btx, pos=[0.02, 0.02])



