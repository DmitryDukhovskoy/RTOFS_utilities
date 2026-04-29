"""
  Plot interpolated RTOFS CICE4 onto mesh025 grid
  gmapi:
  get_gmapi_RTOFS_to_mesh025.py

  interp RTOFS:
  interp_RTOFS_CICE4_iconc_mesh025.py

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

init_hr = 0

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", required=True, type=str)
parser.add_argument("--init", help="RTOFS init date", 
                   choices=[20250704, 20251231], required=True, type=int)
parser.add_argument("--field", help="Field to plot",
                    choices=['iconc','ithkn','hsnow'], required=True, type=str)
args = parser.parse_args()
  
regn = args.regn 
init_date = args.init
fldnm = args.field
  
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
LON, LAT  = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape

# Read RTOFS - CICE4
# Init date:
dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
yrI, mmI, ddI, hrI = mtime.datevec(dnmbI)[:4]

if dnmbI < mtime.datenum([2025,8,1]):
  rtofs_vers = "2.4"
else:
  rtofs_vers = "2.5" 

pthice = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/cice4_fcast/{init_date}/interp_mesh025"
flice  = f"{fldnm}_RTOFSv{rtofs_vers}_{init_date}_{jdm}x{idm}.nc"
if fldnm == 'ithkn':
  varnm  = "ice_thkn"
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
elif fldnm == 'iconc':
  varnm = "ice_conc"
  clrmp = mclrmps.colormap_conc()
  rmin = 0.
  rmax = 1.
elif fldnm == 'hsnow':
  varnm = 'snow_depth'
  clrmp = mclrmps.colormap_temp()
  rmin = 0.
  rmax = 0.4
clrmp.set_bad(color=[0.2, 0.2, 0.2])

dfliceout = os.path.join(pthice,flice)
print(f'Reading interpolated ice thickness --> {dfliceout}')
with xarray.open_dataset(dfliceout) as ds_hi:
  AI = ds_hi[varnm].data.squeeze()


AI = np.where(np.isnan(AI), 0., AI)
AI[HH >=0] = np.nan


print("Plotting ...")


if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)
xh, yh = m(LON,LAT) # GFS coords


plt.ion()
fig1 = plt.figure(1,figsize=(9,9))

plt.clf()
ax1 = plt.axes([0.1, 0.13, 0.8, 0.8])

m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
img1 = ax1.pcolormesh(xh, yh, AI, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(f'{fldnm} RTOFS interp mesh025, {flice}\n {pthice}')

# Colorbars
ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img1, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


btx = 'plot_interp_RTOFS_ice_snow_mesh025.py'
bottom_text(btx, pos = [0.02,0.02])


