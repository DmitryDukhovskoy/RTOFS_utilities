"""
  Plot original AMSR2 L4 OSI SAF EUMETSAT sea ice concentration 
  on native grid
  daily fields

    string :title = "Level-3 Sea Ice Concentration Analysis (AMSR2) from OSI SAF EUMETSAT" ;
    string :institution = "EUMETSAT OSI SAF" ;
    string :history = "Created 2025-05-13 02:01:07" ;
    string :contact = "osisaf-manager@met.no" ;
    string :references = "Product User Manual and Algorithm Theoretical Basis Document available at https://osi-saf.eumetsat.int/documentation/products-documentation" ;
    string :subset\:source = "ARCO data downloaded from the Marine Data Store using the MyOcean Data Portal" ;

Downloaded from:
https://data.marine.copernicus.eu/product/SEAICE_ARC_PHY_AUTO_L4_MYNRT_011_024/services

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib  
import xarray as xr
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
import mod_misc1 as mmisc

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", type=str, required=True)
parser.add_argument("--rdate", help=f"date to plot YYYYMMDD", required=True,type=int)
args = parser.parse_args()
  
regn  = args.regn if args.regn else None
rdate = args.rdate
  
dnmb0 = mtime.rdate2datenum(rdate)
YR, MM, DD = mtime.datevec(dnmb0)[:3]

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

print(f"Processing {YR}/{MM}/{DD} {regn} ...")
pthamsr = f'/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/AMSR2_OSISAF_L4/{regn}/{YR}'    
flamsr = f"osisaf_{regn}_nrt_amsr2_L4_{MM:02d}{YR}.nc"
dfl = os.path.join(pthamsr,flamsr)
assert os.path.isfile(dfl), f"Missing data: {dfl}"

with xr.open_dataset(dfl) as ds_amsr:
  AA = ds_amsr['ice_conc'].isel(time=DD-1).values.squeeze() * 0.01  # % to fractions
  LON = ds_amsr['longitude'].values
  LAT = ds_amsr['latitude'].values  

hlon, hlat = np.meshgrid(LON, LAT)


clrmp = mclrmps.colormap_conc()
rmin = 0.
rmax = 1.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

# Interpolated iconc
if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)

xh, yh = m(hlon,hlat) 

print("Plotting ...")

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.13, 0.8, 0.8])

m.drawparallels(parallels, labels=[0,0,0,0])
m.drawmeridians(meridians, labels=[0,0,0,0])
img = ax1.pcolormesh(xh, yh, AA, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(f"AMSR2 OSI SAF iconc {YR}/{MM:02d}/{DD:02d}")

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_origAMSR2_daily_iconc.py'
bottom_text(btx)

