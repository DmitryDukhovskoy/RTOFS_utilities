"""
  Plot interpolated AVHRR ice thickness data

  Interpolate N-day average albedo (varies for summer and winter seasons)
  for south / north regions

    NOAA Climate Data Record (CDR) of AVHRR Polar Pathfinder Extended (APP-X) Cryosphere, Version 2

  Albedo exists only when solar senith angle > 0
  For winter months, most of the polar region = NaN
  and 2 am albedo = NaN for most months except summer

  NOAA Climate Data Record (CDR) of the eXtended AVHRR Polar Pathfinder (APP-X) 
  cryosphere contains 19 geophysical variables over the Arctic and Antarctic for the period 1982 - present. 

  https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C01580


The Polar Pathfinder - Extended Climate Data Record (CDR), utilizes the Advanced Very High Resolution Radiometer (AVHRR) and Visible Infrared Imaging Radiometer Suite (VIIRS) instruments and contains several geophysical variables over the Arctic and Antarctic from 1982–present. The data products are mapped to a 25 km Equal-Area Scalable Earth (EASE) grid at two local solar times: 04:00 and 14:00 for the Arctic, and 02:00 and 14:00 for the Antarctic. 


Cite as: Key, Jeffrey; Wang, Xuanji; Liu, Yinghui; and NOAA CDR Program (2019). NOAA Climate Data Record of AVHRR Polar Pathfinder Extended (APP-X), Version 2. [indicate subset used]. NOAA National Centers for Environmental Information. doi:10.25921/AE96-0E57 [access date].


  gmapi indices:
  get_gmapi_AVHRR_albedo_to_mesh025.py

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

fld_name = 'ithkn' 

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"Region to process",
                    choices=['north','south'], required=True, type=str)
parser.add_argument("--date", help=f"YYYYMMDD to plot", required=True, type=int) 
args = parser.parse_args()

regn = args.regn if args.regn is not None else None
date_plot = args.date if args.date else None

dnmb_req = mtime.dateint2datenum(date_plot)  # requested date numb
YR, MM, DD = mtime.datevec(dnmb_req)[:3]

    
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
jdim, idim = HH.shape
LMsk = np.where(HH<0, 1, 0)

# Mask out not needed latitudes:
if regn == 'south':
  LMsk = np.where(hlat > -55, 0, LMsk)
else:
  LMsk = np.where(hlat < 50, 0, LMsk)

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthinp  = os.path.join(pthdata,'AVHRR_albedo_ithkn')
flice = f"AVHRR_ithkn_{date_plot}_mesh025_1440x1080_{regn}.nc" 
dflice = os.path.join(pthinp, flice)
with xarray.open_dataset(dflice) as dcice:
  A2d = dcice['ice_thkn'].isel(time=0).data.squeeze()

clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)

xh, yh = m(hlon,hlat) # GFS coords

lon_cntrs = [x for x in range(-180,180,45)]
lat_cntrs = [x for x in range(-80,90,10)]

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.82, 0.82])

# Plot original field on grid:
img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
m.drawparallels(parallels, labels=[0,0,0,0])
m.drawmeridians(meridians, labels=[0,0,0,0])

sttl = f'ithkn AVHRR {date_plot}\n {dflice}'
ax1.set_title(sttl, fontsize=10)    



ax3 = fig1.add_axes([0.1, 0.06, 0.8, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_AVHRR_ithkn_mesh025.py'
bottom_text(btx, pos=[0.08,0.02], fsz=8) 



