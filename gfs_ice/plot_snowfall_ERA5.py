"""
  hourly snowfall from ECMWF ERA5 reanalysis
  https://cds.climate.copernicus.eu/requests?tab=all
  extracted for Southern Ocean region

  Note snowfall is given in "m of water equivalent"
  lwe_thickness_of_snowfall_amount

  m(water) --> kg / m2 *sec = fsnow_m (m) * rho_water (kg/m3) / dt (N secs)

  snow kg/(m2*sec) --> cm / day = fsnow (kg/ m2*sec) *1/rho_snow (m3/kg) * (24*3600) * 100)

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
from mpl_toolkits.basemap import Basemap, cm
import argparse
import pandas as pd

PPTHN = '/home/Dmitry.Dukhovskoy/python'
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob

dayS = 20250110
dayE = 20250111
rho_water = 1000.
rho_snow = 300.  # GFS atmos, snow rho = 200, CICE=300
pltfld = 'mean'  # mean snowfall rate or cumulative fields to plot

# hs_h - grid cell mean (!) snow thickness, m
# snow_ai - snowfall rate cm/day
# dsnow_h - snow formation (cm/day) - can be > or < 0
# snoice_h - snow-ice formation (cm/day)
# melts_h  - top snow melt (cm/day)
parser = argparse.ArgumentParser()
parser.add_argument("--days", help="for avrg start date=20250101", required=True, type=int)
parser.add_argument("--daye", help="for avrg start date=20250110", type=int)
parser.add_argument("--pfld", help=f"Plot mean or accumulated field, default={pltfld}", 
                    choices=['mean','cumul'], 
                    type=str)
args = parser.parse_args()

dayS = args.days if args.days else None
dayE = args.daye if args.daye else dayS
pltfld = args.pfld if args.pfld else pltfld

pthoutp = '/work/Dmitry.Dukhovskoy/data/ERA5_Antarctica'
flinp = 'ecmwf_snowfall_2024_2025_Antarctica.nc'
dflice = os.path.join(pthoutp,flinp)


#  Get output time:
print(f"Reading Time from {dflice}")
with xarray.open_dataset(dflice) as ds:
  time_nc = ds['valid_time'].data  # datetime64 
  lon1d = ds['longitude'].data
  lat1d = ds['latitude'].data

TLON, TLAT = np.meshgrid(lon1d,lat1d)

# Output freq., sec
dTsec = (time_nc[1] - time_nc[0]) / np.timedelta64(1,'s')

TM = []
for irc in range(len(time_nc)):
  t0 = pd.Timestamp(time_nc[irc])
  year = t0.year
  month = t0.month
  day = t0.day
  hour = t0.hour

  dnmb = mtime.datenum([year,month,day,hour])
  TM.append(dnmb)

TM = np.array(TM)

# Start - end date numb:
dnmbS = mtime.dateint2datenum(dayS)
dnmbE = mtime.dateint2datenum(dayE)

# Find start and end of averaging period:
iS = np.argmin(abs(TM-dnmbS))
assert TM[iS]-dnmbS < dTsec/86400., f"index={iS} check start date and TM"

iE = np.argmin(abs(TM-dnmbE))
assert TM[iE]-dnmbE < dTsec/86400., f"indx={iE} check end date and TM"


irec = 0
AAsum = None
for indx in range(iS,iE+1):
  dnmb0 = TM[indx]
  if dnmb0 < dnmbS:
    continue
  elif dnmb0 >= dnmbE+1:
    # Process dnmbE and quit the next day
    break

  print(f"Reading {mtime.datestr(dnmb0)}")
  with xarray.open_dataset(dflice) as ds:
    A2d = ds['sf'].isel(valid_time=indx).data.squeeze()  # m of water equivalent

  # Convert to cm/day - similar to CICE output and average
  Fsnow = A2d *rho_water / dTsec   # kg / m2 * sec
  Fsnow = (Fsnow / rho_snow) * dTsec * 100    # cm of snow / (time step sec.) - per 1 hr

  if AAsum is None:
    AAsum = Fsnow.copy()
  else:
    AAsum = AAsum + Fsnow

  irec += 1

# AAsum - cumulative thickness of snowfall cm  / irec --> cm/hr (mean) * 24 hrs --> cm/day
# Convert AAsum --> cm/day (mean), 
A2d = AAsum * 24./float(irec)   # cumul. cm of hourly snowfall --> mean cm/day

if pltfld == 'cumul':
  A2d = AAsum.copy()    # cumul snowfall in cm over all hourly output


if pltfld == 'mean':
  varnm = 'snowfall'
  units = 'cm/day'
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 1.
elif pltfld == 'cumul':
  varnm = 'cumul snowfall'
  units = 'cm'
  clrmp = mclrmps.colormap_warm()
  rmin = 0.
  rmax = 10.
  
clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.6,0.6,0.6]

sinfo = 'ECMWF ERA5 Reanalysis NRT, snowfall for Southern Ocean\n'
sinfo = sinfo + f"Converted from m of eq. water to cm/day snow using rho_snow={rho_snow:.1f}\n"
sinfo = sinfo + dflice

DV1 = mtime.datevec(dnmbS)
DV2 = mtime.datevec(dnmbE)
sttl = f'ERA5 {varnm} ({units}), avrg: {DV1[0]}/{DV1[1]}/{DV1[2]} - {DV2[0]}/{DV2[1]}/{DV2[2]}'

plt.ion()

m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
#lons, lats = m.makegrid(idim, jdim) # get lat/lons of ny by nx evenly spaced grid.
#x, y = m(lons, lats) # compute map proj coordinates.
xh, yh = m(TLON,TLAT) # GFS coords

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])
m.drawcoastlines()

# draw parallels.
parallels = np.arange(-80,-10,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,359.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

if pltfld == 'mean':
  img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
else:
  img = ax1.pcolormesh(xh, yh, AAsum, cmap=clrmp, vmin=rmin, vmax=rmax)

ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
if rmin < 0:
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
else:
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')

ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.03, 0.8, 0.06])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')

#btx = 'plot_snowfall_ant.py'
btx = 'plot_snowfall_ERA5.py'
bottom_text(btx, pos=[0.2, 0.01])


