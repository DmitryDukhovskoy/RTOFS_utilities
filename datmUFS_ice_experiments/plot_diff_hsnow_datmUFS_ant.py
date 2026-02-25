"""
  Plot difference of snow depth fields from sensitivity tests with
  atm.-forced UFS (datm UFS) vs NASA SSM/I daily clim
  interpolated to mesh025 

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

expt = 'ufs_datm_mx025_v02'
init_date = 20250103
init_hr = 0
regn = 'south'
fld_avrg = "ice"  # SSM/I hsnow over ice area, i.e. m3 / m2_ice
fdaily = 1        # compare with daily (1) or monthly (0) SSM/I climatology

# hs_h - grid cell mean (!) snow thickness, m
# snow_ai - snowfall rate cm/day (in liquid water equivalent !)
# dsnow_h - snow formation (cm/day) - can be > or < 0
# snoice_h - snow-ice formation (cm/day)
# melts_h  - top snow melt (cm/day)
parser = argparse.ArgumentParser()
parser.add_argument("--enmb", help="expt nunmber: 1, ...", type=int, required=True)
parser.add_argument("--init", help=f"init date, default={init_date}", type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--fdaily", help=f"SSM/I climtology 1= daily, 0=mnth, default={fdaily}", type=int)
parser.add_argument("--fday", help=f"forecast day to plot: 1,...,14, =0 - init. cond.", 
                    type=int, required=True)
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--avrg", help=f"plot grid cell or ice area mean: cell or ice, default{fld_avrg}", type=str)
args = parser.parse_args()

enmb      = args.enmb if args.enmb else None
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fday      = args.fday if args.fday is not None else None
regn      = args.regn if args.regn else regn
fld_avrg  = args.avrg if args.avrg else fld_avrg
fdaily    = args.fdaily if args.fdaily is not None else fdaily

TLON = TLAT = LMSK = None
ssmi_daily = fdaily == 1   # compare with daily or monthly clim

pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/ufs_datm_mx025/expt{enmb:02d}/cice6"

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


# Get date:
plot_init = fday == 0  # initial conditions

dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
if plot_init:
  dnmb0 = dnmbI
else:
  dnmb0 = dnmbI + fday-1                              # day to plot

yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
YR,MM,DD = mtime.datevec(dnmb0)[:3]
nsec0 = hr0*3600

if plot_init:
  flinp = f"iceh_ic.{yr0}-{mm0:02d}-{dd0:02d}-{nsec0:05d}.nc"
else:
  flinp = f"iceh.{yr0}-{mm0:02d}-{dd0:02d}.nc"

pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/ufs_datm_mx025/expt{enmb:02d}/cice6"
dflice = os.path.join(pthoutp,flinp)

print(f"Processing {YR}/{MM}/{DD} expt{enmb:02d} {dflice}...")
 # grid cell mean snow thickness, m3/m2 of cell area
with xarray.open_dataset(dflice) as dcice:
  A2d  = dcice['hs_d'].data.squeeze()*100. # m --> cm 
  Aice = dcice['aice_d'].isel(time=0).squeeze().data 
  TLON = dcice['TLON'].data
  TLAT = dcice['TLAT'].data
  LMSK = dcice['tmask'].data

# CICE6 hsnow_cell = m3/m2_cell = sum(vsnon(n)*aicen(n))
# SSM/I hsnow_ssmi = m3/m2_ice, to compare:
# hsnow_ice = hsnow_cell / aice=sum(aicen(n)) m3/m2_cell * m2_cell/m2_ice = m3/m2_ice
if fld_avrg == 'ice':
  # Plot hsnow avrg over ice area:
  A2d = np.divide(A2d, Aice, out=np.zeros_like(A2d), where=Aice > 0)
A2d[LMSK==0] = np.nan

# Read interpolated snow depths:
# Snow depth climatology, Interpolated fields mesh025:
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]

if ssmi_daily:
  # Daily climatology:
  pthsnow = os.path.join(pthdata,'snow_nasa','daily_clim')
  flhsn = f"SSMI_hsnow_mesh025_1440x1080_dailyclim_{MM:02d}_south.nc" 
  fld_freq = 'daily'
else:
  # For CICE6 restart used monthly hsnow:
  pthsnow = os.path.join(pthdata,'snow_nasa','monthly_clim')
  #flhsn = 'SSMI_hsnow_mnthclim_1992_2007_mesh025_1440x1080_south.nc'
  flhsn = 'SSMI_hsnow_mnthclim_1998_2007_mesh025_1440x1080_south.nc'  # used this for creating restarts
  fld_freq = 'monthly'

dflhsn = os.path.join(pthsnow,flhsn)
print(f"Reading interpolated hsnow {dflhsn}")
with xarray.open_dataset(dflhsn) as ds_snow:
  units = ds_snow['snow_depth'].attrs.get('units')
  if units == 'm':
    cff = 100.     # m --> cm
  else:
    cff = 1.
  if ssmi_daily:
    HSi = cff * ds_snow['snow_depth'].isel(time=dd0-1).data.squeeze()
  else:
    HSi = cff * ds_snow['snow_depth'].isel(time=mm0-1).data.squeeze()

  LON = ds_snow['lon'].data
  LAT = ds_snow['lat'].data

# Note in CICE output snow depth is m3_snow/ m2_cell
# Data: m3_snow / m2_ice
#A2d_ice = np.divide(A2d, Aice, out=np.zeros_like(A2d), where=Aice != 0)
dHS = A2d - HSi

# Mask no ice:
dHS = np.where(Aice < 1e-11, 0., dHS)

units = 'cm'
clrmp = mclrmps.colormap_uv()
rmin = -50
rmax = 50
clrmp.set_bad(color=[0.2, 0.2, 0.2])


if fld_avrg == 'cell':
  sttl = f'diff hsnow m3/m2_cell, datmUFS expt{enmb:02d} vs SSM/I {fld_freq} clim\n init:{init_date}/{init_hr} fcast:{YR}/{MM:02d}/{DD:02d}'
  sinfo = 'difference datmUFS-climatology grid cell mean snow thickness, 100*(m3 per m2 od grid cell)\n'
else:
  sttl = f'diff hsnow m3/m2_ice, datmUFS expt{enmb} vs SSM/I {fld_freq} clim\n init:{init_date}/{init_hr} fcast:{YR}/{MM:02d}/{DD:02d}'
  sinfo = 'difference datmUFS-climatology ice area mean snow thickness, 100*(m3 per m2 of ice)\n'
if plot_init:
  sttl = sttl + ' INIT'
sinfo = sinfo + dflhsn + "\n" +  dflice 

plt.ion()

m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
#lons, lats = m.makegrid(idim, jdim) # get lat/lons of ny by nx evenly spaced grid.
#x, y = m(lons, lats) # compute map proj coordinates.
xh, yh = m(TLON,TLAT) # CICE6 coordinates

if regn == 'south':
  xl1 = -8.e6
  xl2 = -1.2e6
  yl1 = xl1
  yl2 = xl2

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

img = ax1.pcolormesh(xh, yh, dHS, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')
ax1.contour(xh, yh, Aice, [0.15], linestyles='solid', colors=[(0.2,0.9,0.2)], linewidths=1)

ax1.set_xlim([xl1, xl2])
ax1.set_ylim([yl1, yl2])
ax1.invert_yaxis()
ax1.invert_xaxis()

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

btx = 'plot_diff_hsnow_datmUFS_ant.py'
bottom_text(btx, pos=[0.2, 0.01])


