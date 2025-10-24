"""
  Plot ice thickn/conc  fields from sensitivity tests with
  atm.-forced UFS (datm UFS)

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

expt = 'ufs_datm_mx025_v02'
init_date = 20250103
init_hr = 0

# hs_h - grid cell mean (!) snow thickness, m
# snow_ai - snowfall rate cm/day (in liquid water equivalent !)
# dsnow_h - snow formation (cm/day) - can be > or < 0
# snoice_h - snow-ice formation (cm/day)
# melts_h  - top snow melt (cm/day)
parser = argparse.ArgumentParser()
parser.add_argument("--enmb", help="expt nunmber: 1, ...", type=int, required=True)
parser.add_argument("--init", help=f"init date, default={init_date}", type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--fday", help=f"forecast day to plot: 1,...,14, =0 - init. cond.", type=int, required=True)
parser.add_argument("--varnm", help="field to plot: iconc or ithkn", type=str)
args = parser.parse_args()

enmb      = args.enmb if args.enmb else None
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fday      = args.fday if args.fday is not None else None
varnm     = args.varnm if args.varnm else None
TLON = TLAT = LMSK = None

dlt_hr = 24.  # delta hours between saved/avrg output 

pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/ufs_datm_mx025/expt{enmb:02d}/cice6"

# Get date:
dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
dnmb0 = dnmbI + fday-1                              # day to plot
yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
nsec0 = hr0*3600

plot_init = dnmb0 == dnmbI


if varnm == 'iconc':
  if plot_init:
    varnc = 'aice'
  else:
    varnc = 'aice_d'
  units = 'fraction'
  clrmp = mclrmps.colormap_uv()
  clrmp = mclrmps.colormap_conc()
  rmin = 0.
  rmax = 1.
  sinfo = 'Ice area daily averaged aggregate\n'
elif varnm == 'ithkn':
  if plot_init:
    varnc = 'hi_d'
  else:
    varnc = 'hi'
  units = 'm'
  clrmp = mclrmps.colormap_temp()
  rmin = 0.
  rmax = 3.
  sinfo = 'grid cell mean ice thickness (ice vol per m2 of grid cell area)\n'
else:
  raise Exception('varnm {varnm} not recognized: iconc or ithkn')

clrmp.set_bad(color=[0.2, 0.2, 0.2])

sinfo = sinfo + pthoutp

if plot_init:
  flinp = f"iceh_ic.{yr0}-{mm0:02d}-{dd0:02d}-{nsec0:05d}.nc"
else:
  flinp = f"iceh.{yr0}-{mm0:02d}-{dd0:02d}.nc"
dflice = os.path.join(pthoutp,flinp)

print(f"Reading {varnm} from {dflice}")
with xarray.open_dataset(dflice) as ds:
  A2d = ds[varnc].isel(time=0).squeeze().data 
  if TLON is None:
    TLON = ds['TLON'].data
    TLAT = ds['TLAT'].data
    LMSK = ds['tmask'].data



clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.6,0.6,0.6]

if plot_init:
  sttl = f'{varnm} {units}, {expt}-expt{enmb:02d} init field {init_date}:{init_hr:02d}hr'
else:
  sttl = f'{varnm} {units}, {expt}-expt{enmb:02d} init:{init_date}:{init_hr:02d}hr daily av.:{yr0}/{mm0:02d}/{dd0:02d}'

print(f'Plotting {sttl} ...')
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

img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
if varnm == 'iconc':
  # Contour ice edge:
  CS = ax1.contour(xh, yh, A2d, [0.15], linestyles='solid', colors=[cntr_clr], linewidths=1)

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
btx = 'plot_datmUFS_ice_ant.py'
bottom_text(btx, pos=[0.2, 0.01])


