"""
  Plot snow depth daily fields from
  SFS runs
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

init_date = 20250701
init_hr = 0
regn = 'south'

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south", required=True, type=str)
parser.add_argument("--init", help=f"init date", choices=[20240701, 20250701], default=init_date, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fday", help=f"forecast day to plot: 0, 1, 2, ... , =0 - init. cond.", type=int, required=True)
parser.add_argument("--enmb", help="experiment number: 0, 1, 2, ...", required=True, type=int)
args = parser.parse_args()

regn      = args.regn if args.regn else regn
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fday      = args.fday if args.fday is not None else fday
enmb      = args.enmb

if fday < 0:
  raise RuntimeError(f"fday should be > 0")
  
# Get date:
plot_init = fday == 0  # initial conditions
  
dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
YRI,MMI,DDI,hrI = mtime.datevec(dnmbI, round_hrs=True)[:4]

if plot_init:
  dnmb0 = dnmbI
else:
  dnmb0 = dnmbI + fday-1                              # day to plot

yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
nsec0 = int(hr0 * 3600)
YR,MM,DD = mtime.datevec(dnmb0)[:3]


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

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/sfs_C192mx025_cice_test/expt{enmb:02d}/cice6"

if plot_init:
  flinp = f"iceh_ic.{yr0}-{mm0:02d}-{dd0:02d}-{nsec0:05d}.nc"
else:
  flinp = f"iceh.{yr0}-{mm0:02d}-{dd0:02d}.nc"
varnm = 'hs_d'

dflice = os.path.join(pthoutp,flinp)

print(f"Processing {YR}/{MM}/{DD} {hr0:02d}:00, SFS init {YRI}/{MMI:02d}/{DDI:02d} {hrI:02d}:00\n{dflice}")
with xarray.open_dataset(dflice) as dcice:
  A2d = dcice[varnm].data.squeeze()

A2d[HH >= 0] = np.nan

clrmp = mclrmps.colormap_temp()
rmin = 0.
rmax = 0.1

clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]

import mod_gfs_cice_anls as mgfscice
importlib.reload(mgfscice)
expt_name = mgfscice.sfs_tests_info(enmb)

sttl = f"hsnow SFS expt{enmb:02d} ({expt_name}) init {init_date}\n lead time={fday:02d}, {YR}/{MM:02d}/{DD:02d}"

sinfo = f'daily average snow depth from CICE6, {varnm}\n'
sinfo = sinfo + dflice


plt.ion()

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)

xh, yh = m(hlon,hlat) # GFS coords

# snow contours:
#cntrs = [x/100 for x in range(1,10,1)]
cntrs=[]

print("Plotting ...")

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])
#m.drawcoastlines()
m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)

img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

if len(cntrs) > 0:
  cs = ax1.contour(xh, yh, A2d, cntrs, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
  ax1.clabel(cs, inline=True, fontsize=10, fmt="%.2f")


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
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=14)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.03, 0.8, 0.06])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')

btx = 'plot_SFS_hsnow_daily.py'
bottom_text(btx, pos=[0.2, 0.01])


f_debug = False
if f_debug:
  plt.clf
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  ax1.pcolormesh(A2d, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

  ax1.set_aspect('equal', adjustable='box')
  ax1.set_ylim(800, 1080)
  ax1.set_xlim(100, 670)

  ii0 = 482
  jj0 = 1062



