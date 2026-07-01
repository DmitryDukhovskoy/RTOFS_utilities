"""
  Plot iconc daily fields from
  GFSv17 or GFSv16 f/casts
  GFSv16 daily avrg and interpolated onto mesh025 grid

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

init_date = 20251231  # 20240715 or 20251231
init_hr = 0
regn = 'south'
fday = 16

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--init", help=f"init date", choices=[20240715,20251231], required=True, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fday", help=f"f/cast day to plot: 1-16, default={fday}", type=int)
parser.add_argument("--gfs", help=f"GFS version to plot", choices=[16,17], required=True, type=int)
args = parser.parse_args()

regn      = args.regn if args.regn else regn
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fday      = args.fday if args.fday is not None else fday
gfsv      = args.gfs if args.gfs else None  

if fday < 0:
  raise RuntimeError(f"fday should be > 0")

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

# Init date:
dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb

match gfsv:
  case 16:
    pthgfs = f'/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/GFSv16/{init_date}/intrp_mesh025'
    flgfs = f'gfsv16_iconc_init{init_date}_days000_016.nc'
    varnm = 'ICEC_surface'  # m
  case 17:
    pthgfs = f'/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/GFSv17/{init_date}/daily'
    flgfs = f'gfsv17_iconc_init{init_date}_days000_016.nc'
    varnm = 'aice_h'


dflice = os.path.join(pthgfs, flgfs)

print(f"Processing {dflice}")
with xarray.open_dataset(dflice) as dcice:
  A2d  = dcice[varnm].isel(time=fday-1).data.squeeze()
  TM  = dcice['time'].isel(time=fday-1).values

A2d[HH >= 0] = np.nan

# Date to plot:
dt = TM.astype('datetime64[ms]').astype(object)
YR = dt.year
MM = dt.month
DD = dt.day


clrmp = mclrmps.colormap_conc()
rmin = 0.
rmax = 1.

clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]

sttl = f"iconc GFSv{gfsv} init {init_date}, lead time={fday:02d} days\n {YR}/{MM:02d}/{DD:02d}"

sinfo = 'daily average ice concentration from CICE6\n'
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


print("Plotting ...")

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])
#m.drawcoastlines()
m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)

img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

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
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=16)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.03, 0.8, 0.06])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')

btx = 'plot_GFS_iconc_daily.py'
bottom_text(btx, pos=[0.2, 0.01])


