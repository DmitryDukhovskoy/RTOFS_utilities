"""
  Plot ithkn instant daily fields from
  RTOFS CICE4 forecasts

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
import mod_read_hycom as mhycom

init_date = 20250704  #
init_hr = 0
regn = 'north'
fhr = 0

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
#parser.add_argument("--init", help=f"init date", choices=[20250103, 20250704], required=True, type=int)
parser.add_argument("--init", help=f"init date", choices=[20251231, 20250704], required=True, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fhr", help=f"f/cast hour to plot", 
                   choices=[0,24,48,72,96,120,144,168,192], required=True, type=int)
args = parser.parse_args()

regn      = args.regn if args.regn else regn
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fhr       = args.fhr if args.fhr is not None else fhr


# Init date:
dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
yrI, mmI, ddI, hrI = mtime.datevec(dnmbI)[:4]


pthice = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/cice4_fcast/{init_date}"
if fhr == 0:
  flice = f"rtofs_glo.t{init_hr:02d}z.n00.cice_inst.nc"
else:
  flice = f"rtofs_glo.t{init_hr:02d}z.f{fhr:02d}.cice_inst.nc"
dflice = os.path.join(pthice, flice)

# Get grid
with xarray.open_dataset(dflice) as dcice:
  LON  = dcice["TLON"].data
  LAT  = dcice["TLAT"].data
  hice = dcice["hi"].data.squeeze()

JDIM, IDIM = LON.shape
JDIM = JDIM + 1   # ocean grid has + 1 row

# Read RTOFS topo:
# Note that RTOFS grid has +1 row at the top compared to CICE6
pthtopo = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/topo_grid/'
ftopo  = 'depth_GLBb0.08_09m11'
HH0 = mhycom.read_topo(pthtopo, ftopo, IDIM, JDIM)
HH0 = HH0[:-1,:]     # discard the extra row
#LMsk = np.where(HH0>=0, 0, 1)

# Date to plot:
dnmbP = dnmbI + fhr // 24
YR, MM, DD = mtime.datevec(dnmbP)[:3]

# Subset hemispheres:
if regn == 'north':
  ilat = np.argmax(np.any(LAT >= 50, axis=1))
  hlat = LAT[ilat:,:]
  hlon = LON[ilat:,:]
  A2d  = hice[ilat:,:]
  HH   = HH0[ilat:,:]  

elif regn == 'south':
  ilat = np.argmax(np.any(LAT >= -50, axis=1))
  hlat = LAT[:ilat,:]
  hlon = LON[:ilat,:]
  A2d  = hice[:ilat,:]
  HH   = HH0[ilat:,:]  

jdm, idm = hlon.shape

#A2d[HH >= 0] = np.nan

clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]


sttl = f"ithkn RTOFS-CICE4 init:{init_date}, lead time={fhr:03d}hrs\n {YR}/{MM:02d}/{DD:02d}"

sinfo = 'instant. ice concentration from operational RTOFS-CICE4\n'
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
#ax1.contour(xh, yh, HH, [0], linestyles='solid', colors=[cntr_clr], linewidths=1)

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

btx = 'plot_RTOFS_ithkn.py'
bottom_text(btx, pos=[0.2, 0.01])


