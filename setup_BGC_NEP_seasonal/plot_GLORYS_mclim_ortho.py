"""
  Plot 2D fields from GLORYS daily climatologies
  for visual comparison with WOA23
  climatologies prepared in derive_GLORYS_clim.py

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import matplotlib.colors as colors
import pickle
import argparse

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
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')
import mod_swstate as msw
import conversions as gsw

from mod_utils_fig import bottom_text
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
importlib.reload(mutob)


parser = argparse.ArgumentParser()
parser.add_argument("--varnm", help="Variable to plot: temp/thetao, salin/so", type=str)
parser.add_argument("--month", help="month to plot: 1,...,12", type=int)
parser.add_argument("--day", help="day of the month: 1,...,31", type=int)
parser.add_argument("--depth", help="depth to plot, abs. values: 0, ..., 6000, default=0", type=float)
parser.add_argument("--interp", help="=1: interpolate to depth, default=0: closest layer", type=int)
args = parser.parse_args()

if args.varnm:
  varnm_in = args.varnm
if args.month:
  MM = args.month
if args.day:
  DD = args.day
finterp = args.interp or 0
zz0 = -abs(args.depth or 0)

YRS = 1993
YRE = 2024

if varnm_in == 'temp' or varnm_in == 'thetao':
  varnm = 'thetao'
elif varnm_in == 'salin' or varnm_in == 'so':
  varnm = 'so'

pthclim = '/archive/Dmitry.Dukhovskoy/datasets_NEP/GLORYS_clim'
flnm = f"GLORYS_CLIM_1993-2024_{MM:02d}{DD:02d}.nc"
dfclm = os.path.join(pthclim,flnm)

dset = xarray.open_dataset(dfclm)
ZM   = -abs(dset['depth'].values)
lon = dset['longitude'].values
lat = dset['latitude'].values
hlon, hlat = np.meshgrid(lon,lat)

assert(zz0 > ZM[-1]), f"Deepest level to plot should be shallower than {ZM[-1]:.2f}"

DLT = abs(ZM-zz0)
iz0 = np.argmin(DLT)
if finterp == 0:
  A2d = dset[varnm].isel(depth=iz0,  time=0).data
  zz0 = ZM[iz0]
else:
  if zz0 > ZM[0]:
    iz0 = 0
    A2d = dset[varnm].isel(depth=iz0, time=0).data
  else:
    if zz0>ZM[iz0]:
      iz1 = iz0-1
      iz2 = iz0
    else:
      iz1 = iz0
      iz2 = iz0+1

    zz1 = ZM[iz1]
    zz2 = ZM[iz2]
    A1 = dset[varnm].isel(depth=iz1, time=0).data
    A2 = dset[varnm].isel(depth=iz2, time=0).data
    assert(zz2 <= zz0 <= zz1), f"{zz0:.2f} must be between zz1={zz1} and zz2={zz2}" 
    A2d = A1*(zz2-zz0)/(zz2-zz1) + A2*(zz0-zz1)/(zz2-zz1)


if varnm == 'so':
  clrmp = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 28.0
  rmax = 38.0
elif varnm == 'thetao':
  clrmp = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[1,1,1])
  rmin = -2.
  rmax = 23.
elif varnm == 'zos':
  clrmp = mutil.colormap_ssh(nclrs=200)
  rmin = -0.5
  rmax = 0.5


# Get GLORYS topo proxi
# See: anls_seasonal_NEP/derive_GLORYS_topo_NEP.py
pthtopo   = '/work/Dmitry.Dukhovskoy/data/glorys_topo_NEP/'
dftopo = os.path.join(pthtopo,'GLORYS12_topoNEP_865x1321.pkl')
with open(dftopo,'rb') as fid:
  HH = pickle.load(fid)


# Set up orthographic projection
from mpl_toolkits.basemap import Basemap, cm


lon0 = 220.
lat0 = 50.
res  = 'l'
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
xR, yR = m(hlon, hlat)
PMsk = ( (xR > 1e20) | (yR > 1e20) )
AA = A2d.copy()
AA = np.where(HH >= 0, np.nan, AA)
#AA = np.insert(AA, 0, AA[:,-1], axis=1)
#AA = np.insert(AA, -1, AA[-1,:], axis=0)
AA[PMsk] = np.nan
xR[PMsk]   = 1.e30
yR[PMsk]   = 1.e30

ny, nx = A2d.shape
#AA = AA[0:ny, 0:nx]
#AA = np.where(HHG >= 0., np.nan, AA)
#xBND, yBND = m(IBND, JBND)

#import mod_colormaps as mclrmp

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
im1 = m.pcolormesh(xR, yR, AA, cmap=clrmp, vmin=rmin, vmax=rmax)

m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

#m.plot(xBND, yBND, 'r.')

sttl = f"GLORYS12v1, {varnm} climatology {YRS}-{YRE} {MM}/{DD}, z={zz0:.1f}m"
ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

sinfo=dfclm
ax3 = fig1.add_axes([0.1,0.05,0.8,0.02])
ax3.text(0,0, sinfo, fontsize=8)
ax3.axis('off')


btx = 'plot_GLORYS_mclim_ortho.py'
bottom_text(btx, fsz=6, pos=[0.05, 0.03])














