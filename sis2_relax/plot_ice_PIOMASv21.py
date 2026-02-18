"""
  Plot ice thickness from PIOMAS ice reanalysis v2.1
  monthly fields
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
import argparse

PPTHN = '/Users/ddmitry/python'
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

varnm = 'ithck' # ithck or iconc
YR0 = 1994
MM0 = 2

parser = argparse.ArgumentParser()
parser.add_argument("--YR", help=f"Year to plot, default={YR0}", type=int)
parser.add_argument("--MM", help=f"Month to plot, default={MM0}", type=int)
parser.add_argument("--varnm", help=f"Variable name, ithck or iconc, default={varnm}", type=str)
args = parser.parse_args()


if args.YR:
  YR0 = args.YR
if args.MM:
  MM0 = args.MM
if args.varnm:
  varnm=args.varnm

pthdata = '/Users/ddmitry/DATA/PIOMASv21'
if varnm == 'ithck':
  flinp = f'piomas_heff{YR0}_v21.nc'
  varnc = 'heff'
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 4.
else:
  flinp = f'piomas_area{YR0}_v21.nc'
  varnc = 'area'
  clrmp = mclrmps.colormap_conc()
  rmin = 0.
  rmax = 1.

dflithkn = os.path.join(pthdata, flinp)

dset = xarray.open_dataset(dflithkn)
LAT  = dset['lat_scaler'].data
LON  = dset['lon_scaler'].data

# Find record #: days since 1901-01-01 = day=1, index=0
dnmb0 = mtime.datenum([YR0,MM0,1])
dnmbR = mtime.datenum([1901,1,1])
ndays = int(dnmb0-dnmbR) + 1
#Time  = dset['time'].data  # np datetime array
Month = dset['month'].data
Year  = dset['year'].data
D     = np.sqrt((Month-MM0)**2 + (Year-YR0)**2)
rindx = np.argmin(D)
A2d   = dset[varnc].data[rindx,:].squeeze()
A2d[0,:]  = 0.
A2d[-1,:] = 0.
A2d[:,0]  = 0.
A2d[:,-1] = 0

A2d = np.where(A2d>9999., np.nan, A2d)

clrmp.set_bad(color=[0.2, 0.2, 0.2])

Nclrs = clrmp.N


# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3300*1.e3, resolution='l',\
            projection='stere', lat_ts=55, lat_0=62, lon_0=-175)

# North Polar stereographic projection
#m = Basemap(projection='npstere',boundinglat=50,lon_0=0,resolution='l')

LON = np.where(LON<-900, np.nan, LON)
LAT = np.where(LAT<-900, np.nan, LAT)

xR, yR = m(LON, LAT)

ss1  = 'PIOMAS: ice edge 0.15 assimilated, no ice thickn assimilated\n'
ss2  = 'https://psc.apl.uw.edu/research/projects/piomas-20c/'
sinfo = ss1 + ss2



plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)

sttl = f"PIOMAS {varnm} {YR0}/{MM0}" 
ax1.set_title(sttl)


ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.025, 0.8, 0.05])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')


btx = 'plot_ice_piomas.py'
bottom_text(btx, pos=[0.2, 0.01])

