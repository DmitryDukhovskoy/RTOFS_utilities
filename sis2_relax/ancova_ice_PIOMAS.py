"""
  Predict seasonal ice thickness from ice concentration 
  and previous season

  Use PIOMAS monthly fields 

  Y(i) = mu + tau(i) + a1*(ciconc(i)) + a2*(ciconc(i-1))
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
YRS = 1993
YRE = 2023
MMS = 2

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help=f"Year to start, default={YRS}", type=int)
parser.add_argument("--YRE", help=f"Year to end, default={YRE}", type=int)
parser.add_argument("--MMS", help=f"Month start, default={MMS}", type=int)
parser.add_argument("--MME", help=f"Month end, default={MMS}", type=int)
args = parser.parse_args()

if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMS:
  MMS = args.MMS
  MME = MMS
if args.MME:
  MME = args.MME

pthdata = '/Users/ddmitry/DATA/PIOMASv21'
varH = 'heff'
clrmpH = mclrmps.colormap_ice_thkn()
hmin = 0.
hmax = 4.
varC = 'area'
clrmpC = mclrmps.colormap_conc()
cmin = 0.
cmax = 1.

if MMS <= MME:
  MNTHS = np.arange(MMS,MME+1)
else:
  MNTHS = np.append(np.arange(MMS,13),np.arange(1,MME+1))

LON = LAT = None
icc = 0
for YR in range(YRS,YRE+1):
  flH = f'piomas_heff{YR}_v21.nc'
  dflH = os.path.join(pthdata, flH)
  print(f'Loading {dflH}')
  dsetH = xarray.open_dataset(dflH)
  flC = f'piomas_area{YR}_v21.nc'
  dflC = os.path.join(pthdata, flC)
  print(f'Loading {dflC}')
  dsetC = xarray.open_dataset(dflC)

  if LAT is None or LON is None:
    LAT  = dsetH['lat_scaler'].data
    LON  = dsetH['lon_scaler'].data

  for MM in MNTHS:
    # Find record #:
    Month = dsetH['month'].data
    Year  = dsetH['year'].data
    D     = np.sqrt((Month-MM)**2 + (Year-YR)**2)
    rindx = np.argmin(D)
    H2d   = dsetH['heff'].data[rindx,:].squeeze()
    H2d[0,:]  = 0.
    H2d[-1,:] = 0.
    H2d[:,0]  = 0.
    H2d[:,-1] = 0
    H2d = np.where(H2d>9999., np.nan, H2d)

    C2d = dsetC['area'].data[rindx,:].squeeze()
    C2d = np.where(np.isnan(H2d), np.isnan, C2d)

    H2d = np.expand_dims(H2d, axis=0)
    C2d = np.expand_dims(C2d, axis=0)
    if icc == 0:
      HI = H2d.copy()
      CI = C2d.copy()
    else:
      HI = np.append(HI, H2d, axis=0)
      CI = np.append(CI, C2d, axis=0)

    # Keep previous year:
    H2d_prev = H2d.copy()

    icc += 1 

i0 = 188
j0 = 57

Y = HI[:,j0,i0]
X = CI[:,j0,i0]

# Scale:
Xs = X/np.std(X)
Xs = Xs - np.mean(Xs)

Ys = Y/np.std(Y)
Ys = Ys - np.mean(Ys)

A = STOP 

clrmpH.set_bad(color=[0.2, 0.2, 0.2])
clrmpC.set_bad(color=[0.2, 0.2, 0.2])

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

sttl = f"PIOMAS {varnm} {YR}/{MM}" 
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

