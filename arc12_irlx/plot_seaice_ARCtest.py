"""
  Plot sea ice conc/thickness 
  from test simulations

  Plot daily ice fields (icem*.nc) 

 
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
parser.add_argument("--yr", help="year to plot: 1993, ..., 2020", type=int)
parser.add_argument("--jday", help="year day to plot: 1, ..., 366", type=int)
parser.add_argument("--mo", help="month to plot: 1,..., 12", type=int)
parser.add_argument("--day", help="day to plot: 1, ..., 31", type=int)
parser.add_argument("--varnm", help="field to plot: ithkn or iconc", type=str)
parser.add_argument("--ntest", help="test run nmb: 1, ..., ", type=str)
args = parser.parse_args()

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
varnm  = 'ithkn'  # iconc or ithkn
test_name = 'irlx3'
jday_plt = 0

if args.varnm:
  varnm = args.varnm
if args.yr:
  yr_plt = args.yr
if args.jday:
  jday_plt = args.jday
if args.ntest:
  test_nmb = args.ntest
if args.mo:
  mo_plt = args.mo
if args.day:
  day_plt = args.day


test_name = f'ARC_test{test_nmb}'
pthtest = f'/work/Dmitry.Dukhovskoy/tmp/{test_name}'
outfld  = 'icem'
prfx = ''  # 19930401 - time stamp used in SIS2 output in file names, note that find
           # closest archive file does not work for 19930401.icem*.nc file names
           # rename files using ./rename_archive_v0.sh 0 in the output dir

if jday_plt >0 and jday_plt <=366:
  dnmbR  = mtime.jday2dnmb(yr_plt,jday_plt)
else:
  dnmbR  = mtime.datenum([yr_plt,mo_plt,day_plt])  # day to plot


dvR     = mtime.datevec(dnmbR)
dnmb0   = dnmbR
ndav = 1
jF0 = 377
iF0 = 24
i0  = iF0-1
j0  = jF0-1

# Find closest output:
YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthtest, dnmbR, fld=outfld)
dv0  = mtime.datevec(dnmb0)
YR0, MM0, DD0 = dv0[:3]
jday0   = int(mtime.date2jday([YR0,MM0,DD0]))


print(f'Test  Run: {test_name} Plot date: {dvR[0]}/{dvR[1]}/{dvR[2]}')

if len(prfx) > 0:
  flice_name = f'{prfx}.icem_{YR0}_{jday0:03d}.nc'
else:
  flice_name  = f'icem_{YR0}_{jday0:03d}.nc'
dfsis2 = os.path.join(pthtest, flice_name)

print(f'Reading {dfsis2}')

dset   = xarray.open_dataset(dfsis2)

HIce = dset['sithick'].isel(time=0).data
CIce = dset['siconc'].isel(time=0).data
if varnm == 'iconc':
  A2d = CIce
elif varnm == 'ithkn':
  A2d = CIce*HIce


# -------------------
#
# Plot ice fields
#
# -------------------
if varnm == 'iconc':
  clrmp = mclrmps.colormap_conc()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 1.
elif varnm == 'ithkn':
  clrmp = mclrmps.colormap_ice_thkn()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 5.

# Averaging period:
dnmb_av1 = dnmb0 - np.floor(ndav/2)
#if dnmb_av1 < dnmbI: dnmb_av1=dnmbI
dnmb_av2 = dnmb_av1 + ndav-1

dv_av1 = mtime.datevec(dnmb_av1)
yrs, mms, dds = dv_av1[:3]
dv_av2 = mtime.datevec(dnmb_av2)
yre, mme, dde = dv_av2[:3]

sttl = f"Test run {test_name} {varnm} avrg: {yrs}/{mms}/{dds}-{yre}/{mme}/{dde}"
if j0 >= 0 and i0 >= 0:
  sttl = sttl + f"\n Test pnt iF0/jF0 = {iF0}/{jF0} {varnm}={A2d[j0,i0]:.6f}"

# Stereographic Map projection:
#from mpl_toolkits.basemap import Basemap, cm
#m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
#            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

#xR, yR = m(hlon, hlat)

sinfo = f'SIS2: {dfsis2}\n'

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
#m.drawcoastlines()
#m.drawparallels(np.arange(-90.,120.,10.))
#m.drawmeridians(np.arange(-180.,180.,10.))

img = ax1.pcolormesh(A2d, cmap=clrmp, vmin=rmin, vmax=rmax)

ax1.set_title(sttl)

if varnm == 'ithkn':
  ax1.contour(A2d,[1,2,3,4,5], linestyles='solid', colors=[(0.95, 0.95, 0.95)])

# Show test pnt:
if j0 >= 0 and i0 >=0:
  #xTst, yTst = m(hlon[j0,i0],hlat[j0,i0])
  #ax1.plot(xTst,yTst,'o')
  ax1.plot(i0,j0,'o')

ax1.axis('scaled')


ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.025, 0.8, 0.05])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')


btx = 'plot_seaice_ARCtest.py'
bottom_text(btx, pos=[0.2, 0.01])


