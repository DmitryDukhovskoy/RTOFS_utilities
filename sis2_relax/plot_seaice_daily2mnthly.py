"""
  Plot sea ice conc/thickness 
  from test simulations
  Monthly mean fields
  Wei's hindcast 

  Derive from daily ice_daily.nc 
 
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
import datetime as dt

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

derive_mean = True

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="year to plot: 1993, ..., 2020", type=int)
parser.add_argument("--mo", help="month to plot: 1,..., 12", type=int)
parser.add_argument("--varnm", help="field to plot: ithkn or iconc", type=str)
parser.add_argument("--test", help="test run name: irlx1, irlx2, ", type=str)
args = parser.parse_args()

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
varnm  = 'ithkn'  # iconc or ithkn
jday_plt = 0

varnm = args.varnm if args.varnm else None
yr_plt = args.yr if args.yr else 2012
test_name = args.test if args.test else 0
mo_plt = args.mo if args.mo else None

test_name = 'Wei_hindcast'

#pthtest = f'/work/Dmitry.Dukhovskoy/tmp/test_{test_name}'
#outfld  = 'icem'
pthtest = '/work/Dmitry.Dukhovskoy/tmp'
prfx = '20120101'  # 19930401 - time stamp used in SIS2 output in file names, note that find
           # closest archive file does not work for 19930401.icem*.nc file names
           # rename files using ./rename_archive_v0.sh 0 in the output dir

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

expt       = 'seasonal_daily'
pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
#ndav       = pthseas['MOM6_NEP'][expt]['ndav']  # # of days output averaged

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')


print(f'Test  Run: {test_name} Monthly mean Plot date: {yr_plt}/{mo_plt}')

YYR = yr_plt
MMR = mo_plt
mday1 = int(mtime.datenum([YYR,MMR,1]))
mday2 = mday1 + int(mtime.month_days(MMR,YYR))-1


# derive_mean:
Asum = None
icc  = 0
for dnmb0 in range(mday1,mday2+1):
  flice_name = f'{prfx}.ice_daily.nc'
  dfsis2 = os.path.join(pthtest, flice_name)
  YR0,MM0,DD0 = mtime.datevec(dnmb0)[:3]

  print(f'Reading {YR0}/{MM0:02d}/{DD0:02d}: {dfsis2}')

  with xarray.open_dataset(dfsis2, decode_times=False) as dset:
    dnmb_ref = mtime.datenum([1900,1,1])   # ref date in netcdf file
    TM = np.floor(dset['time'].data + dnmb_ref)
    dday = 0.5
    kday = np.where(np.abs(TM - dnmb0) < dday)[0]
    if kday.size == 0:
      raise ValueError(f"No time match of {dnmb0}")
    itime = kday[0] 
    HIce = dset['sithick'].isel(time=itime).data
    CIce = dset['siconc'].isel(time=itime).data

  if varnm == 'iconc':
    A2d = CIce
  elif varnm == 'ithkn':
    A2d = CIce*HIce

  if Asum is None:
    Asum = A2d.copy()
  else:
    Asum = Asum + A2d

  icc += 1

A2d = None
if icc>1:
  A2d = Asum / icc
else:
  A2d = Asum.copy()


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
  rmax = 4.


sttl = f"Test run {test_name} {varnm} monthly mean {YYR}/{MMR}"

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

sinfo = f'SIS2: {dfsis2}\n'

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)

ax1.set_title(sttl)

if varnm == 'ithkn':
  ax1.contour(xR, yR, A2d,[6, 10, 14], linestyles='solid', colors=[(0.95, 0.95, 0.95)])


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

try:
  btx = os.path.basename(__file__)
except:
  btx = 'plot_seaice_daily2mnthly.py'

bottom_text(btx, pos=[0.2, 0.01])


