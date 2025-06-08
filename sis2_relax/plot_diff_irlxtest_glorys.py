"""
  Plot difference in ocean temp fields from ice relax experiments
  and GLORYS
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
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')
import mod_swstate as msw
import conversions as gsw

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
importlib.reload(mutob)

YR0  = 1993
MM0  = 10
lr0  = 1 
expt = 'irlx3'

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help=f"year to plot: 1993, ..., 2020, default={YR0}", type=int)
parser.add_argument("--mo", help=f"month to plot: 1,...,12, default={MM0}", type=int)
parser.add_argument("--lr", help=f"layer to plot: 1,...,75, default={lr0}", type=int)
parser.add_argument("--irlx", help=f"irlx test expert nmb, 3,...,9 note not all have ocean outp", type=int)
args = parser.parse_args()

if args.yr:
  YR0= args.yr
if args.mo:
  MM0 = args.mo
if args.lr:
  lr0 = args.lr
if args.irlx:
  expt = f'irlx{args.irlx}'

dnmb0 = mtime.datenum([YR0,MM0,15])

# 
varnm = 'temp'  # temp or salin

dv0 = mtime.datevec(dnmb0)

print(f'{expt} {varnm} lr={lr0} {YR0}/{MM0}')

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

pthtopo    = gridfls['MOM6_NEP']['seasonal_fcst']['pthgrid']
fgrid      = gridfls['MOM6_NEP']['seasonal_fcst']['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"]["seasonal_fcst"]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

HH = dstopo_nep['depth'].data
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

match varnm:
  case "temp":
    varnm_mom = 'potT'
    varnm_glorys = 'thetao'
  case "salin":
    varnm_mom = 'salt'
    varnm_glorys = 'so'

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

# Not all test experiments have ocean fields saved:
pthtest=f'/work/Dmitry.Dukhovskoy/tmp/test_{expt}'
ndav = 5 # N-day average fields

YYR = dv0[0]
MMR = dv0[1]
mday1 = int(mtime.datenum([YYR,MMR,1]))
mday2 = mday1 + int(mtime.month_days(MMR,YYR))-1
Asum = None
icc = 0
for dnmbR in range(mday1+1,mday2+1,ndav):
  # Find closest output:
  YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthtest, dnmbR, fld='oceanm')
  dv0  = mtime.datevec(dnmb0)
  YR0, MM0, DD0 = dv0[:3]
  jday0   = int(mtime.date2jday([YR0,MM0,DD0]))

  if MM0 != MMR:
    print(f' found month {MM0} requested {MMR}, skipping ...')
    continue

  flname = f'oceanm_{YR0}_{jday0:03d}.nc'
  dfmom = os.path.join(pthtest, flname)

  print(f'Reading {YR0}/{MM0:02d}/{DD0:02d}: {dfmom}')

  dset = xarray.open_dataset(dfmom)
  A2d  = dset[varnm_mom].isel(time=0,zl=lr0-1).data

  if Asum is None:
    Asum = A2d.copy()
  else:
    Asum = Asum + A2d

  icc += 1

Tmom = None
if icc>1:
  Tmom = Asum / icc
else:
  Tmom = Asum.copy()

# 
# Get GLORYS:
pthglorys = '/archive/e1n/mom6/NEP/sponge/glorys/nep_10k'
flglorys  = f'glorys_monthly_NEP_sponge_{YR0}_clim.nc'
dflgl = os.path.join(pthglorys,flglorys)
print(f'Reading {dflgl}')
dset = xarray.open_dataset(dflgl)
Time = dset['time'].data
TM = mmisc.convert_nptime_to_datenum(Time)
D = abs(TM-dnmb0)
itime = np.argmin(D)
dv0 = mtime.datevec(TM[itime])
assert dv0[0]==YR0, f'Requested YR={YR0}, year in rlx file={dv0[0]}'
assert dv0[1]==MM0, f'Requested month={MM0}, month in rlx file={dv0[1]}'

Tglorys = dset[varnm_glorys].isel(time=itime,depth=lr0-1).data
Tglorys = np.where(HH>=0, np.nan, Tglorys)


dltT = Tmom - Tglorys

clrmp = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBu_r')
clrmp.set_bad(color=[0.2,0.2,0.2])
rmin = -0.5
rmax = 0.5


# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

sttl = f'(MOM6-SIS2_{expt}-GLORYS)  diff{varnm} lr={lr0} {YR0}/{MM0}'

tscntrs = [-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5]
tslabels = tscntrs


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))
img = m.pcolormesh(xR, yR, dltT, cmap=clrmp, vmin=rmin, vmax=rmax)
#CS = ax1.contour(xR,yR,dltT,tscntrs, linestyles='solid', linewidths=1, colors=[(0., 0., 0.)])
#ax1.clabel(CS, tslabels,inline=1, fontsize=10)
ax1.set_title(sttl)

# extend: min, max, both
ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')

ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


btx = 'plot_diff_irlxtest_glorys.py'
bottom_text(btx, fsz=6, pos=[0.08, 0.03])



