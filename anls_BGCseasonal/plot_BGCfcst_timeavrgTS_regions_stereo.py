"""
  Plot 2D mean T/S fields to analyze seasonal water mass structure
  BGC seasona forecasts

  stereographic projection
  in different regions: # CalCur - Calif Current region, Alaska - Alaska region, BerSea - Bering
  Following Stoke et al., 2015

  in the Calif. Current region, see analysis:
  Auad et al., 2011
  The California Current System in relation to the Northeast Pacific Ocean circulation
  https://www.sciencedirect.com/science/article/pii/S0079661111001157

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
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_colormaps as mclrmps
importlib.reload(mutob)
importlib.reload(manseas)

parser = argparse.ArgumentParser()
parser.add_argument("--varnm", help="salin or temp", type=str, required=True)
parser.add_argument("--regnm", help="CalCur (default) or BeringChuk", type=str)
parser.add_argument("--mmi", help="f/cast init Month: 1,4,7,10", type=int, required=True)
parser.add_argument("--ens", help="ensemble run number: 1,...,10", type=int, required=True)
parser.add_argument("--yrs", help="Year to start averaging", type=int, required=True)
parser.add_argument("--yre", help="Year to end averaging", type=int, required=True)
parser.add_argument("--seas", help="WOA season: 13-JFM, 14-AMJ, 15-JJS, 16-OND, 0-annual", \
                    type=int, required=True)
parser.add_argument("--lr", help="MOM6 NEP vert layer 1,...,75", type=int, required=True)
args = parser.parse_args()

varnm = args.varnm if args.varnm else None
regn_name = args.regnm if args.regnm else 'CalCur'
MMI   = args.mmi if args.mmi else None
YRS   = args.yrs if args.yrs else None
YRE   = args.yre if args.yre else None
seas  = args.seas if args.seas else None
ens_nmb = args.ens if args.ens else None
lr0   = args.lr if args.lr else None # ocean layers from 1, ..., 75
                                     # lr 22 = -49.9 m, lr 31 =-102 m, lr 37 = -192 m

YAVRG = [x for x in range(YRS,YRE+1)]
match seas:
  case 13:
    MAVRG=[1,2,3]
  case 14:
    MAVRG=[4,5,6]
  case 15:
    MAVRG=[7,8,9]
  case 16:
    MAVRG=[10,11,12]
  case 0:
    MAVRG = [x for x in range(1,13)]
  case mm if 1 <= mm <= 12:
    MAVRG = [mm]
  case _:
    raise Exception(f"{seas} is not a valid value for season")

expt_nmb = 1
expt_name = f'NEPbgc_fcst_dailyOB{expt_nmb:02d}'
run_info = f'{expt_name} MMI={MMI} e{ens_nmb:02d}, avrg {varnm}: {min(YAVRG)}-{max(YAVRG)} Mo: {min(MAVRG)}-{max(MAVRG)}'

print(f'Plotting {varnm} {expt_name} ')
print(f'{run_info}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP']['seasonal_daily']['pthgrid']
fgrid      = pthseas['MOM6_NEP']['seasonal_daily']['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"]['seasonal_daily']["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

#mo_fcsts = manseas.yrmo_seasonal_fcst(YRS, MMI)

# Read vertical layers from oceanm archive file:
ptharch = f'/archive/Dmitry.Dukhovskoy/fre/NEP/forecast_bgc/NEPbgc_fcst_dailyOB01'

def read_mnth_mom6(YRI, MMI, yr0, mm0, varnc, isel1=None):
  pthfull = os.path.join(ptharch,f'{YRI}-{MMI:02d}-e{ens_nmb:02d}','history')
  flnm = f'oceanm_{yr0}_{mm0:02d}.nc'
  dflnm = os.path.join(pthfull,flnm)
  with xarray.open_dataset(dflnm) as ds_mom:
    if isel1 is not None:
      AA = ds_mom[varnc].isel(zl=isel1).data.squeeze()
    else:
      AA = ds_mom[varnc].data

  return AA

ZM = read_mnth_mom6(1994,1,1994,1,'zl')
ZM = -abs(ZM)
zz0 = ZM[lr0-1]

if varnm == 'salin':
  varnc = 'salt'
elif varnm == 'temp':
  varnc = 'potT'

Time = []
irec = 0
# 1995 - 2004
for YR0 in YAVRG:
  for MM0 in MAVRG:
    print(f'Reading {YR0}/{MM0}')
    dnmb0 = mtime.datenum([YR0,MM0,15])
    # Find init year for requested year and month given init month:
    YRI = manseas.yr_init_fcst_from_datenum(dnmb0, MMI,)
    AA  = read_mnth_mom6(YRI, MMI, YR0, MM0, varnc, isel1=lr0-1)
    Time.append(dnmb0)

    if irec == 0:
      A2d = AA.copy()
    else:
      A2d = A2d + AA

    irec += 1

A2d = A2d/irec
DV = mtime.datevec1D(Time)

# Mask ocean > zmin depth:
#A2d = np.where( (np.isnan(A2d)) & (HH<0), -1.e3, A2d)
if lr0 > 2:
  A2d = np.where(HH>=zz0, np.nan, A2d)

II = pthseas['ANLS_NEP'][regn_name]['II']
JJ = pthseas['ANLS_NEP'][regn_name]['JJ']
xlim1 = min(II)
xlim2 = max(II)
ylim1 = min(JJ)
ylim2 = max(JJ)

# Plot boundaries of the region:
lon_s = hlon[ylim1, xlim1:xlim2+1]
lat_s = hlat[ylim1, xlim1:xlim2+1]

lon_n = hlon[ylim2, xlim1:xlim2+1]
lat_n = hlat[ylim2, xlim1:xlim2+1]

lon_w = hlon[ylim1:ylim2+1, xlim1]
lat_w = hlat[ylim1:ylim2+1, xlim1]

lon_e = hlon[ylim1:ylim2+1, xlim2]
lat_e = hlat[ylim1:ylim2+1, xlim2]

Xreg, Yreg = mmisc.connect_segments([lon_w, lon_n, lon_e, lon_s], \
                                    [lat_w, lat_n, lat_e, lat_s])

rmin, rmax, tscntrs, tslabels = manseas.colormap_params(regn_name, varnm, zz0=zz0)

if varnm == 'salin' or varnm == 'salt': 
  clrmp = mclrmps.colormap_haline2()
  clrmp.set_bad(color=[0., 0., 0.])
#  clrmp.set_under(color=[0.6, 0.6, 0.6])
elif varnm == 'temp' or varnm == 'potT': 
  clrmp = mclrmps.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[0.,0.,0.])
elif varnm == 'ssh':
  clrmp = mclrmps.colormap_ssh(nclrs=200)
  rmin = -0.5
  rmax = 0.5

btx = 'plot_BGCtimeavrgTS_regions_stereo.py'
sttl = f"{run_info} z={zz0:8.1f} m"

# Stereographic projection:
from mpl_toolkits.basemap import Basemap, cm
match regn_name:
  case 'CalCur':
    width  = 4000*1.e3
    height = 4000*1.e3
    lat0   = 33.5
    lon0   = -128.
  case 'BeringChuk':
    width  = 3300*1.e3
    height = 3700*1.e3
    lat0   = 65.
    lon0   = -175.

m = Basemap(width=width, height=height, resolution='l',\
            projection='stere', lat_ts=55, lat_0=lat0, lon_0=lon0)

xR, yR = m(hlon, hlat)

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()

ax1 = manseas.plot_stereogr_axis(fig1, m, xR, yR, A2d, clrmp, rmin, rmax, \
                       btx=btx, tscntrs=tscntrs, tslabels=tslabels, sttl=sttl)

# Plot NEP domain:
plt.sca(ax1)
xdom, ydom = m(Xreg, Yreg)
m.plot(xdom, ydom, 'w-')



