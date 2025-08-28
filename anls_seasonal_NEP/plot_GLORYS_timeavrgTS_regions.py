"""
  Plot 2D mean T/S fields to analyze seasonal water mass structure
  in different regions: # CalCur - Calif Current region, Alaska - Alaska region, BerSea - Bering
  Following Stoke et al., 2015

  Use monthlg GLORYS interpolated onto NEP region

  in the Calif. Current region, see analysis:
  Auad et al., 2011
  The California Current System in relation to the Northeast Pacific Ocean circulation
  https://www.sciencedirect.com/science/article/pii/S0079661111001157

Data fields from Liz:
location for monthly GLORYS means for NEP region:
/archive/e1n/datasets/GLORYS/monthly_means/

location for padded monthly GLORYS means, concatenated by year used to generate NEP clim nudging files:
/archive/e1n/datasets/GLORYS/monthly_climatologies/

location for monthly GLORYS means, regridded to NEP for nudging as individual months:
/archive/e1n/mom6/NEP/sponge/monthly_sponge_files/

location for padded monthly GLORYS means, regridded to NEP for nudging and concatenated by year:
/archive/e1n/mom6/NEP/sponge/clims/

The last directory contains the files I used for nudging the solution to GLORYS. 

Daily GLORYS reanalysis for NEP domain prepared by Liz:
/archive/e1n/datasets/GLORYS/YYYY/nep_10

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import pickle
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
importlib.reload(mutob)

dnmb = mtime.datenum([1994,4,1])
expt = 'GLORYS_NEP'  # GLORYS extracted for NEP domain
varnm = 'salin'  # thetao, so, zos 
# Averaging time period:
YRS = 1995
YRE = 2004
YAVRG = [x for x in range(YRS,YRE+1)]
MMS = 1
MME = 3
#YAVRG = [1995]
MAVRG = [x for x in range(MMS,MME+1)]
regn_name = 'CalCur' # CalCur - Calif Current region, Alaska, BeringChuk
                     # Following Stoke et al., 2015
lr0  = 31  # ocean layers from 1, ..., 75
          # lr 22 = -47.5m, lr 31 =-102 m, lr 38 = -216 m


parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help=f"Start year to average, default={YRS}", type=int)
parser.add_argument("--YRE", help=f"End Year to average, default={YRS}", type=int)
parser.add_argument("--MMS", help=f"Start Month to average, default={MMS}", type=int)
parser.add_argument("--MME", help=f"End Month to average, default={MMS}", type=int)
parser.add_argument("--regn", help=f"Region: CalCur, BeringChuk, Alaska, default={regn_name}", type=str)
parser.add_argument("--lr", help=f"Model layer to plot: 1,...,52, default={lr0}", type=int)
parser.add_argument("--varnm", help=f"Variable name to plot: salin, temp, default={varnm}", type=str)
args = parser.parse_args()

if args.YRS:
  YRS = args.YRS
  YRE = YRS
if args.YRE:
  YRE = args.YRE
if args.regn:
  regn_name = args.regn
if args.lr:
  lr0 = args.lr
if args.MMS:
  MMS = args.MMS
  MME = MMS
if args.MME:
  MME = args.MME
if args.varnm:
  varnm = args.varnm

YAVRG = np.arange(YRS,YRE+1)
MAVRG = np.arange(MMS,MME+1)

if varnm == 'salin':
  ncvar = 'so'
elif varnm == 'temp':
  ncvar = 'thetao'


dv0 = mtime.datevec(dnmb)
YR0, MM0, DD0 = dv0[:3]

run_info = f'{expt} avrg {varnm}: {min(YAVRG)}-{max(YAVRG)} Mo: {min(MAVRG)}-{max(MAVRG)}'

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

# GLORYS monthly fields interpolated onto MOM6 NEP grid, use MOM topo
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


# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

icnt = 0
for YRS in (YAVRG):
  pthdata = gridfls[expt]["monthly"]['pthoutp']
  flnm = gridfls[expt]["monthly"]['fdata'].format(year=YRS)
  dfl_glorys = os.path.join(pthdata,flnm)
  print(f'Reading {ncvar} <-- {dfl_glorys}')
  dset = xarray.open_dataset(dfl_glorys)

  for MM in (MAVRG):
    itime = MM-1
    idepth = lr0-1

    AA = dset[ncvar].isel(time=itime, depth=idepth).data.squeeze()

    if icnt == 0:
      A2d = AA.copy()
      ZM = dset['depth'].data.squeeze()
      ZM = -abs(ZM)
    else:
      A2d = A2d + AA

    icnt += 1

A2d = A2d/icnt
zz0 = ZM[lr0-1]

A2d = np.where(HH > -.1, np.nan, A2d)
# Mask bottom:
if zz0 < -5.:
  A2d = np.where(HH>=zz0, np.nan, A2d)

II = pthseas['ANLS_NEP'][regn_name]['II']
JJ = pthseas['ANLS_NEP'][regn_name]['JJ']
xlim1 = min(II)
xlim2 = max(II)
ylim1 = min(JJ)
ylim2 = max(JJ)

rmin, rmax, tscntrs, tslabels = manseas.colormap_params(regn_name, ncvar, zz0=zz0)

if varnm == 'salin' or varnm == 'salt' or varnm == 'so':
  clrmp = mclrmps.colormap_haline2()
  clrmp.set_bad(color=[0., 0., 0.])
#  clrmp.set_under(color=[0.6, 0.6, 0.6])
elif varnm == 'temp' or varnm == 'potT' or varnm == 'thetao':
  clrmp = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[0.,0.,0.])
elif varnm == 'ssh':
  clrmp = mutil.colormap_ssh(nclrs=200)
  rmin = -0.5
  rmax = 0.5

btx = 'plot_GLORYS_timeavrgTS_regions.py'
sttl = f"{run_info} z={zz0:8.1f} m"
#match regn_name:
#  case 'CalCur':
#    manseas.plot2D_CalCur(A2d, clrmp, rmin, rmax, xlim1, xlim2, ylim1, ylim2, \
#                  fgnmb=1, btx=btx, tscntrs=tscntrs, tslabels=tslabels, sttl=sttl, \
#                  hlon=hlon, hlat=hlat, HH=HH)


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


