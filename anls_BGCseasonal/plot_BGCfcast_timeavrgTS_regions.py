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
parser.add_argument("--regnm", help="CalCur BeringChuk GulfAlaska", type=str)
parser.add_argument("--yrs", help="Year to start averaging", type=int, required=True)
parser.add_argument("--yre", help="Year to end averaging", type=int)
parser.add_argument("--mms", help="Calendar month, start of avrg", type=int, required=True)
parser.add_argument("--mme", help="Calendar month, end of avrg, default=yrs", type=int)
parser.add_argument("--minit", help="Initialization month 1,4,7,10", type=int, required=True)
parser.add_argument("--ensS", help="Start ens avergaing: Ensemb numb, 1,...,10", type=int, required=True)
parser.add_argument("--ensE", help="End ens avergaing, default=ensS no averaging", type=int)
parser.add_argument("--zz", help="approximate depth to plot: 0, ..., 5000 m", type=float, required=True)
args = parser.parse_args()

varnm = args.varnm if args.varnm else None
regn_name = args.regnm if args.regnm else 'CalCur'
YRS    = args.yrs if args.yrs else None
YRE    = args.yre if args.yre else YRS
MMS    = args.mms if args.mms else None
MME    = args.mme if args.mme else MMS
MMI    = args.minit if args.minit else None
ensS   = args.ensS if args.ensS else None
ensE   = args.ensE if args.ensE else ensS
zz_plt = args.zz if args.zz is not None else None # ocean depth to plot

zz_plt = -abs(zz_plt)

YAVRG = [x for x in range(YRS,YRE+1)]
MAVRG = [x for x in range(MMS,MME+1)]
ENSMB = [x for x in range(ensS,ensE+1)]

expt_name = 'NEPbgc_fcst_dailyOB01'
run_info = f'{expt_name} avrg {varnm} MI={MMI} {min(YAVRG)}-{max(YAVRG)} '+\
           f'Mo: {min(MAVRG)}-{max(MAVRG)} e{ensS:02d}-e{ensE:02d}'

print(f'Plotting {varnm} {expt_name} e-{ensS:02d}-e{ensE:02d} ')
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
ptharch = f'/archive/Dmitry.Dukhovskoy/fre/NEP/forecast_bgc/{expt_name}'

def read_mnth_mom6(YRI, MMI, yr0, mm0, ens_nmb, varnc, iselZ=None):
  """
    3 month runs are assumed for the hindcasts
  """
  pthfull = os.path.join(ptharch,f'{YRI}-{MMI:02d}-e{ens_nmb:02d}/history')
  flnm = f'oceanm_{yr0}_{mm0:02d}.nc'
  dflnm = os.path.join(pthfull,flnm)
  print(f'Reading {dflnm}')
  with xarray.open_dataset(dflnm) as ds_mom:
    if (iselZ is not None):
      AA = ds_mom[varnc].isel(zl=iselZ).data.squeeze()
    else:
      AA = ds_mom[varnc].data

  return AA

ZM = read_mnth_mom6(YRS, MMI, YRS, MMI, ensS, 'zl')
ZM = -abs(ZM)
dZ = np.abs(ZM-zz_plt)
ilr0 = np.argmin(dZ)
lr0  = ilr0+1
zz0 = ZM[ilr0]  # actual depth to be plotted

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

    Time.append(dnmb0)

    for ens_nmb in ENSMB:
      print(f'Reading ensmb e{ens_nmb:02d}')
      AA  = read_mnth_mom6(YRI, MMI, YR0, MM0, ens_nmb, varnc, iselZ=ilr0)

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

btx = 'plot_BGCfcast_timeavrgTS_regions.py'
sttl = f"{run_info} z={zz0:8.1f} m"

# Stereographic projection:
from mpl_toolkits.basemap import Basemap, cm

lon0, lat0, height, width = manseas.stereogr_params_regions(regn_name)

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



