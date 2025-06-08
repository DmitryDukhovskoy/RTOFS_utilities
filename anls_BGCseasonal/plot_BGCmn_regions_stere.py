"""
  BGC runs
  default output files, monthly means 

  Plot 2D mean BGC fields to analyze seasonal water mass structure
  stereographic projection
  in different regions: # CalCur - Calif Current region, Alaska - Alaska region, BerSea - Bering


  Variables:
  no3 - Nitrate, mol/kg
  o2  - Oxygen, mol/kg
  nlg - Large phytoplankton nitrogen, mol/kg
  nmd - Medium -"-   -"-, mol/kg
  nsm - Small -"-    -"-, mol/kg
  pdi - Diazotroph Phosphorus, mol/kg
  plg - Large phytoplankton (diatoms) phosphorus, mol/kg
  pmd - Medium -"-    -"-, mol/kg
  psm - Small -"-  -"-, mol/kg
  silg - Large phytoplankton silicon, mol/kg
  simd - Medium  -"-   -"-, mol/kg
  chl  - Chlorophyll, ug/kg - micrograms (1e-6) of chlor. per 1 kg of sea water
  nlgz - Large zooplankton nitrogen, mol/kg 
  nmdz - Medium zooplankton nitrogen, mol/kg
  nsmz - Small zooplankton nitrogen, mol/kg
  po4  - phospate, mol/kg
  nh4  - ammonia
  fed  - dissolved iron, mol/kg
  fedet - Detrital iron, iron contained in non-living particulate organic matter, mol/kg
  dissic - Dissolved inorganic Carbon Concentration, mol/m3 DIC=[CO2]+[HCO3-]+[CO2^3-]
  dissoc - Dissolved Organic Carbon Concentration, mol/m3,  DOC - fuels microbial activity 
  


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
import mod_colormaps as mclrmps
importlib.reload(mutob)
importlib.reload(manseas)

# Initial date
# Look at ens run #1 - the only ens. that has 5-day av. output fields
expt      = 'spinup'  # seasonal forecasts with dailyOB from SPEAR
varnm     = 'salin'  # temp (potential) / salin
YR        = 1995
MMS       = 12
regn_name = 'CalCur' # CalCur - Calif Current region, Alaska - Alaska region, 
                     # BeringChuk - Bering Sea and Chukchi Shelf
                     # Following Stoke et al., 2015
lr0  = 10 # ocean layers from 1, ..., 52 - reduced vertical grid
          # lr 10 = -47.5 m, lr 16 =-105 m, lr 23 = -212 m

parser = argparse.ArgumentParser()
parser.add_argument("--YR", help=f"Year to plot, default={YR}", type=int)
parser.add_argument("--MMS", help=f"Start Month to average, default={MMS}", type=int)
parser.add_argument("--MME", help=f"End Month to average, default={MMS}", type=int)
parser.add_argument("--regn", help=f"Region: CalCur, BeringChuk, Alaska, default={regn_name}", type=str)
parser.add_argument("--lr", help=f"Model layer to plot: 1,...,52, default={lr0}", type=int)
parser.add_argument("--varnm", help=f"Variable name to plot: salin, temp, default={varnm}", type=str)
args = parser.parse_args()

if args.YR:
  YR = args.YR
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
  
MAV = np.arange(MMS,MME+1)

# Not needed for spinup/hindcast:
nensR    = 1
expt_nmb = 2   # 2 - seas f/casts with dailyOB

if varnm == 'salin':
  ncvar = 'so'
elif varnm == 'temp':
  ncvar = 'thetao'

expt_name = f'NEPBGC_spinup'
run_info = f'{expt_name}, {varnm} {YR} M={MMS}:{MME}'

print(f'Plotting {varnm} {expt_name} ')
print(f'{run_info}')

fyaml = 'bgc_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
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

#mo_fcsts = manseas.yrmo_seasonal_fcst(YRS, MMS)

YR=1993
pthfcst0 = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YR=YR)
floceanm = 'ocean_month_z.nc'
dflin = os.path.join(pthfcst0,floceanm)
ds = xarray.open_dataset(dflin)
ZM = ds['z_l'].data
ZM = -abs(ZM)
nlrs = len(ZM)
assert(lr0<=nlrs), f"lr0={lr0} > Max layers in {floceanm} {nlrs}"
zz0 = ZM[lr0-1]

Time = []
iyr  = 0
icc = 0
for MM in MAV:
  dnmbS    = mtime.datenum([YR,MM,15])
  pthfcst0 = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YR=YR)
  ds = xarray.open_dataset(dflin)
  iz = lr0-1
  itime = MM-1
  A2d = ds[ncvar].isel(z_l=iz,time=itime).data

  if icc == 0:
    Asum = A2d.copy()
  else:
    Asum = Asum + A2d

  icc += 1

if icc > 1:
  A2d = Asum/icc

# Mask ocean > zmin depth:
#A2d = np.where( (np.isnan(A2d)) & (HH<0), -1.e3, A2d)
if lr0 > 2:
  A2d = np.where(HH>=zz0, np.nan, A2d)

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pths = safe_load(ff)
II = pths['ANLS_NEP'][regn_name]['II']
JJ = pths['ANLS_NEP'][regn_name]['JJ']
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

btx = 'plot_TSmn_regions_stere.py'
sttl = f"{run_info} z={zz0:.1f} m"

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



