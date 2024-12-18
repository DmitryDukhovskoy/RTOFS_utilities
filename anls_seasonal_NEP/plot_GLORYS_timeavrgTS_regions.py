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
varnm = 'so'  # thetao, so, zos 
# Averaging time period:
YAVRG = [x for x in range(2011,2021)]
MAVRG = [1,2,3]  # months to average:
regn_name = 'CalCur' # CalCur - Calif Current region, Alaska - Alaska region, BerSea - Bering
                     # Following Stoke et al., 2015
lr0  = 1  # ocean layers from 1, ..., 75
          # lr 31 =-102 m, lr 38 = -216 m


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
  print(f'Reading {varnm} <-- {dfl_glorys}')
  dset = xarray.open_dataset(dfl_glorys)

  for MM in (MAVRG):
    itime = MM-1
    idepth = lr0-1

    AA = dset[varnm].isel(time=itime, depth=idepth).data.squeeze()

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

rmin, rmax, tscntrs, tslabels = manseas.colormap_params(regn_name, varnm, zz0=zz0)

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
match regn_name:
  case 'CalCur':
    manseas.plot2D_CalCur(A2d, clrmp, rmin, rmax, xlim1, xlim2, ylim1, ylim2, \
                  fgnmb=1, btx=btx, tscntrs=tscntrs, tslabels=tslabels, sttl=sttl, \
                  hlon=hlon, hlat=hlat, HH=HH)

