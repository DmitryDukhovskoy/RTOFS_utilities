"""
  Derive NSIDC sea ice concentration climatology
  interpolated to NEP grid

  calc_NSIDC_iconc_clim.py --YRS 2010 --YRE 2014
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
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="start year for deriving climat.: 1993, ..., 2022", type=int)
parser.add_argument("--YRE", help="end year for deriving climat.: 1993, ..., 2022", type=int)
args = parser.parse_args()

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
varnm  = 'iconc'

# Default values that can be modified by keywords
# Day to plot either in year days or actual date:
# Years in the relax file also used in the rlx file name:
YRS = 2011  # init yr
YRE = 2024
MMI = 1     # init month
f_save = True

if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE

# 5-yr climatologies:
#ICLIM=[[1993,1997],[1995,1999],[2000,2004],[2005,2009],[2010, 2014],[2015,2019],[2016,2020]]
#ICLIM=np.array(ICLIM)

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

expt     = "seasonal_daily"
pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

jdm, idm = hlon.shape
Asum = np.zeros((12,jdm,idm))
icc = 0
for YR in range(YRS,YRE+1):
  #pthnsidc = f'/work/Dmitry.Dukhovskoy/data/NRT_NOAA_NSIDC_seaconc/{YR}_mnth'
  pthnsidc = pthseas['ALL']['dirnsidc_intrp'].format(YR=YR)
  fliceout = f'NSIDC_iconc_mnth_interpNEP816x342_{YR}.nc'
  dfliceout = os.path.join(pthnsidc,fliceout)
  print(f'Loading {dfliceout}')
  dset = xarray.open_dataset(dfliceout)
  for MM in range(1,13):
    imo = MM-1
    CI = dset['ice_conc'].isel(time=imo).data
    Asum[imo,:,:] = Asum[imo,:,:] + CI

  icc += 1

Asum = Asum / icc

mcal = np.arange(1,13)
if f_save:
  darr_cice = xarray.DataArray(Asum, dims=("months","jdim","idim"),\
                     coords={"months": np.arange(12), \
                             "jdim": np.arange(jdm), \
                             "idim": np.arange(idm)})
  darr_months = xarray.DataArray(mcal, dims=("months"), \
                       coords={"months": np.arange(12)})
  dset = xarray.Dataset({"calend_months": darr_months, "ice_conc": darr_cice})

  # Add global attributes:
  dset.attrs.update({
    "info": "NSIDC NRT ice concentration climatology interpolated to NEP SIS2 grid",
    "code": "calc_NSIDC_iconc_clim.py"
  })

  pthclim = pthseas['ALL']['dirnsidc_clim']
  flclim  = f'NSIDC_NRT_interpNEP_iconc_clim_{YRS}_{YRE}.nc'
  dflclim = os.path.join(pthclim,flclim)
  print(f'Dumping climtology --> {dflclim}')
  dset.to_netcdf(dflclim, format='NETCDF3_64BIT', engine='netcdf4')


f_plt = False
if f_plt:
  clrmp = mclrmps.colormap_conc()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 1.
  if varnm == 'ithkn':
    clrmp = mclrmps.colormap_ice_thkn()
    clrmp.set_bad(color=[0.2, 0.2, 0.2])
    rmin = 0.
    rmax = 4.

  sttl = f"NSIDC NRT ice conc, {YR}/{MM:02d}"
  # Stereographic Map projection:
  from mpl_toolkits.basemap import Basemap, cm
  m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
              projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

  xR, yR = m(hlon, hlat)


  plt.ion()

  MM = 9
  imo = MM-1
  A2d = Asum[imo,:,:].squeeze()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  #m.contour(xR, yR, HH, [-1000], colors=[(0,0,0)], linestyles='solid')
  #ax1.axis('scaled')
  sttl(f'Ice Conc. Climatology M={MM}')
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

  sinfo='NSIDC NRT sea ice conc, NASA retrieval algorithm'
  ax3 = fig1.add_axes([0.02, 0.025, 0.8, 0.05])
  ax3.text(0, 0, sinfo, fontsize=8)
  ax3.axis('off')


  btx = 'plot_NSIDC_ice_stere.py'
  bottom_text(btx, pos=[0.2, 0.01])


