"""
  Plot monthly ice conce from seas f/cast experiments
  Specify months (calendar numbering!) to average statistics by seasons

  Usage: plot_fcst_iconc_mnthly.py --MMI=4 --YRS=1993 --YRE=1993 --MMS=10 --MME=10
  use --help for more information on keywargs

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
import pickle

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
import mod_sis2_relax as msisrlx
importlib.reload(mutob)
importlib.reload(manseas)

parser = argparse.ArgumentParser()
parser.add_argument("--MMI", help="Forecast initialization month ", type=int)
parser.add_argument("--expt", help="f/cast experiment number: 2, 3", type=int)
parser.add_argument("--ensnmb", help="ensemble run number, 1,...,10", type=int)
parser.add_argument("--YRS", help="Calendar (! not init.) year to start stat. averaging", type=int)
parser.add_argument("--YRE", help="Calendar year to end averaging of statistics", type=int)
parser.add_argument("--MMS", help="Calendar month to start averaging of statistics", type=int)
parser.add_argument("--MME", help="Calendar month to end averaging of statistics", type=int)
args = parser.parse_args()

f_cntrobs = True   # Plot observation-derived ice edge

# experiments: 2 - daily OB seasonal forecasts, 3 - same as 2 but with sea ice relaxation
# Default values: 
expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
ensnmb = 1
# Default Averaging time period:
MMI   = 4    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YRS = 1993   # Start: f/cast init. year to use for monthly averaging
YRE = YRS   # End
MMS = 10
MME = MMS
expt_nmb = 3

plot_fld = True # True - plot original fields from the expts, False - show only difference fld

if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
else:
  YRE=YRS
if args.ensnmb:
  ensnmb=args.ensnmb
if args.MMI:
  MMI = args.MMI
if args.MMS:
  MMS = args.MMS
if args.MME:
  MME = args.MME
else:
  MME = MMS
if args.expt:
  expt_nmb=args.expt

# Determine init year for request averaging time period:
#dnmbS = mtime.datenum([YRS,MMS,15])
#YRS = manseas.yr_init_fcst_from_datenum(dnmbS, MMI)
#dnmbE = mtime.datenum([YRE,MME,15])
#YRE = manseas.yr_init_fcst_from_datenum(dnmbE, MMI)


if YRS == 1993 and MMI == 1:
  print(f"Requested start date of averaging: {YRS}/{MMS}")
  raise Exception("First initial month should be 4 for 1993, given MMI={MMI}")

YAVRG = [x for x in range(YRS,YRE+1)]

expt_name = f'NEPphys_frcst_dailyOB-expt{expt_nmb:02d}'
run_info = f'{expt_name} init MM={MMI} e{ensnmb:02d}, ice conc: {min(YAVRG)}-{max(YAVRG)}'


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
ndav       = pthseas['MOM6_NEP'][expt]['ndav']  # # of days output averaged

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape

pthoutp = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=expt_nmb)


TM = []
icc = 0
Nrec = 0
for YRA in (YAVRG):
  for MMA in range(MMS,MME+1):
    print(f"Processing {YRA}/{MMA} IceConc ...")

    dnmb0 = mtime.datenum([YRA, MMA, 15])
    YRI = manseas.yr_init_fcst_from_datenum(dnmb0, MMI)  # init year
    if YRI < 1993:
      continue    # cycle, outside the f/cast time period
    elif YRI > 2020:
      continue

    pthfcst = os.path.join(pthoutp,f'{YRI}-{MMI:02d}-e01','history')
    dcice = os.path.join(pthfcst,f'ice_month.nc')
    print(f'Reading {dcice}')

    MMF = manseas.mofcst_from_mocalend(YRI,MMI,MMA) # forecast month #
    imo = MMF-1

    ds = xarray.open_dataset(dcice)
    A2d = ds['siconc'].isel(time=imo).data.squeeze()

    if icc == 0:
      AMN = A2d
    else:
      AMN = AMN + A2d

    icc += 1

if icc > 1:
  AMN  = AMN.squeeze()/icc

# Mask out deep region:
#AMN1 = np.where(HH<-500, 1.e3, AMN1)
#AMN2 = np.where(HH<-500, 1.e3, AMN2)

#if varnm == 'iconc':
clrmp = mclrmps.colormap_conc()
clrmp.set_bad(color=[0.2, 0.2, 0.2])
rmin = 0.
rmax = 1.


def plot_field(fgnmb, m, xR, yR, A2d, clrmp, rmin, rmax, sttl=[]):
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(sttl)

  cntr_clr = [0.2, 0.7, 1.0]
  if fgnmb==1 : 
    clevel=0.
  else:
    clevel=0.15

  #CS = ax1.contour(xR, yR, A2d, [clevel], linestyles='solid', colors=[cntr_clr], linewidths=1)
  #ax1.clabel(CS, inline=1, fontsize=10)
    
  # extend: min, max, both
  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  if fgnmb>1:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')
  else:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')

  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'plot_fcst_iconc_mnthly.py'
  bottom_text(btx, pos=[0.2, 0.01])

  return ax1

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

mo_fcast = msisrlx.cal_mo_to_fcast(MMI,MMS)
mE_fcast = msisrlx.cal_mo_to_fcast(MMI,MME)
if MMS==MME and YRS==YRE:
  sttl = (f'{expt_name} IceConc init M={MMI} e{ensnmb:02d}\n' + \
         f'{YRS} mo={MMS} Lead: {mo_fcast}')
else:
  sttl = (f'{expt_name} IceConc init M={MMI} \n' + \
         f'avrg {YRS}-{YRE} mo={MMS}-{MME} Lead: {mo_fcast}-{mE_fcast}')

plt.ion()
fgnmb=1
ax1 = plot_field(1, m, xR, yR, A2d, clrmp,rmin,rmax,sttl=sttl)

if f_cntrobs:
  # Not averaged, for 1 year/month
  # Change if averaged contour is needed
  ci0=0.15
  cntr_clr = [1, 0.4, 0]
  ax1.contour(xR,yR,A2d,[ci0],linestyles='solid', colors=[cntr_clr], linewidths=1.2)
  # Use interpolated fields
  YR0 = YRS
  MM0 = MMS
  pthnsidc = f'/work/Dmitry.Dukhovskoy/data/NRT_NOAA_NSIDC_seaconc/{YR0}_mnth'
  fliceout = f'NSIDC_iconc_mnth_interpNEP816x342_{YR0}.nc'
  dfliceout = os.path.join(pthnsidc,fliceout)
  dset_nsidc = xarray.open_dataset(dfliceout)
  imo = MM0-1
  ICnrt = dset_nsidc['ice_conc'].isel(time=imo).data

  # Hgrid lon. lat:
  fyaml = 'paths_seasfcst.yaml'
  with open(fyaml) as ff:
    pthseas = safe_load(ff)

  expt = 'seasonal_daily'
  pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
  fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
  ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
  hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
  hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
  dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
  dfgrid_mom = os.path.join(pthtopo, fgrid)
  hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

  xRm, yRm = m(hlon,hlat)

  cntr_nsidc = [0.,0.2,0.6]
  ax1.contour(xRm,yRm,ICnrt,[ci0],linestyles='solid', colors=[cntr_nsidc], linewidths=1.5)

  ax3 = plt.axes([0.02, 0.12, 0.1, 0.1])
  dx = 0.15
  ax3.plot([0,dx],[0.1,0.1],'-',color=cntr_clr)
  ax3.text(2*dx, 0.1, 'F/cast', va='center')
  ax3.plot([0,dx],[0.2,0.2],'-',color=cntr_nsidc)
  ax3.text(2*dx,0.2, 'NSIDC NRT', va='center')
  ax3.set_xlim([0,0.7])
  ax3.set_ylim([0,0.3])
  ax3.axis('off')


