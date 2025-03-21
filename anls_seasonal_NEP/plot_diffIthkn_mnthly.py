"""
  Plot difference maps of monthly ice conce from different experiments
  Specify months (calendar numbering!) to average statistics by seasons
  Currently only 2 experiments are compared: expt-02 and expt-03 for daily forecasts

  If more than 1 year is given in keyargs than statistics are averaged over these years
  
   monthly mean and StDev  bottom T derived in calc_mnthlyTSbtm.py
  Save by years

  Usage: plot_diffIthkn_mnthly.py --MMI=4 --YAS=1993 --YAE=1993 --MAS=10 --MAE=10
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
import mod_rtofs as mrtofs
importlib.reload(mutob)
importlib.reload(manseas)

parser = argparse.ArgumentParser()
parser.add_argument("--MMI", help="Forecast initialization month ", type=int)
parser.add_argument("--YAS", help="Calendar (! not init.) year to start stat. averaging", type=int)
parser.add_argument("--YAE", help="Calendar year to end averaging of statistics", type=int)
parser.add_argument("--MAS", help="Calendar month to start averaging of statistics", type=int)
parser.add_argument("--MAE", help="Calendar month to end averaging of statistics", type=int)
args = parser.parse_args()


# experiments: 2 - daily OB seasonal forecasts, 3 - same as 2 but with sea ice relaxation
# Default values: 
# Look at ens run #1 - the only ens. that has 5-day av. output fields
expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
nensR    = 1
# Default Averaging time period:
MMI   = 4    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YAS = 1993   # Start: f/cast init. year to use for monthly averaging
YAE = YAS   # End
MAS = 10
MAE = MAS

plot_fld = True    # True - plot original fields from the expts, False - show only difference fld
plot_iconc = True  # True - show ice edge
plot_obsice = True # True - show NRT NSIDC ice edge

if args.YAS:
  YAS = args.YAS
if args.YAE:
  YAE = args.YAE
else:
  YAE=YAS
if args.MMI:
  MMI = args.MMI
if args.MAS:
  MAS = args.MAS
if args.MAE:
  MAE = args.MAE
else:
  MAE = MAS

# Determine init year for request averaging time period:
dnmbS = mtime.datenum([YAS,MAS,15])
YRS = manseas.yr_init_fcst_from_datenum(dnmbS, MMI)
dnmbE = mtime.datenum([YAE,MAE,15])
YRE = manseas.yr_init_fcst_from_datenum(dnmbE, MMI)


if YRS == 1993 and MMI == 1:
  print(f"Requested start date of averaging: {YAS}/{MAS}")
  raise Exception("First initial month should be 4 for 1993, given MMI={MMI}")

YAVRG = [x for x in range(YAS,YAE+1)]

enmb1=3
enmb2=2
expt_name1 = f'NEPphys_frcst_dailyOB-expt{enmb1:02d}'
expt_name2 = f'NEPphys_frcst_dailyOB-expt{enmb2:02d}'
run_info = f'{expt_name1} init MM={MMI} e{nensR:02d}, T/S bottom: {min(YAVRG)}-{max(YAVRG)}'


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

pthoutp1 = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=enmb1)
pthoutp2 = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=enmb2)


TM = []
icc = 0
Nrec = 0
for YRA in (YAVRG):
  for MMA in range(MAS,MAE+1):

    dnmb0 = mtime.datenum([YRA, MMA, 15])
    YRI = manseas.yr_init_fcst_from_datenum(dnmb0, MMI)  # init year
    if YRI < 1993:
      continue    # cycle, outside the f/cast time period
    elif YRI > 2020:
      continue

    pthfcst1 = os.path.join(pthoutp1,f'{YRI}-{MMI:02d}-e01','history')
    pthfcst2 = os.path.join(pthoutp2,f'{YRI}-{MMI:02d}-e01','history')
    dcice1 = os.path.join(pthfcst1,f'ice_month.nc')
    dcice2 = os.path.join(pthfcst2,f'ice_month.nc')

    MMF = manseas.mofcst_from_mocalend(YRI,MMI,MMA) # forecast month #
    imo = MMF-1

    print(f"Processing {YRA}/{MMA} IceThkn, init mnth={MMI} f/cast mnth={imo}")

    ds = xarray.open_dataset(dcice1)
    CM1 = ds['siconc'].isel(time=imo).data.squeeze()
    TM1 = ds['sithick'].isel(time=imo).data.squeeze()
    AM1 = CM1*TM1      # ice thickness average in a grid cell
    ds = xarray.open_dataset(dcice2)
    CM2 = ds['siconc'].isel(time=imo).data.squeeze()
    TM2 = ds['sithick'].isel(time=imo).data.squeeze()
    AM2 = CM2*TM2      # ice thickness average in a grid cell

    if icc == 0:
      AMN1 = AM1
      AMN2 = AM2
      CMN1 = CM1
      CMN2 = CM2
    else:
      AMN1 = AMN1 + AM1
      AMN2 = AMN2 + AM2
      CMN1 = CMN1 + CM1
      CMN2 = CMN2 + CM2

    icc += 1


AMN1  = AMN1.squeeze()/icc
AMN2  = AMN2.squeeze()/icc
CMN1  = CMN1.squeeze()/icc
CMN2  = CMN2.squeeze()/icc
dltA  = AMN1-AMN2  # new expt - old expt

if plot_obsice:
  CMobs,Xobs,Yobs = manseas.avrg_cice_NSIDC(YAS, YAE, MAS, MAE)

  # Get ice edge contour in the Bering Sea
  # get rid of ice in unneeded part of the domain
  CMobs[:,150:] = np.nan
  CMobs[:200,:] = np.nan

  CNTR = manseas.derive_ice_contour(CMobs, nmin=10)



# Mask out deep region:
#AMN1 = np.where(HH<-500, 1.e3, AMN1)
#AMN2 = np.where(HH<-500, 1.e3, AMN2)

#if varnm == 'iconc':
#clrmp = mclrmps.colormap_conc()
clrmp = mclrmps.colormap_ice_thkn()
clrmp.set_bad(color=[0.2, 0.2, 0.2])
rmin = 0.
rmax = 5.

clrmp_dlt = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='GnBu_r')
clrmp_dlt.set_bad(color=[0.2,0.2,0.2])
dmin = -2.
dmax = 2.

def plot_field(fgnmb, m, xR, yR, AMN1, CMN1, clrmp, rmin, rmax, \
               sttl=[], tscntrs=[], tslabels=[], tconc=[], CNTR=[], LON=[], LAT=[]):
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = m.pcolormesh(xR, yR, AMN1, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(sttl)

  cntr_clr = [0.2, 0.7, 1.0]
  if len(tscntrs) > 0:
    CS = ax1.contour(xR, yR, AMN1, tscntrs, linestyles='solid', colors=[cntr_clr], linewidths=1)
    ax1.clabel(CS, tslabels, inline=1, fontsize=10)

  iconc_clr = [1.,0.7,0.6]
  if tconc:
    ax1.contour(xR,yR,CMN1, [tconc], linestyles='solid', colors=[iconc_clr], linewidths=1)
 
  # Plot observed contour if provided:
  cobc_clr = [1.,0.,1]
  ncc = len(CNTR)
  for kcc in range(ncc):
    Ic = CNTR[kcc][:,0]
    Jc = CNTR[kcc][:,1]
# Get geodetic coordinates:
# use exact (float) indices to interpolate exact geodetic coordinate
    nic = len(Ic)
    Xc = np.zeros((nic))
    Yc = np.zeros((nic))
    for ipp in range(nic):
      ii0 = Ic[ipp]
      jj0 = Jc[ipp]
      xc0, yc0 = mrtofs.interp_indx2lonlat(ii0, jj0, LON, LAT)
      Xc[ipp] = xc0
      Yc[ipp] = yc0

    Xcm, Ycm = m(Xc,Yc)
    ax1.plot(Xcm, Ycm, linewidth=1, color=cobc_clr, linestyle='solid')

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

  if ncc>0:
    ax3 = plt.axes([0.02, 0.05, 0.1, 0.1])
    x1 = 0.1
    y1 = 0.2
    y2 = 0.1
    ax3.plot([0,x1],[y1,y1],'-',linewidth=1, color=cobc_clr)
    ax3.plot([0,x1],[y2,y2],'-',linewidth=1, color=iconc_clr)
    ax3.text(x1+0.04,y1-0.02,'iconc NSIDC')
    ax3.text(x1+0.04,y2-0.02,'iconc SIS2')
    ax3.set_xlim([0,0.4])
    ax3.set_ylim([0.0,0.3])
    ax3.axis('off')

  btx = 'plot_diffIthkn_mnthly.py'
  bottom_text(btx, pos=[0.2, 0.01])


# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

sttl1 = (f'{expt_name1} IceThkn init M={MMI} \n' + \
         f'avrg {YAS}-{YAE} mo={MAS}-{MAE}, cntrs=StDev')

sttl2 = (f'{expt_name2} IceThkn init M={MMI} \n' + \
         f'avrg {YAS}-{YAE} mo={MAS}-{MAE}, cntrs=StDev')

sttl3 = (f'diff IceThkn: {expt_name1}-{expt_name2}  init M={MMI} \n' + \
         f'avrg {YAS}-{YAE} mo={MAS}-{MAE}, cntrs=IceThkn')

tscntrs = [x for x in range(4,16,2)]
tslabels = [x for x in range(4,16,2)]

tconc = []
if plot_iconc:
  tconc = 0.15

plt.ion()

fgnmb=1
plot_field(1, m, xR, yR, dltA, CMN1, clrmp_dlt, dmin, dmax, sttl=sttl3, tconc=tconc)
if plot_fld:
  plot_field(2, m, xR, yR, AMN1, CMN1, clrmp, rmin, rmax, \
       sttl=sttl1, tscntrs=tscntrs, tslabels=tslabels, tconc=tconc, CNTR=CNTR, LON=Xobs, LAT=Yobs)
  plot_field(3, m, xR, yR, AMN2, CMN2, clrmp, rmin, rmax, \
       sttl=sttl2, tscntrs=tscntrs, tslabels=tslabels, tconc=tconc, CNTR=CNTR, LON=Xobs, LAT=Yobs)



