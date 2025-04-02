"""
  Plot difference maps of monthly T/S bottom from different experiments
  Specify months (calendar numbering!) to average statistics by seasons
  Currently only 2 experiments are compared: expt-02 and expt-03 for daily forecasts

  If more than 1 year is given in keyargs than statistics are averaged over these years
  
   monthly mean and StDev  bottom T derived in calc_mnthlyTSbtm.py
  Save by years

  Usage: plot_diffTSbtm_mnthly.py --MMI=4 --YAS=1993 --YAE=1993 --MAS=10 --MAE=10 --varnm=temp
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
importlib.reload(mutob)
importlib.reload(manseas)

parser = argparse.ArgumentParser()
parser.add_argument("--MMI", help="Forecast initialization month ", type=int)
parser.add_argument("--YAS", help="Calendar (! not init.) year to start stat. averaging", type=int)
parser.add_argument("--YAE", help="Calendar year to end averaging of statistics", type=int)
parser.add_argument("--MAS", help="Calendar month to start averaging of statistics", type=int)
parser.add_argument("--MAE", help="Calendar month to end averaging of statistics", type=int)
parser.add_argument("--varnm", help="Variable to plot: temp / salin", type=str)
args = parser.parse_args()


# experiments: 2 - daily OB seasonal forecasts, 3 - same as 2 but with sea ice relaxation
# Default values: 
# Look at ens run #1 - the only ens. that has 5-day av. output fields
expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
nensR    = 1
# Averaging time period:
MMI   = 4    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YAS = 1993   # Start: f/cast init. year to use for monthly averaging
YAE = YAS   # End
MAS = 10
MAE = MAS

plot_fld = True # True - plot original fields from the expts, False - show only difference fld

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
if args.varnm:
  varnm = args.varnm

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

#DX, DY = mmom6.dx_dy(hlon, hlat)
#Acell  = DX*DY

TM = []
icc = 0
Nrec = 0
for YRA in (YAVRG):
  for MMA in range(MAS,MAE+1):
    print(f"Processing {YRA}/{MMA} {varnm} ...")

    dnmb0 = mtime.datenum([YRA, MMA, 15])
    YRI = manseas.yr_init_fcst_from_datenum(dnmb0, MMI)  # init year
    if YRI < 1993:
      continue    # cycle, outside the f/cast time period
    elif YRI > 2020:
      continue

    pthanls1 = pthseas['MOM6_NEP'][expt]['pthanls'].format(expt_nmb=enmb1)
    pthanls2 = pthseas['MOM6_NEP'][expt]['pthanls'].format(expt_nmb=enmb2)
    dflt1 = os.path.join(pthanls1,f'mnthly_Tbtm_expt{enmb1:02d}_{YRI}{MMI:02d}.pkl')
    dflt2 = os.path.join(pthanls2,f'mnthly_Tbtm_expt{enmb2:02d}_{YRI}{MMI:02d}.pkl')
    dfls1 = os.path.join(pthanls1,f'mnthly_Sbtm_expt{enmb1:02d}_{YRI}{MMI:02d}.pkl')
    dfls2 = os.path.join(pthanls2,f'mnthly_Sbtm_expt{enmb2:02d}_{YRI}{MMI:02d}.pkl')

    if varnm == 'temp':
      print(f'Loading T bottom stat  <-- {dflt1}')
      with open(dflt1, 'rb') as fid:
        Amn1, Astd1, TM = pickle.load(fid)
      print(f'Loading T bottom stat  <-- {dflt2}')
      with open(dflt2, 'rb') as fid:
        Amn2, Astd2, TM = pickle.load(fid)
    elif varnm == 'salin':
      print(f'Loading S bottom stat  <-- {dfls1}')
      with open(dfls1, 'rb') as fid:
        Amn1, Astd1, TM = pickle.load(fid)
      print(f'Loading S bottom stat  <-- {dfls2}')
      with open(dfls2, 'rb') as fid:
        Amn2, Astd2, TM = pickle.load(fid)
         
    MMF = manseas.mofcst_from_mocalend(YRI,MMI,MMA) # forecast month #
    imo = MMF-1

    AM1 = Amn1[imo,:,:]
    AM1 = np.expand_dims(AM1, axis=0)
    AM2 = Amn2[imo,:,:]
    AM2 = np.expand_dims(AM2, axis=0)

    AV1 = Astd1[imo,:,:]**2
    AV1 = np.expand_dims(AV1, axis=0)
    AV2 = Astd2[imo,:,:]**2
    AV2 = np.expand_dims(AV2, axis=0)

    # To get oevrall std from monthly STD: need to convert to var
    # If n(imo) of data in each months differs then:
    # var_tot = sum(n(imo)/n*(var(imo)) ), where sum(n(imo)) = n - total # of data in all months being processed
    # if nmb of data in processed months is the same, then var_tot = 1/nmonths*sum(var(imo))
    # For 5-day mean:
    ndays_mean = 5
    nmdays = mtime.month_days(MMA,YRA)
    nreci = np.floor((nmdays-1)/ndays_mean)  # nmb of records in ith month to derive std
    if icc == 0:
      AMN1 = AM1
      AMN2 = AM2
      AVAR1 = AV1*nreci 
      AVAR2 = AV2*nreci 
    else:
      AMN1 = AMN1 + AM1
      AMN2 = AMN2 + AM2
      AVAR1 = AVAR1 + AV1*nreci 
      AVAR2 = AVAR2 + AV2*nreci 

    icc += 1
    Nrec = Nrec + nreci                      # total nmb of recs


AMN1  = AMN1.squeeze()/icc
AMN2  = AMN2.squeeze()/icc
AVAR1 = AVAR1.squeeze()/Nrec
AVAR2 = AVAR2.squeeze()/Nrec

STD1 = np.sqrt(AVAR1)
STD2 = np.sqrt(AVAR2)

dltA = AMN1-AMN2  # new expt - old expt

# Mask out deep region:
AMN1 = np.where(HH<-500, 1.e3, AMN1)
AMN2 = np.where(HH<-500, 1.e3, AMN2)
STD1 = np.where(HH<-500, np.nan, STD1)
STD2 = np.where(HH<-500, np.nan, STD2)

#if varnm == 'iconc':
#  clrmp = mclrmps.colormap_conc()
#  clrmp.set_bad(color=[0.2, 0.2, 0.2])
#  rmin = 0.
#  rmax = 1.

#clrmp = mclrmps.colormap_ice_thkn()
#clrmp.set_bad(color=[0.2, 0.2, 0.2])
#rmin = -2.
#rmax = 4.

def plot_field(fgnmb, m, xR, yR, AMN1, STD1, clrmp, rmin, rmax, sttl=[], tscntrs=[], tslabels=[]):
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
    CS = ax1.contour(xR, yR, STD1, tscntrs, linestyles='solid', colors=[cntr_clr], linewidths=1)
    ax1.clabel(CS, tslabels, inline=1, fontsize=10)

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

#  ax3 = fig1.add_axes([0.02, 0.025, 0.8, 0.05])
#  ax3.text(0, 0, sinfo, fontsize=8)
#  ax3.axis('off')


  btx = 'plot_diffTSbtm_mnthly.py'
  bottom_text(btx, pos=[0.2, 0.01])


CLRS = [[0.6, 0.02, 0.6],
        [0.2, 0.38, 1],
        [0., 0.8, 0.5],
        [0.9, 0.6, 0],
        [1, 1, 1]]

clrmp_dlt = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBuGn_r')
clrmp_dlt.set_bad(color=[0.6,0.6,0.6])

if varnm == 'temp':
  clrmp = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[0.6,0.6,0.6])
  rmin = -2.
  rmax = 8.

  dmin = -0.8
  dmax = 0.8
else:
  clrmp = mutil.colormap_salin(clr_ramp=[0,0.2,0.5])
  clrmp.set_bad(color=[0.6, 0.6, 0.6])
  rmin = 31.0
  rmax = 35.5

  dmin = -1.2
  dmax = 1.2


tscntrs = [x/100 for x in range(0,800,30)]
tslabels = [x/100 for x in range(0,800,30)]


# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

sttl1 = (f'{expt_name1} {varnm}_btm init M={MMI} \n' + \
         f'avrg {YAS}-{YAE} mo={MAS}-{MAE}, cntrs=StDev')

sttl2 = (f'{expt_name2} {varnm}_btm init M={MMI} \n' + \
         f'avrg {YAS}-{YAE} mo={MAS}-{MAE}, cntrs=StDev')

sttl3 = (f'diff {varnm}_btm: {expt_name1}-{expt_name2}  init M={MMI} \n' + \
         f'avrg {YAS}-{YAE} mo={MAS}-{MAE}, cntrs=StDev')

plt.ion()

fgnmb=1
plot_field(1, m, xR, yR, dltA, STD1,clrmp_dlt,dmin,dmax,sttl=sttl3,tscntrs=tscntrs, tslabels=tslabels)
if plot_fld:
  plot_field(2, m, xR, yR, AMN1, STD1, clrmp, rmin, rmax, sttl=sttl1, tscntrs=tscntrs, tslabels=tslabels)
  plot_field(3, m, xR, yR, AMN2, STD2, clrmp, rmin, rmax, sttl=sttl2, tscntrs=tscntrs, tslabels=tslabels)



