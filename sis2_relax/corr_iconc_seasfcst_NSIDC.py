"""
  Calc corr for NEP seasonal f/casts
 
  usage: corr_iconc_seasfcst_NSIDC.py --ensmb 1 --MMI 1 --YRS 1993 --YRE 2020 --MMS 1 --MME 3

"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
from yaml import safe_load
import argparse
import pandas as pd

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
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc
from mod_utils_fig import bottom_text
import mod_anls_seas as manseas
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="start year for corr. anls", type=int)
parser.add_argument("--YRE", help="end year for corr. anls", type=int)
parser.add_argument("--MMS", help="start month for corr. anls", type=int)
parser.add_argument("--MME", help="end month for corr. anls", type=int)
parser.add_argument("--MMI", help="init month of NEP f/cast: 1, ..., 12", type=int)
parser.add_argument("--ensmb", help="ensemble number, 1,..., 15", type=int)
parser.add_argument("--exptnmb", help="experiment number, 2,3", type=int)
args = parser.parse_args()

f_cntrobs = True   # Plot observation-derived ice edge
# Years in the relax file also used in the rlx file name:
YRS = 1993  # init yr
YRE = 2020
MMS = 1
MME = 12
MMI = 1     # init month
ifld = 'siconc' # partial area only
ens_nmb = 1  # NEP ensemble #
expt     = "seasonal_daily"
expt_nmb = 2


if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI
if args.ensmb:
  ens_nmb = args.ensmb
if args.MMS:
  MMS = args.MMS
  MME = MMS
if args.MME:
  MME = args.MME
if args.exptnmb:
 expt_nmb=args.exptnmb

varnm = ifld
runname  = f"NEPphys_frcst_dailyOB-expt{expt_nmb:02d}"

if MMI==1 and YRS==1993:
  YRS=1994

if MMI>1 and YRE==2020:
  YRE = 2019

nyrs_clim = 5
if nyrs_clim == 5:
  if MMI == 1:
    ICLIM=[[1993,1997],[1995,1999],[2000,2004],[2005,2009],[2010, 2014],[2015,2019],[2016,2020]]
  else:
    ICLIM=[[1993,1997],[1995,1999],[2000,2004],[2005,2009],[2010, 2014],[2015,2019]]

ICLIM=np.array(ICLIM)

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthoutp1 = pthseas['MOM6_NEP']['seasonal_daily']['pthoutp'].format(expt_nmb=expt_nmb)

expt     = "seasonal_daily"
pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape
LMsk = np.where(HH<0, 1, 0)

DX, DY = mmom6.dx_dy(hlon, hlat)
Acell  = DX*DY*1.e-6  # km2

# Mask out southern lats:
LMsk = np.where(hlat<53.,0,LMsk)
LMsk[:567,:] = 0
LMsk[:,:39] = 0
LMsk[:595,177:] = 0
# Mask for Bering Sea + Ber. Str. + S. Chukchi Shelf
BMsk = LMsk.copy()
BMsk[750:,189:] = 0
BMsk[:750,239:] = 0
JB,IB = np.where(BMsk==1)
# Mask for the Arctic Oc. part of the domain:
AMsk = LMsk.copy()
AMsk = np.where(BMsk==1, 0, AMsk)
JA,IA = np.where(AMsk==1)

#mcal = np.arange(mmi,mmi+12)
#mcal = np.where(mcal>12, mcal-12, mcal)
#d = abs(mcal-mm0)
#itime = np.argmin(d)

if MMS <= MME:
  cal_months = np.arange(MMS,MME+1)
else:
  cal_months = np.concatenate((np.arange(MMS,13), np.arange(1,MME + 1)))

# Select correct init year for given YR/MM
TM      = []
dInep   = []
dInsidc = []
for YR in range(YRS,YRE+1):
  iclm = np.where((ICLIM[:,0] <= YR) & (ICLIM[:,1] >= YR))[0][0]
  YRC1,YRC2 = ICLIM[iclm,:]

  if MMI==1 and YRC1==1993:
    YRC1=1994
    YRC2=1998

  # NEP clim:
  pthdump = pthseas['MOM6_NEP']['seasonal_daily']['pthsis2'].format(expt_nmb=expt_nmb)
  flclim = f'NEPseasfcast_siconc_clim_{YRC1}_{YRC2}_MI{MMI:02d}e{ens_nmb:02d}.nc'
  dfclim = os.path.join(pthdump,flclim)
  print(f'Opening {dfclim}')
  dclim_nep = xarray.open_dataset(dfclim)
  months_clim = dclim_nep['calend_months'].data
  # Climatology by cal. months 1, ..., 12
  # NSIDC clim
  dclim_nsidc = msisrlx.read_NSIDC_iconc_clim_interp(YR)
  for MM in cal_months:
    itime = MM-1
    # NEP seas f/casts
    # Find model init. year for given calendar month:
    dnmb     = mtime.datenum([YR,MM,15])
    YRI      = manseas.yr_init_fcst_from_datenum(dnmb, MMI)
    if YRI < 1993:
      continue
    pthfcst  = os.path.join(pthoutp1,f'{YRI}-{MMI:02d}-e01','history')
    dsis_ice = os.path.join(pthfcst,f'ice_month.nc')

    dfnep  = os.path.join(pthfcst,f'ice_month.nc')
    ds_nep = xarray.open_dataset(dfnep)
    tm_nep = ds_nep['time'].data
    tmP = pd.to_datetime(tm_nep)
    nrec = len(tmP)
    TNEP = np.zeros((nrec,3))
    for irec in range(nrec):
      yr0 = tmP.year[irec]
      mo0 = tmP.month[irec]
      dd0 = tmP.day[irec]
      TNEP[irec,:] = [yr0,mo0,dd0]

    mfcst = TNEP[:,1].astype(int)
    yrfcst = TNEP[:,0].astype(int)
    ifcst  = np.where((mfcst==MM) & (yrfcst==YR))[0][0]
    iclim  = np.where(months_clim==MM)[0][0]
    print(f'NEP: {dfnep}')
    print(f'Processing MMI={MMI} {YR}/{MM} f/cast mo={ifcst+1}, clim mo={iclim+1}')
    CInep = ds_nep['siconc'].isel(time=ifcst).data
    Clim_nep = dclim_nep['siconc'].isel(months=iclim).data

    # NSIDC interpolated
    pthnsidc = pthseas['ALL']['dirnsidc_intrp'].format(YR=YR)
    flnsidc = f'NSIDC_iconc_mnth_interpNEP{jdm}x{idm}_{YR}.nc'
    ds_nsidc = xarray.open_dataset(os.path.join(pthnsidc,flnsidc))
    CInsidc = ds_nsidc['ice_conc'].isel(time=itime).data
    Clim_nsidc = dclim_nsidc['ice_conc'].isel(months=itime).data

    # Anomalies:
    dInep.append(CInep-Clim_nep)
    dInsidc.append(CInsidc-Clim_nsidc)

    TM.append(dnmb)

dInep = np.array(dInep)
dInsidc = np.array(dInsidc)

TM = np.array(TM)
DV = mtime.datevec1D(TM, fHR=False)
DV = np.array(DV).transpose()

nrc = len(dInep)
#nyr = int(nrc/len(cal_months))
#assert(nyr*12==nrc), f'Check record length dInep={nrc}'

# Test pnt:
# Ber. Str
#i0=227
#j0=705
# W. Ber. Sea
#i0=165
#j0=691
i0=134
j0=708
a1 = dInep[:,j0,i0]
a2 = dInsidc[:,j0,i0]
# Plot ACF
#from statsmodels.graphics.tsaplots import plot_acf
#aplot_acf(ts, lags=40)  # You can adjust lags


def corr_stat(dInep, dInsidc):
  nrc = len(dInep)
  # Corr:
  sgm_sp = np.std(dInep, axis=0)
  sgm_ns = np.std(dInsidc, axis=0)
  mean_sp = np.mean(dInep, axis=0)
  mean_ns = np.mean(dInsidc, axis=0)

  sgm_sp = np.where(sgm_sp < 1e-10, 1.e-10, sgm_sp)
  sgm_ns = np.where(sgm_ns < 1e-10, 1.e-10, sgm_ns)

  A3d = dInep.copy()*0.
  for kk in range(nrc):
    Asp = dInep[kk,:,:]
    Ans = dInsidc[kk,:,:]
    A3d[kk,:,:] = (Asp-mean_sp)*(Ans-mean_ns)/(sgm_sp*sgm_ns)

  RR = 1./nrc*np.sum(A3d, axis=0)

  # Compute a p-value of t-statistics for correlation coeff. 
  # Compute t-statistic
  from scipy import stats
  t_stat = RR * np.sqrt(nrc - 2) / np.sqrt(1 - RR*RR)
  df = nrc - 2

  # Two-tailed p-value
  PVAL = 2 * (1 - stats.t.cdf(abs(t_stat), df))
  PVAL = np.where(RR==0, np.nan, PVAL)

  # Check for auto-correlation using 
  # Durbin-Watson test 
  # 0< DW <4
  # ~2	No autocorrelation
  # < 2	Positive autocorrelation
  # > 2	Negative autocorrelation
  from statsmodels.stats.stattools import durbin_watson
  DW_STAT = durbin_watson(dInep, axis=0)

  return RR, PVAL, DW_STAT

def plot_timeser(RR,dInep,dInsidc,j0,i0,sttl1):
  TX = np.arange(YRS,YRE+1)
  fig2 = plt.figure(2,figsize=(9,8))
  plt.clf()
  a1 = dInep[:,j0,i0]
  a2 = dInsidc[:,j0,i0]
  clr1 = [0.,0.5,0.8]
  clr2 = [0.9,0.4,0]
  ax21 = plt.axes([0.1, 0.5, 0.8, 0.4])
  ln1, = ax21.plot(TX, a1, '-', linewidth=2, color=clr1, label='NEP f/cast')
  ln2, = ax21.plot(TX, a2, '-', linewidth=2, color=clr2, label='NSIDC NRT')
  
  ax21.grid('on')
  
  ax21.set_title(sttl1)
  lgd = plt.legend(handles=[ln1,ln2], loc='upper left')

  btx = 'corr_iconc_seasfcst_NSIDC.py'
  bottom_text(btx)

  return ax21 


RR, PVAL, DW_STAT = corr_stat(dInep,dInsidc)

DW_STAT = np.where(np.isnan(RR), np.nan, DW_STAT)
DW_STAT = np.where(RR==0, 2., DW_STAT)


# Remove auto-correlation:
# t.s. are weekly correlated, this does not change much
#dInep_ind = np.diff(dInep, axis=0)
#dInsidc_ind = np.diff(dInsidc, axis=0)
#RR, PVAL, DW_STAT = corr_stat(dInep_ind,dInsidc_ind)

f_timeser = True
if f_timeser:
  # W Ber sea
#  i0=134
#  j0=708
  # E Ber Sea
  i0=157
  j0=650
  sttl1 = f'i={i0}, j={j0} Ice conc. anomalies, Months={MMS}-{MME}, R={RR[j0,i0]:.2f}\n'
  sttl1 = sttl1 + f'{runname} MI:{MMI:02d}e{ens_nmb:02d}'

  ax21 = plot_timeser(RR,dInep,dInsidc,j0,i0,sttl1)


from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

clrmp_dlt = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBu_r')
clrmp_dlt.set_bad(color=[0.2,0.2,0.2])
rmin = -1.
rmax = 1.

sinfo = 'Corr. from monthly anomlies of sea ice conc. wrt to 5-yr climatologies\n'
sinfo = sinfo + f'NEP MI{MMI:02d}e{ens_nmb:-2d} seas. f/casts vs NSIDC NRT v3 ice conc'

plt.ion()
fgnmb=1
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()

ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
#m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

cntr_clr = [0.2, 0.7, 1.0]
img = m.pcolormesh(xR, yR, RR, cmap=clrmp_dlt, vmin=rmin, vmax=rmax)
alf=0.05
ax1.contour(xR,yR,PVAL,[alf], linestyles='solid', colors=[cntr_clr], linewidths=1)

sttl = f'Xcorr & signif. {alf:.2f}  {runname} MI:{MMI:02d}e{ens_nmb:02d} vs NSIDC iconc anom\n'
sttl = sttl + f'{YRS}-{YRE} {MMS:02d}-{MME:02d}'
ax1.set_title(sttl)

# extend: min, max, both
ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(img, cax=ax2, orientation='vertical')

ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'corr_iconc_seasfcst_NSIDC.py'
bottom_text(btx)

plt_dw = False
if plt_dw:
  plt.ion()
  fig2 = plt.figure(2,figsize=(9,8))
  plt.clf()

  ax21 = plt.axes([0.1, 0.1, 0.8, 0.8])

  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  dw_clr = [0.2, 0.7, 1.0]
  dmin = 0
  dmax = 4
  img2 = m.pcolormesh(xR, yR, DW_STAT, cmap=clrmp_dlt, vmin=dmin, vmax=dmax)
  #dw_check = 2.  # no auto-correlation
  dw1 = 1.6
  dw2 = 2.4
  ax21.contour(xR,yR,DW_STAT,[dw1,dw2], linestyles='solid', colors=[dw_clr], linewidths=1.1)

  sttl2 = f'Durbin-Watson a/corr Test NEP iconc anom {YRS}-{YRE} {MMS}-{MME}\n'
  sttl2 = sttl2 + f'=2: no a/corr, <2 - pos. a/corr, >2 - neg. a/corr'
  ax21.set_title(sttl2)

  # extend: min, max, both
  ax22 = fig2.add_axes([ax21.get_position().x1+0.025, ax21.get_position().y0,
                     0.02, ax21.get_position().height])
  clb = plt.colorbar(img2, cax=ax22, orientation='vertical')

  ax22.yaxis.set_ticks(list(np.linspace(dmin,dmax,11)))
  ax22.set_yticklabels(ax22.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax2.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  bottom_text(btx)

  





