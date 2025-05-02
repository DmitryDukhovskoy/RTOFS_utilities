"""
  Calc RMSE, ice extent and ice area 
  Using mean ensmb iconc anomalies
  see: calc_SPEAR_iconc_anom_ensmean.py
 
  usage: corr_iconc_SPEARmnensmb_NSIDC.py --MMI 1 --YRS 1993 --YRE 2020 --MMS 1 --MME 3

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
parser.add_argument("--MMI", help="init month of SPEAR f/cast: 1, ..., 12", type=int)
args = parser.parse_args()

# Years in the relax file also used in the rlx file name:
YRS = 1993  # init yr
YRE = 2020
MMS = 1
MME = 12
MMI = 1     # init month
ifld = 'siconc' # partial area only

if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI
if args.MMS:
  MMS = args.MMS
  MME = MMS
if args.MME:
  MME = args.MME

varnm = ifld

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

dirspear = pthseas['ALL']['dirspear_anls']

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
TM       = []
dIspear  = []
dInsidc  = []
for YR in range(YRS,YRE+1):
  dclim_nsidc = msisrlx.read_NSIDC_iconc_clim_interp(YR)   # clim by calend months
  for MM in cal_months:
    itime = MM-1
    # Find model init. year for given calendar month:
    dnmb = mtime.datenum([YR,MM,15])
    YRI = manseas.yr_init_fcst_from_datenum(dnmb, MMI)
    if YRI < 1993:
      continue
    dirnsidc = pthseas['ALL']['dirnsidc_intrp'].format(YR=YR)  #NSIDC interpolated monthly iconc
    # SPEAR interpolated mean ensmb anomalies:
    flanom = f'spear_siconc_monanom_ensmean_{YRI}{MMI:02d}.nc'
    dflanom = os.path.join(dirspear,flanom)
    print(f'Loading SPEAR ens.mean  siconc anomalies --> {dflanom}')
    ds_spear = xarray.open_dataset(dflanom)
    mfcast  = ds_spear['months'].data
    ifcst  = np.where(mfcast==MM)[0][0]  # lead fcast = ifcst+1
    print(f'SPEAR: {dirspear}')
    print(f'Processing  {YR}/{MM} f/cast:  MMI={MMI} YRI={YRI} fcast lead mo={ifcst+1}')
    CIspear_anom = ds_spear['siconc_anom'].isel(months=ifcst).data

    # NSIDC interpolated
    pthnsidc = pthseas['ALL']['dirnsidc_intrp'].format(YR=YR)
    flnsidc = f'NSIDC_iconc_mnth_interpNEP{jdm}x{idm}_{YR}.nc'
    ds_nsidc = xarray.open_dataset(os.path.join(pthnsidc,flnsidc))
    CInsidc = ds_nsidc['ice_conc'].isel(time=itime).data
    Clim_nsidc = dclim_nsidc['ice_conc'].isel(months=itime).data

    # Anomalies:
    dIspear.append(CIspear_anom)
    dInsidc.append(CInsidc-Clim_nsidc)

    TM.append(dnmb)

dIspear = np.array(dIspear)
dInsidc = np.array(dInsidc)

TM = np.array(TM)
DV = mtime.datevec1D(TM, fHR=False)
DV = np.array(DV).transpose()

nrc = len(dIspear)
#nyr = int(nrc/len(cal_months))
#assert(nyr*12==nrc), f'Check record length dIspear={nrc}'

# Test pnt:
# Ber. Str
#i0=227
#j0=705
# W. Ber. Sea
i0=165
j0=691
a1 = dIspear[:,j0,i0]
a2 = dInsidc[:,j0,i0]
# Plot ACF
#from statsmodels.graphics.tsaplots import plot_acf
#aplot_acf(ts, lags=40)  # You can adjust lags


def corr_stat(dIspear, dInsidc):
  nrc = len(dIspear)
  # Corr:
  sgm_sp = np.std(dIspear, axis=0)
  sgm_ns = np.std(dInsidc, axis=0)
  mean_sp = np.mean(dIspear, axis=0)
  mean_ns = np.mean(dInsidc, axis=0)

  sgm_sp = np.where(sgm_sp < 1e-10, 1.e-10, sgm_sp)
  sgm_ns = np.where(sgm_ns < 1e-10, 1.e-10, sgm_ns)

  A3d = dIspear.copy()*0.
  for kk in range(nrc):
    Asp = dIspear[kk,:,:]
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
  DW_STAT = durbin_watson(dIspear, axis=0)

  return RR, PVAL, DW_STAT

RR, PVAL, DW_STAT = corr_stat(dIspear,dInsidc)

DW_STAT = np.where(np.isnan(RR), np.nan, DW_STAT)
DW_STAT = np.where(RR==0, 2., DW_STAT)

# Remove auto-correlation:
# t.s. are weekly correlated, this does not change much
#dIspear_ind = np.diff(dIspear, axis=0)
#dInsidc_ind = np.diff(dInsidc, axis=0)
#RR, PVAL, DW_STAT = corr_stat(dIspear_ind,dInsidc_ind)


from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

clrmp_dlt = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBu_r')
clrmp_dlt.set_bad(color=[0.2,0.2,0.2])
rmin = -1.
rmax = 1.

sinfo = 'Corr. from monthly anomlies of sea ice conc. wrt to 5-yr climatologies\n'
sinfo = sinfo + f'SPEAR MI{MMI:02d} mean ensmb seas. f/casts vs NSIDC NRT v3 ice conc'

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

sttl = f'Correlation & significance {alf:.2f} SPEAR MI:{MMI:02d} mnensmb vs NSIDC iconc anom\n'
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

btx = 'corr_iconc_SPEAR_NSIDC.py'
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

  sttl2 = f'Durbin-Watson a/corr Test SPEAR iconc anom {YRS}-{YRE} {MMS}-{MME}\n'
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

  





