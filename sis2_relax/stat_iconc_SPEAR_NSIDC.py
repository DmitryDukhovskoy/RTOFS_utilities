"""
  Calc RMSE, ice extent and ice area 
 
  usage: plot_SPEAR_ice_month_stere.py --YRI=1993 --MMI=1

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
parser.add_argument("--YRS", help="init year of SPEAR f/cast: 1993, ..., 2020", type=int)
parser.add_argument("--YRE", help="init year of SPEAR f/cast: 1993, ..., 2020", type=int)
parser.add_argument("--MMI", help="init month of SPEAR f/cast: 1, ..., 12", type=int)
parser.add_argument("--ensmb", help="ensemble number, 1,..., 15", type=int)
args = parser.parse_args()

f_cntrobs = True   # Plot observation-derived ice edge
# Years in the relax file also used in the rlx file name:
YRS = 1993  # init yr
YRE = 2020
MMI = 1     # init month
ifld = 'siconc' # partial area only
ens_nmb = 1  # SPEAR ensemble #


if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI
if args.ensmb:
  ens_nmb = args.ensmb

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

# Select correct init year for given YR/MM
TM       = []
RMSE_Ber = []
RMSE_Arc = []
BIAS_Ber = []
BIAS_Arc = []
IAS_Ber  = []
IAN_Ber  = []
IAS_Arc  = []
IAN_Arc  = []
for YR in range(YRS,YRE+1):
  for MM in range(1,13):
    itime = MM-1
    dnmb = mtime.datenum([YR,MM,15])
    YRI = manseas.yr_init_fcst_from_datenum(dnmb, MMI)
    dirnsidc = pthseas['ALL']['dirnsidc_intrp'].format(YR=YRI)  #NSIDC interpolated monthly iconc fields
    # SPEAR interpolated
    flout = f'spear_interpNEP_siconc_mnth_{YR}{MMI:02d}.nc'
    dfspear = os.path.join(dirspear,flout)
    ds_spear = xarray.open_dataset(dfspear)
    mfcast  = ds_spear['months'].data
    ifcst  = np.where(mfcast==MM)[0][0]
    CIspear = ds_spear['siconc'].isel(nmonths=ifcst).data

    print(f'SPEAR: {dirnsidc}')
    print(f'Processing MMI={MMI} {YR}/{MM} f/cast mo={ifcst+1}')
    # NSIDC interpolated
    pthnsidc = pthseas['ALL']['dirnsidc_intrp'].format(YR=YR)
    flnsidc = f'NSIDC_iconc_mnth_interpNEP{jdm}x{idm}_{YR}.nc'
    ds_nsidc = xarray.open_dataset(os.path.join(pthnsidc,flnsidc))
    CInsidc = ds_nsidc['ice_conc'].isel(time=itime).data

    Rsq = (CIspear - CInsidc)**2
    nB = np.count_nonzero(~np.isnan(Rsq[JB, IB]))
    nA = np.count_nonzero(~np.isnan(Rsq[JA, IA]))
    rmseB = np.sqrt(np.nansum(Rsq[JB,IB])/nB)
    rmseA = np.sqrt(np.nansum(Rsq[JA,IA])/nA)
    biasB = np.nanmean(CIspear[JB,IB] - CInsidc[JB,IB])
    biasA = np.nanmean(CIspear[JA,IA] - CInsidc[JA,IA])

    # Ice Extent Ber. Sea:
    Cber_nsidc = CInsidc[JB,IB]
    Cber_spear = CIspear[JB,IB]
    Aber = Acell[JB,IB]
    knan = np.where(np.isnan(Cber_nsidc)) # mismatch near coasts, Aleutian islands
    Cber_spear[knan] = np.nan
    IAreaB_spear = np.nansum(Cber_spear*Aber)  # ice area, km2, Ber. Sea reg
    IAreaB_nsidc = np.nansum(Cber_nsidc*Aber)

    # Ice Extent Arctic part:
    Carc_nsidc = CInsidc[JA,IA]
    Carc_spear = CIspear[JA,IA]
    Aarc = Acell[JA,IA]
    knan = np.where(np.isnan(Carc_nsidc)) # mismatch near coasts, Aleutian islands
    Carc_spear[knan] = np.nan
    IAreaA_spear = np.nansum(Carc_spear*Aarc)  # ice area, km2, Ber. Sea reg
    IAreaA_nsidc = np.nansum(Carc_nsidc*Aarc)

    TM.append(dnmb)
    RMSE_Ber.append(rmseB)
    RMSE_Arc.append(rmseA)
    BIAS_Ber.append(biasB)
    BIAS_Arc.append(biasA)
    IAS_Ber.append(IAreaB_spear) # Ber S. ice area, spear
    IAN_Ber.append(IAreaB_nsidc) # -"- -"- -"- , nsidc
    IAS_Arc.append(IAreaA_spear) # Arctic ice area, spear
    IAN_Arc.append(IAreaA_nsidc) # -"- -"- -"- , nsidc

RMSE_Ber = np.array(RMSE_Ber)
RMSE_Arc = np.array(RMSE_Arc)
BIAS_Ber = np.array(BIAS_Ber)
BIAS_Arc = np.array(BIAS_Arc)
IAS_Ber  = np.array(IAS_Ber)*1e-3  # 1e3 km2
IAN_Ber  = np.array(IAN_Ber)*1e-3
IAS_Arc  = np.array(IAS_Arc)*1e-3
IAN_Arc  = np.array(IAN_Arc)*1e-3

TM = np.array(TM)
DV = mtime.datevec1D(TM, fHR=False)
DV = np.array(DV).transpose()

nyr = int(len(RMSE_Ber)/12)
assert(nyr*12==len(RMSE_Ber)), f'Check record length RMSE_Ber={len(RMSE_Ber)}'

R2dB  = RMSE_Ber.reshape(nyr,12)
R2dA  = RMSE_Arc.reshape(nyr,12)
B2dB  = BIAS_Ber.reshape(nyr,12)
B2dA  = BIAS_Arc.reshape(nyr,12)
I2dSB = IAS_Ber.reshape(nyr,12) # ice area Ber. Sea spear
I2dNB = IAN_Ber.reshape(nyr,12) # ice area Ber. Sea NSIDC
I2dSA = IAS_Arc.reshape(nyr,12) # ice area Arc, spear
I2dNA = IAN_Arc.reshape(nyr,12) 

# Bering Sea
R2dB_md = np.median(R2dB, axis=0)
R2dB_pu = np.percentile(R2dB, 90, axis=0)
R2dB_pl = np.percentile(R2dB, 10, axis=0)

B2dB_md = np.median(B2dB, axis=0)
B2dB_pu = np.percentile(B2dB, 90, axis=0)
B2dB_pl = np.percentile(B2dB, 10, axis=0)

# Arctic
R2dA_md = np.median(R2dA, axis=0)
R2dA_pu = np.percentile(R2dA, 90, axis=0)
R2dA_pl = np.percentile(R2dA, 10, axis=0)

B2dA_md = np.median(B2dA, axis=0)
B2dA_pu = np.percentile(B2dA, 90, axis=0)
B2dA_pl = np.percentile(B2dA, 10, axis=0)

# Ice area:
I2dSB_md = np.median(I2dSB, axis=0)
I2dSB_pu = np.percentile(I2dSB, 90, axis=0)
I2dSB_pl = np.percentile(I2dSB, 10, axis=0)
I2dNB_md = np.median(I2dNB, axis=0)
I2dNB_pu = np.percentile(I2dNB, 90, axis=0)
I2dNB_pl = np.percentile(I2dNB, 10, axis=0)

I2dSA_md = np.median(I2dSA, axis=0)
I2dSA_pu = np.percentile(I2dSA, 90, axis=0)
I2dSA_pl = np.percentile(I2dSA, 10, axis=0)
I2dNA_md = np.median(I2dNA, axis=0)
I2dNA_pu = np.percentile(I2dNA, 90, axis=0)
I2dNA_pl = np.percentile(I2dNA, 10, axis=0)

plt.ion()
def plot_rmse_bias(fgnmb,time_yrs, rmse, rmse_md,rmse_pu,rmse_pl,\
                   bias, bias_md, bias_pu, bias_pl,clr1,clr2,regn):
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  # Time series RMSE
  time_yrs = DV[:,0]+(DV[:,1]-1)/12
  ax1 = plt.axes([0.06, 0.56, 0.4, 0.38])
  ax1.plot(time_yrs,rmse, '-', linewidth=2, color=clr1)
  sttl = f'SPEAR vs NSIDC RMSE Ice Area \n {regn} {YRS}-{YRE}'
  ax1.set_title(sttl)

  # Monthly RMSE
  time_mnth = np.arange(1,13)
  ax2 = plt.axes([0.06,0.08,0.4,0.38])
  ax2.plot(time_mnth,rmse_md,'-',linewidth=2, color=clr1)
  ax2.plot(time_mnth,rmse_md, marker='o', markersize=7, color=clr1)
  ax2.plot(time_mnth,rmse_pu,'-',linewidth=1, color=clr2)
  ax2.plot(time_mnth,rmse_pl,'-',linewidth=1, color=clr2)
  ax2.set_xticks(time_mnth)
  ax2.grid('on')
  ax2.set_title(f'RMSE {regn}, Median & IDR')
  ax2.set_xlabel('Months')

  # Time Series bias:
  ax3 = plt.axes([0.56, 0.56, 0.4, 0.38])
  ax3.plot(time_yrs,bias, '-', linewidth=2, color=clr1)
  sttl = f'SPEAR vs NSIDC Bias Ice Area \n {regn} {YRS}-{YRE}'
  ax3.set_title(sttl)

  # Monthly Bias
  ax4 = plt.axes([0.56,0.08,0.4,0.38])
  ax4.plot(time_mnth,bias_md,'-',linewidth=2, color=clr1)
  ax4.plot(time_mnth,bias_md, marker='o', markersize=7, color=clr1)
  ax4.plot(time_mnth,bias_pu,'-',linewidth=1, color=clr2)
  ax4.plot(time_mnth,bias_pl,'-',linewidth=1, color=clr2)
  ax4.set_xticks(time_mnth)
  ax4.grid('on')
  ax4.set_title(f'Bias {regn}, Median & IDR')
  ax4.set_xlabel('Months')


  btx = 'stat_iconc_SPEAR_NSIDC.py'
  bottom_text(btx, pos=[0.02,0.01])

  return fig1, ax1, ax2, ax3, ax4

# plt.figure(fig1)

clr_ber = [0,0.4,0.8]
clr2_ber = [0.7,0.8,1]
clr_arc = [0.,0.7,0.2]
clr2_arc = [0.7,1,0.9]

# Plot time series and monthly stat for RMSE and bias, Bering Sea:
fgnmb=1
time_yrs = DV[:,0]+(DV[:,1]-1)/12
fig1,ax11,ax12,ax13,ax14 = plot_rmse_bias(fgnmb,time_yrs,RMSE_Ber,R2dB_md,R2dB_pu,R2dB_pl,\
               BIAS_Ber,B2dB_md,B2dB_pu,B2dB_pl,clr_ber,clr2_ber,"Bering Sea")

# RMSE and bias for Arct.
fig3,ax31,ax32,ax33,ax34 = plot_rmse_bias(3,time_yrs,RMSE_Arc,R2dA_md,R2dA_pu,R2dA_pl,\
               BIAS_Arc,B2dA_md,B2dA_pu,B2dA_pl,clr_arc,clr2_arc,"Arctic Ocean")

# Time Ser. of ice area:
# Ber. Sea
clr_spear = [0.,0.6,0.8]
clr2_spear = [0.8,0.9,1]
clr_nsidc = [0.8,0.3,0]
clr2_nsidc = [1,0.9,0.8]
fig2 = plt.figure(2,figsize=(9,8))
plt.clf()
# Time series 
ax21 = plt.axes([0.06, 0.56, 0.4, 0.38])
ln1, = ax21.plot(time_yrs,IAS_Ber, '-', linewidth=2, color=clr_spear, label='SPEAR')
ln2, = ax21.plot(time_yrs,IAN_Ber, '-', linewidth=1.6, color=clr_nsidc, label='NSIDC NRT')
sttl21 = f'SPEAR & NSIDC Ice Area*1.e3 km2 \n Bering Sea Reg. {YRS}-{YRE}'
ax21.set_title(sttl21)

ax22 = plt.axes([0.06,0.08,0.4,0.38])
# SPEAR
ax22.plot(time_mnth,I2dSB_md,'-',linewidth=2, color=clr_spear)
ax22.plot(time_mnth,I2dSB_md, marker='o', markersize=7, color=clr_spear)
ax22.plot(time_mnth,I2dSB_pu,'-',linewidth=1, color=clr2_spear)
ax22.plot(time_mnth,I2dSB_pl,'-',linewidth=1, color=clr2_spear)
# NDISC
ax22.plot(time_mnth,I2dNB_md,'-',linewidth=2, color=clr_nsidc)
ax22.plot(time_mnth,I2dNB_md, marker='o', markersize=7, color=clr_nsidc)
ax22.plot(time_mnth,I2dNB_pu,'-',linewidth=1, color=clr2_nsidc)
ax22.plot(time_mnth,I2dNB_pl,'-',linewidth=1, color=clr2_nsidc)
ax22.set_xticks(time_mnth)
ax22.grid('on')
ax22.set_title('Ice Area Ber. Sea, Median & IDR')
ax22.set_xlabel('Months')

# Arctic:
ax23 = plt.axes([0.56, 0.56, 0.4, 0.38])
ax23.plot(time_yrs,IAN_Arc, '-', linewidth=2, color=clr_nsidc, label='NSIDC NRT')
ax23.plot(time_yrs,IAS_Arc, '-', linewidth=1.6, color=clr_spear, label='SPEAR')
sttl23 = f'SPEAR & NSIDC Ice Area*1.e3 km2 \n Arctic {YRS}-{YRE}'
ax23.set_title(sttl23)

ax24 = plt.axes([0.56,0.08,0.4,0.38])
# SPEAR
ax24.plot(time_mnth,I2dSA_md,'-',linewidth=2, color=clr_spear)
ax24.plot(time_mnth,I2dSA_md, marker='o', markersize=7, color=clr_spear)
ax24.plot(time_mnth,I2dSA_pu,'-',linewidth=1, color=clr2_spear)
ax24.plot(time_mnth,I2dSA_pl,'-',linewidth=1, color=clr2_spear)
# NDISC
ax24.plot(time_mnth,I2dNA_md,'-',linewidth=2, color=clr_nsidc)
ax24.plot(time_mnth,I2dNA_md, marker='o', markersize=7, color=clr_nsidc)
ax24.plot(time_mnth,I2dNA_pu,'-',linewidth=1, color=clr2_nsidc)
ax24.plot(time_mnth,I2dNA_pl,'-',linewidth=1, color=clr2_nsidc)
ax24.set_xticks(time_mnth)
ax24.grid('on')
ax24.set_title('Ice Area Arctic, Median & IDR')
ax24.set_xlabel('Months')

# Legend
ax25 = plt.axes([0.48, 0.44, 0.1, 0.1])
lgd = plt.legend(handles=[ln1,ln2], loc='upper right')
ax25.axis('off')

btx = 'stat_iconc_SPEAR_NSIDC.py'
bottom_text(btx, pos=[0.02,0.01])




