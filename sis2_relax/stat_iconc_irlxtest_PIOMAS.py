"""
  Calc RMSE, ice extent and ice area, ice volume 
  from ice relaxation test runs and PIOMAS 
  monthly fields

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
#parser.add_argument("--regn", help="Region: NEP or ARC", type=str)
parser.add_argument("--yr", help="year to plot, default 2001", type=int)
parser.add_argument("--intrp", help=" =1: interp PIOMAS mnth to daily for better accur., default=1", \
                    type=int)
parser.add_argument("--nocntr", help=" =1: do not show the control no irlx run, default=0", \
                    type=int)
args = parser.parse_args()

# Test runs were performed for only 1 year
#regn = args.regn if args.regn else None
YRS = args.yr if args.yr else 2001
interp = args.intrp if args.intrp else 1
nocntr = args.nocntr if args.nocntr else 0

interp_mnthly = interp > 0  # for more accurate comparison, do time interpolation of PIOMAS 
                      # to get mnthly mean values, similar to how it is done in SIS2
                      # when deriving iconc ithkn for day=d0 from PIOMAS target fields
skip_noirlx = nocntr == 1  # Do not show control run with no irlx


#use_mnth = not (args.fday and args.fday > 0)
use_mnth = True
plt_rgn = False # Show Arc and Ber regions

pthrlx  = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax'

ifld = 'siconc' # partial area only
varnm = ifld

# Number of test runs:
Nexpts = 5 # Change to 2 to plot no irxl with and without ice ridging

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

DX, DY = mmom6.dx_dy(hlon, hlat)
Acell  = DX*DY*1.e-6  # km2

# Bering Sea - Chukchi Sea :
hsh = -5000.
LMsk = np.where((HH>=hsh) & (HH<0), 1, 0)
# Mask out southern lats:
LMsk = np.where(hlat<55.,0,LMsk)
LMsk[:567,:] = 0
LMsk[:,:39] = 0
LMsk[:595,177:] = 0
LMsk[:579,:129] = 0
LMsk[:575,:143] = 0
#LMsk[748:,:143] = 0

# Remove near-boundary points:
LMsk[810:,:] = 0
LMsk[:,338:] = 0

# 
# Mask for Bering Sea
# Bounded by the Bering Strait 
BMsk = LMsk.copy()
BMsk = np.where(hlat>66,0,BMsk)
JB,IB = np.where(BMsk==1)
# Mask for the Arctic Oc. part of the domain:
# Ber. Str. + S. Chukchi Shelf
AMsk = LMsk.copy()
AMsk = np.where(BMsk==1, 0, AMsk)
AMsk[:,:192] = 0
JA,IA = np.where(AMsk==1)

Aber = Acell[JB,IB]
Aarc = Acell[JA,IA]

#mcal = np.arange(mmi,mmi+12)
#mcal = np.where(mcal>12, mcal-12, mcal)
#d = abs(mcal-mm0)
#itime = np.argmin(d)

RMSE_Ber_ai = np.zeros((12,Nexpts))
RMSE_Arc_ai = np.zeros((12,Nexpts))
RMSE_Ber_hi = np.zeros((12,Nexpts))
RMSE_Arc_hi = np.zeros((12,Nexpts))
BIAS_Ber_ai = np.zeros((12,Nexpts))
BIAS_Arc_ai = np.zeros((12,Nexpts))
BIAS_Ber_hi = np.zeros((12,Nexpts))
BIAS_Arc_hi = np.zeros((12,Nexpts))
IAREA_Ber_nep = np.zeros((12,Nexpts))  # ice area Bering Sea, NEP test run
IAREA_Arc_nep = np.zeros((12,Nexpts))
IVOL_Ber_nep  = np.zeros((12,Nexpts))  # ice vol Ber. Sea, NEP test run
IVOL_Arc_nep  = np.zeros((12,Nexpts))
IAREA_Ber_pms = np.zeros((12))
IAREA_Arc_pms = np.zeros((12))
IVOL_Ber_pms  = np.zeros((12))
IVOL_Arc_pms  = np.zeros((12))
TM = []

ECOLR = np.array([[0.,0.4,0.9],
                  [0.9,0.5,0],
                  [0.,0.9,0.7],
                  [1.,0.9,0],
                  [0.8,0.,0.5]])

if Nexpts == 2:
  ECOLR[1,:] = np.array([1.,0.4,0])

itot = -1
for kexpt in range(Nexpts):
  expt_nmb = kexpt+1
  # Modified for no irlx with no ridging
  if Nexpts == 2 and kexpt == 1:
    expt_nmb = 11

  pthtest = f'/archive/Dmitry.Dukhovskoy/fre/NEP/test_ice_relax/NEPphys_expt{expt_nmb:02d}/{YRS}-01'
    
  if skip_noirlx and kexpt==0:
    print('Skipping noirlx control run')
    RMSE_Ber_ai[:,kexpt] = np.nan
    RMSE_Arc_ai[:,kexpt] = np.nan
    RMSE_Ber_hi[:,kexpt] = np.nan
    RMSE_Arc_hi[:,kexpt] = np.nan
    BIAS_Ber_ai[:,kexpt] = np.nan
    BIAS_Arc_ai[:,kexpt] = np.nan
    BIAS_Ber_hi[:,kexpt] = np.nan
    BIAS_Arc_hi[:,kexpt] = np.nan
    IAREA_Ber_nep[:,kexpt] = np.nan  # ice area Bering Sea, NEP test run
    IAREA_Arc_nep[:,kexpt] = np.nan
    IVOL_Ber_nep [:,kexpt] = np.nan  # ice vol Ber. Sea, NEP test run
    IVOL_Arc_nep [:,kexpt] = np.nan

    continue

  itot += 1
  YR = YRS
  for MM in range(1,13):
    print(f'Reading expt={kexpt+1} {YR}/{MM}')
    itime = MM-1
    dnmb0 = mtime.datenum([YR,MM,15])

    # irlx test experiments:
    prfx=''
    CInep = msisrlx.read_sis2_testrun(dnmb0, pthtest, prfx, 'iconc', use_mnth)
    HInep = msisrlx.read_sis2_testrun(dnmb0, pthtest, prfx, 'ithkn', use_mnth)

    # PIOMAS rlx fields (interpoalted onto NEP):
    YR1 = YR
    YR2 = YR+1
    flout = f'PIOMASv21_ithkn_iconc_{YR1}_{YR2}_monthly.nc'
    dfpiomas = os.path.join(pthrlx, flout)
    if interp_mnthly:
      CIpms = msisrlx.mnthly_PIOMAS_linear_daily(dfpiomas,dnmb0,'iconc')
      HIpms = msisrlx.mnthly_PIOMAS_linear_daily(dfpiomas,dnmb0,'ithkn')
    else:
      CIpms,_ = msisrlx.read_relax_piomas(dnmb0, pthrlx, 'iconc') 
      HIpms,_ = msisrlx.read_relax_piomas(dnmb0, pthrlx, 'ithkn') 

    # RMSE ice conc:
    Rsq = (CInep - CIpms)**2
    #nB = np.count_nonzero(~np.isnan(Rsq[(JB, IB)])) # this counts 0 and non-0 after np.isnan check
    nB = np.sum(~np.isnan(Rsq[JB, IB]))
    nA = np.sum(~np.isnan(Rsq[JA, IA]))
    rmseB_ai = np.sqrt(np.nansum(Rsq[JB,IB])/nB)
    rmseA_ai = np.sqrt(np.nansum(Rsq[JA,IA])/nA)
    biasB_ai = np.nanmean(CInep[JB,IB] - CIpms[JB,IB])
    biasA_ai = np.nanmean(CInep[JA,IA] - CIpms[JA,IA])

    # RMSE ice thickness:
    Rsq = (HInep - HIpms)**2
    nB = np.sum(~np.isnan(Rsq[JB, IB]))
    nA = np.sum(~np.isnan(Rsq[JA, IA]))
    rmseB_hi = np.sqrt(np.nansum(Rsq[JB,IB])/nB)
    rmseA_hi = np.sqrt(np.nansum(Rsq[JA,IA])/nA)
    biasB_hi = np.nanmean(HInep[JB,IB] - HIpms[JB,IB])
    biasA_hi = np.nanmean(HInep[JA,IA] - HIpms[JA,IA])
    
    # Ice Extent Ber. Sea:
    Cber_pms = CIpms[JB,IB]
    Cber_nep = CInep[JB,IB]
    IAreaB_nep = np.nansum(Cber_nep*Aber)  # ice area, km2, Ber. Sea reg
    IAreaB_pms = np.nansum(Cber_pms*Aber)

    # Ice Extent Arctic part:
    Carc_pms = CIpms[JA,IA]
    Carc_nep = CInep[JA,IA]
    IAreaA_nep = np.nansum(Carc_nep*Aarc)  # ice area, km2, Ber. Sea reg
    IAreaA_pms = np.nansum(Carc_pms*Aarc)

    # Ice volume Ber. Sea:
    VolIce_nep  = HInep*CInep*Acell*1e-3   # km3, HI - m, Acell - km
    VolB_nep = np.nansum(VolIce_nep[JB,IB])
    VolA_nep = np.nansum(VolIce_nep[JA,IA])

    VolIce_pms  = HIpms*CIpms*Acell*1e-3   # km3, HI - m, Acell - km
    VolB_pms = np.nansum(VolIce_pms[JB,IB])
    VolA_pms = np.nansum(VolIce_pms[JA,IA])

    # Register:
    # RMSE and biases:
    imm = MM-1
    RMSE_Ber_ai[imm,kexpt] = rmseB_ai
    RMSE_Arc_ai[imm,kexpt] = rmseA_ai
    RMSE_Ber_hi[imm,kexpt] = rmseB_hi
    RMSE_Arc_hi[imm,kexpt] = rmseA_hi

    BIAS_Ber_ai[imm,kexpt] = biasB_ai
    BIAS_Arc_ai[imm,kexpt] = biasA_ai
    BIAS_Ber_hi[imm,kexpt] = biasB_hi
    BIAS_Arc_hi[imm,kexpt] = biasA_hi

    # Ice area:
    IAREA_Ber_nep[imm,kexpt] = IAreaB_nep # Ber S. ice area, nep
    IAREA_Arc_nep[imm,kexpt] = IAreaA_nep # Arctic ice area, nep
    if itot == 0:
      IAREA_Ber_pms[imm] = IAreaB_pms # -"- -"- -"- , pms
      IAREA_Arc_pms[imm] = IAreaA_pms # -"- -"- -"- , pms

    # Ice volume:
    IVOL_Ber_nep[imm,kexpt] = VolB_nep
    IVOL_Arc_nep[imm,kexpt] = VolA_nep
    if itot == 0:
      IVOL_Ber_pms[imm] = VolB_pms
      IVOL_Arc_pms[imm] = VolA_pms

    # Time
    TM.append(dnmb0)

IAREA_Ber_nep = IAREA_Ber_nep*1e-3  # 1e3 km2
IAREA_Ber_pms = IAREA_Ber_pms*1e-3  # 1e3 km2
IAREA_Arc_nep = IAREA_Arc_nep*1e-3  # 1e3 km2
IAREA_Arc_pms = IAREA_Arc_pms*1e-3  # 1e3 km2

TM = np.array(TM)
DV = mtime.datevec1D(TM, fHR=False)
DV = np.array(DV).transpose()


pms_clr = [0,0.,0.4]


plt.ion()
def plot_ice_stat(fgnmb, iarea_nep, iarea_pms, ivol_nep, ivol_pms, 
                  rmse_ai, rmse_hi, bias_ai, bias_hi, regn):

  time_m = np.arange(1,13)

  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  # Ice area:
  ax1 = plt.axes([0.06, 0.72, 0.4, 0.25])
  for ik in range(Nexpts):
    clr1 = ECOLR[ik]
    ax1.plot(time_m,iarea_nep[:,ik], '-', linewidth=2, color=clr1)
  ax1.plot(time_m,iarea_pms, '-', linewidth=2, color=pms_clr)

  ax1.set_xticks(time_m)
  ax1.grid('on')
  ax1.set_xlim([0.5,12.5])
  sttl = f'IceArea x10^3 km2, {regn}, irlx tests'
  ax1.set_title(sttl)

  # Ice volume:
  ax2 = plt.axes([0.53, 0.72, 0.4, 0.25])
  for ik in range(Nexpts):
    clr1 = ECOLR[ik]
    ax2.plot(time_m,ivol_nep[:,ik], '-', linewidth=2, color=clr1)
  ax2.plot(time_m,ivol_pms, '-', linewidth=2, color=pms_clr)

  ax2.set_xticks(time_m)
  ax2.grid('on')
  ax2.set_xlim([0.5,12.5])
  sttl2 = f'IceVol km3, {regn}'
  ax2.set_title(sttl2)

  # Rmse ice conc
  ax3 = plt.axes([0.06, 0.41, 0.4, 0.25])
  for ik in range(Nexpts):
    clr1 = ECOLR[ik]
    ax3.plot(time_m, rmse_ai[:,ik], '-', linewidth=2, color=clr1)

  ax3.set_xticks(time_m)
  ax3.grid('on')
  ax3.set_xlim([0.5,12.5])
  sttl3 = f'RMSE_ai, {regn}'
  ax3.set_title(sttl3)

  # Rmse ice thkn
  ax4 = plt.axes([0.53, 0.41, 0.4, 0.25])
  for ik in range(Nexpts):
    clr1 = ECOLR[ik]
    ax4.plot(time_m, rmse_hi[:,ik], '-', linewidth=2, color=clr1)

  ax4.set_xticks(time_m)
  ax4.grid('on')
  ax4.set_xlim([0.5,12.5])
  sttl4 = f'RMSE_hi (m), {regn}'
  ax4.set_title(sttl4)

  # Bias iconc
  ax5 = plt.axes([0.06, 0.1, 0.4, 0.25])
  for ik in range(Nexpts):
    clr1 = ECOLR[ik]
    ax5.plot(time_m, bias_ai[:,ik], '-', linewidth=2, color=clr1)

  ylm = np.nanmax(abs(bias_ai))
  if 0.01 <= ylm < 1:
    ylm = np.ceil(ylm/0.1)*0.1
  elif 0.001 <= ylm < 0.1:
    ylm = np.ceil(ylm/0.01)*0.01
  elif 1<= ylm < 10:
    ylm = np.floor(ylm/0.1)*0.1
  else:
    ylm = np.ceil(ylm)

  ax5.set_xticks(time_m)
  ax5.grid('on')
  ax5.set_xlim([0.5,12.5])
  #ax5.set_ylim([-ylm,ylm])
  sttl5 = f'bias_ai (m), {regn}'
  ax5.set_title(sttl5)
 
  # Bias ice thickness    
  ax6 = plt.axes([0.53, 0.1, 0.4, 0.25])
  for ik in range(Nexpts):
    clr1 = ECOLR[ik]
    ax6.plot(time_m, bias_hi[:,ik], '-', linewidth=2, color=clr1)

  ylm = np.nanmax(abs(bias_hi))
  if 0.01 <= ylm < 1:
    ylm = np.ceil(ylm/0.1)*0.1
  elif 0.001 <= ylm < 0.1:
    ylm = np.ceil(ylm/0.01)*0.01
  elif 1<= ylm < 10:
    ylm = np.ceil(ylm/0.1)*0.1
  else:
    ylm = np.ceil(ylm)
    
  ax6.set_xticks(time_m)
  ax6.grid('on')
  ax6.set_xlim([0.5,12.5])
  #ax6.set_ylim([-ylm,ylm])
  sttl6 = f'bias_hi (m), {regn}'
  ax6.set_title(sttl6)

  # Legend:
  ax7 = plt.axes([0.6,0.01,0.39,0.08])
  x0 = 0.
  y0 = 0.35
  dy = 0.1
  xl = 0.12
  nlns = 3
  for ik in range(Nexpts):
    clr1 = ECOLR[ik]
    expt_nmb = ik+1
    if expt_nmb <= nlns:
      xS = x0
    else:
      xS = x0 + 0.5
      yS = y0 - dy*(ik-3)
      
    yS = y0 - dy*(ik%nlns)
    xE = xS+xl
    xT = xE+0.1*xl
    ax7.plot([xS,xE],[yS, yS],'-', linewidth=2, color=clr1)
    ax7.text(xT,yS,f'expt{expt_nmb:02d}', va='center')
  ik += 1
  yS = y0 - dy*(ik%nlns)
  ax7.plot([xS,xE],[yS, yS],'-', linewidth=2, color=pms_clr)
  ax7.text(xT,yS,f'PIOMAS',va='center')
  ax7.set_xlim([0,0.8])
  ax7.set_ylim([0.12,0.5])
  ax7.axis('off')
 
  btx = 'stat_iconc_irlxtest_PIOMAS.py'
  bottom_text(btx, pos=[0.02,0.01])

  return fig1, ax1, ax2, ax3, ax4, ax5, ax6, ax7

# plt.figure(fig1)

clr_ber = [0,0.4,0.8]
clr2_ber = [0.7,0.8,1]
clr_arc = [0.,0.7,0.2]
clr2_arc = [0.7,1,0.9]
time_mnth = np.arange(1,13)

f_debug = False
if f_debug:
  iarea_nep = IAREA_Ber_nep
  iarea_pms = IAREA_Ber_pms
  ivol_nep  = IVOL_Ber_nep
  ivol_pms  = IVOL_Ber_pms
  rmse_ai   = RMSE_Ber_ai
  rmse_hi   = RMSE_Ber_hi
  bias_ai   = BIAS_Ber_ai
  bias_hi   = BIAS_Ber_hi
  regn      = 'BerS'

# Plot time series and monthly stat for RMSE and bias, Bering Sea:
fgnmb=1
time_yrs = DV[:,0]+(DV[:,1]-1)/12
fig1,ax11,ax12,ax13,ax14,ax15, ax16, ax17 = plot_ice_stat(fgnmb, IAREA_Ber_nep, IAREA_Ber_pms, \
                                         IVOL_Ber_nep, IVOL_Ber_pms, \
                                         RMSE_Ber_ai, RMSE_Ber_hi, BIAS_Ber_ai, BIAS_Ber_hi, 'BerS')

fig2,ax21,ax22,ax23,ax24,ax25, ax26, ax27 = plot_ice_stat(2, IAREA_Arc_nep, IAREA_Arc_pms, \
                                         IVOL_Arc_nep, IVOL_Arc_pms, \
                                         RMSE_Arc_ai, RMSE_Arc_hi, BIAS_Arc_ai, BIAS_Arc_hi, 'ArcOc')

if plt_rgn:
  fig3 = plt.figure(3,figsize=(9,8))
  plt.clf()
  # Ice area:
  ax31 = plt.axes([0.1, 0.1, 0.8, 0.8])
  LMSK = np.where(HH<0,1,0)
  cmp_lmsk = mclrmps.colormap_landmask()
  ax31.pcolormesh(LMSK, cmap=cmp_lmsk)
  ax31.axis('scaled')
  ax31.set_ylim([550,816])
  ax31.plot(IA,JA,'.',color=[0.8,0.8,0.8])
  ax31.plot(IB,JB,'.',color=[0.4,0.4,0.4])

  loncntrs = [x for x in range(140,360,10)]
  latcntrs = [x for x in range(40,89,10)]
  ax31.contour(hlon, levels=loncntrs, linstyles='-', linewidths=1, colors=[[0.9,0.9,0.9]])
  ax31.contour(hlat, levels=latcntrs, linstyles='-', linewidths=1, colors=[[0.9,0.9,0.9]])

  sttl = 'Ber and Arc regions for stat anls'
  ax31.set_title(sttl)

  btx = 'stat_iconc_irlxtest_PIOMAS.py'
  bottom_text(btx, pos=[0.1,0.08])
 


  

