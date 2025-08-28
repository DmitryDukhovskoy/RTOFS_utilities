"""
  Calc RMSE, ice extent and ice area, ice volume 
  from ice relaxation test runs and PIOMAS 
  monthly fields

  ARC12 domain

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
parser.add_argument("--yr", help="year to plot, default 2001", type=int)
parser.add_argument("--intrp", help=">1: interp PIOMAS mnth to daily for better accur., default=0", type=int)
args = parser.parse_args()

# Test runs were performed for only 1 year
YRS = args.yr if args.yr else 1995
interp = args.intrp if args.intrp else 0  
interp_mnthly = interp > 0  # for more accurate comparison, do time interpolation of PIOMAS 
                      # to get mnthly mean values, similar to how it is done in SIS2
                      # when deriving iconc ithkn for day=d0 from PIOMAS target fields
#use_mnth = not (args.fday and args.fday > 0)
regn = 'ARC'
use_mnth = True
plt_rgn = False # Show Arc and Ber regions

pthrlx  = '/work/Dmitry.Dukhovskoy/ARC12/irlx'

# Arctic domain:
ptharc  = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'
dflarc  = os.path.join(ptharc,'ocean_hgrid.nc')
dfltopo = os.path.join(ptharc,'ocean_topog.nc')
hlon, hlat = mmom6.read_mom6grid(dflarc, grdpnt='hgrid')

ds_topo = xarray.open_dataset(dfltopo)
HH = -(ds_topo['depth'].data)
jdm, idm = HH.shape

assert HH[300,200] < 0., f'Check sign of topography, ocean pnts should be < 0'

DX, DY = mmom6.dx_dy(hlon, hlat)
Acell  = DX*DY*1.e-6  # km2

# Perform stat analysis only north of 60N, where relaxation is applied:
hsh = -5000.
LMsk = np.where((HH>=hsh) & (HH<0), 1, 0)
LMsk = np.where(hlat<60., 0, LMsk)

JA,IA = np.where(LMsk == 1)
Aarc = Acell[JA,IA]


Nexpts = 5
RMSE_Arc_ai = np.zeros((12,Nexpts))
RMSE_Arc_hi = np.zeros((12,Nexpts))
BIAS_Arc_ai = np.zeros((12,Nexpts))
BIAS_Arc_hi = np.zeros((12,Nexpts))
IAREA_arc = np.zeros((12,Nexpts))
IVOL_arc  = np.zeros((12,Nexpts))
IAREA_pms = np.zeros((12))
IVOL_pms  = np.zeros((12))
TM = []

for kexpt in range(5):
  expt_nmb = kexpt+1
  pthtest = f'/archive/Dmitry.Dukhovskoy/fre/ARC12/test_ice_relax/ARCphys_expt{expt_nmb:02d}/{YRS}-01'

  YR = YRS
  for MM in range(1,13):
    print(f'Reading expt={kexpt+1} {YR}/{MM}')
    itime = MM-1
    dnmb0 = mtime.datenum([YR,MM,15])

    # irlx test experiments:
    prfx=''
    CIarc = msisrlx.read_sis2_testrun(dnmb0, pthtest, prfx, 'iconc', use_mnth)
    HIarc = msisrlx.read_sis2_testrun(dnmb0, pthtest, prfx, 'ithkn', use_mnth)

    # PIOMAS rlx fields (interpoalted onto NEP):
    YR1 = YR
    YR2 = YR+1
    flout = f'PIOMASv21_ARC12_ithkn_iconc_{YR1}_{YR2}_monthly.nc'
    dfpiomas = os.path.join(pthrlx, flout)
    if interp_mnthly:
      CIpms = msisrlx.mnthly_PIOMAS_linear_daily(dfpiomas,dnmb0,'iconc')
      HIpms = msisrlx.mnthly_PIOMAS_linear_daily(dfpiomas,dnmb0,'ithkn')
    else:
      CIpms,_ = msisrlx.read_relax_piomas(dnmb0, pthrlx, 'iconc') 
      HIpms,_ = msisrlx.read_relax_piomas(dnmb0, pthrlx, 'ithkn') 

    # RMSE ice conc:
    Rsq = (CIarc - CIpms)**2
    #nB = np.count_nonzero(~np.isnan(Rsq[(JB, IB)])) # this counts 0 and non-0 after np.isnan check
    nA = np.sum(~np.isnan(Rsq[JA, IA]))
    rmseA_ai = np.sqrt(np.nansum(Rsq[JA,IA])/nA)
    biasA_ai = np.nanmean(CIarc[JA,IA] - CIpms[JA,IA])

    # RMSE ice thickness:
    Rsq = (HIarc - HIpms)**2
    nA = np.sum(~np.isnan(Rsq[JA, IA]))
    rmseA_hi = np.sqrt(np.nansum(Rsq[JA,IA])/nA)
    biasA_hi = np.nanmean(HIarc[JA,IA] - HIpms[JA,IA])
    
    # Ice Extent Arctic part:
    Cpms = CIpms[JA,IA]
    Carc = CIarc[JA,IA]
    IAreaA_arc = np.nansum(Carc*Aarc)  # ice area, km2, Ber. Sea reg
    IAreaA_pms = np.nansum(Cpms*Aarc)

    # Ice volume:
    VolIce_arc  = HIarc*CIarc*Acell*1e-3   # km3, HI - m, Acell - km
    VolA_arc = np.nansum(VolIce_arc[JA,IA])

    VolIce_pms  = HIpms*CIpms*Acell*1e-3   # km3, HI - m, Acell - km
    VolA_pms = np.nansum(VolIce_pms[JA,IA])

    # Register:
    # RMSE and biases:
    imm = MM-1
    RMSE_Arc_ai[imm,kexpt] = rmseA_ai
    RMSE_Arc_hi[imm,kexpt] = rmseA_hi

    BIAS_Arc_ai[imm,kexpt] = biasA_ai
    BIAS_Arc_hi[imm,kexpt] = biasA_hi

    # Ice area:
    IAREA_arc[imm,kexpt] = IAreaA_arc # Arctic ice area, nep
    if kexpt == 0:
      IAREA_pms[imm] = IAreaA_pms # -"- -"- -"- , pms

    # Ice volume:
    IVOL_arc[imm,kexpt] = VolA_arc
    if kexpt == 0:
      IVOL_pms[imm] = VolA_pms

    # Time
    TM.append(dnmb0)

IAREA_arc = IAREA_arc*1e-3  # 1e3 km2
IAREA_pms = IAREA_pms*1e-3  # 1e3 km2

TM = np.array(TM)
DV = mtime.datevec1D(TM, fHR=False)
DV = np.array(DV).transpose()

ECOLR = [[0.,0.4,0.9],
         [0.9,0.5,0],
         [0.,0.9,0.7],
         [1.,0.9,0],
         [0.8,0.,0.5]]

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

  ylm = np.max(abs(bias_ai))
  if 0.01 <= ylm < 1:
    ylm = np.ceil(ylm/0.1)*0.1
  elif 0.001 <= ylm < 0.1:
    ylm = np.ceil(ylm/0.01)*0.01
  elif 1<= ylm < 100:
    ylm = np.ceil(ylm)

  ax5.set_xticks(time_m)
  ax5.grid('on')
  ax5.set_xlim([0.5,12.5])
  ax5.set_ylim([-ylm,ylm])
  sttl5 = f'bias_ai (m), {regn}'
  ax5.set_title(sttl5)
 
  # Bias ice thickness    
  ax6 = plt.axes([0.53, 0.1, 0.4, 0.25])
  for ik in range(Nexpts):
    clr1 = ECOLR[ik]
    ax6.plot(time_m, bias_hi[:,ik], '-', linewidth=2, color=clr1)

  ylm = np.max(abs(bias_hi))
  if 0.01 <= ylm < 1:
    ylm = np.ceil(ylm/0.1)*0.1
  elif 0.001 <= ylm < 0.1:
    ylm = np.ceil(ylm/0.01)*0.01
  elif 1<= ylm < 100:
    ylm = np.ceil(ylm)
    
  ax6.set_xticks(time_m)
  ax6.grid('on')
  ax6.set_xlim([0.5,12.5])
  ax6.set_ylim([-ylm,ylm])
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

fig1,ax11,ax12,ax13,ax14,ax15, ax16, ax17 = plot_ice_stat(1, IAREA_arc, IAREA_pms, \
                                         IVOL_arc, IVOL_pms, \
                                         RMSE_Arc_ai, RMSE_Arc_hi, BIAS_Arc_ai, BIAS_Arc_hi, 'ArcOc')

if plt_regn:
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
 


  

