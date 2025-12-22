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

yr_run = 2001
regn = 'NEP'

parser = argparse.ArgumentParser()
#parser.add_argument("--regn", help="Region: NEP or ARC", type=str)
parser.add_argument("--yr", help=f"year to plot, default {yr_run}", type=int)
parser.add_argument("--intrp", help=" =1: interp PIOMAS mnth to daily for better accur., default=1", \
                    type=int)
parser.add_argument("--warea", help="T: weight RMSE and bias by cell area, default=N", choices=["T","F"], type=str)
parser.add_argument(
    "--enmb",
    help="List of experiment numbers (e.g., 1 3 5 32 33)",
    type=int,
    nargs="+",             # <-- allows one or more integers and will generate a list
    required=True
)
args = parser.parse_args()

# Test runs were performed for only 1 year
#regn = args.regn if args.regn else None
YRS = args.yr if args.yr else yr_run
interp = args.intrp if args.intrp else 1
warea = args.warea if args.warea else "F"
ENMBS = args.enmb if args.enmb else None

# Number of test runs:
Nexpts = len(ENMBS)

interp_mnthly = interp > 0  # for more accurate comparison, do time interpolation of PIOMAS 
                      # to get mnthly mean values, similar to how it is done in SIS2
                      # when deriving iconc ithkn for day=d0 from PIOMAS target fields

use_mnth = True
plt_rgn = False # Show Arc and Ber regions
if warea.lower() == "t":
    wt_area = True
elif warea.lower() == "f":
    wt_area = False
else:
    raise Exception(f"Unknown option for warea: {warea}")


pthrlx  = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax'

ifld = 'siconc' # partial area only
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

ECOLR = msisrlx.irlx_tests_colors()

itot = -1
for kexpt in range(Nexpts):
  expt_nmb = ENMBS[kexpt]

  pthtest = f'/archive/Dmitry.Dukhovskoy/fre/NEP/test_ice_relax/NEPphys_expt{expt_nmb:02d}/{YRS}-01'
    
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
    maskB = ~np.isnan(Rsq[JB, IB])
    maskA = ~np.isnan(Rsq[JA, IA])
    if wt_area:
      # Do cell area weighted RMSE & Bias
      nB = np.sum(Acell[JB, IB][maskB])
      nA = np.sum(Acell[JA, IA][maskA])
      rmseB_ai = np.sqrt(np.nansum(Acell[JB, IB][maskB] * Rsq[JB, IB][maskB]) / nB)
      rmseA_ai = np.sqrt(np.nansum(Acell[JA, IA][maskA] * Rsq[JA, IA][maskA]) / nA)
      biasB_ai = np.nansum(Acell[JB, IB][maskB]* (CInep[JB,IB][maskB] - CIpms[JB,IB][maskB])) / nB
      biasA_ai = np.nansum(Acell[JA, IA][maskA]* (CInep[JA,IA][maskA] - CIpms[JA,IA][maskA])) / nA
    else:
      # No area weighting: 
      nB = np.sum(~np.isnan(Rsq[JB, IB]))
      nA = np.sum(~np.isnan(Rsq[JA, IA]))
      rmseB_ai = np.sqrt(np.nansum(Rsq[JB,IB])/nB)
      rmseA_ai = np.sqrt(np.nansum(Rsq[JA,IA])/nA)
      biasB_ai = np.nanmean(CInep[JB,IB] - CIpms[JB,IB])
      biasA_ai = np.nanmean(CInep[JA,IA] - CIpms[JA,IA])

    # RMSE ice thickness:
    Rsq = (HInep - HIpms)**2
    maskB = ~np.isnan(Rsq[JB, IB])
    maskA = ~np.isnan(Rsq[JA, IA])
    if wt_area:
      # Do cell area weighted RMSE & Bias
      nB = np.sum(Acell[JB, IB][maskB])
      nA = np.sum(Acell[JA, IA][maskA])
      rmseB_hi = np.sqrt(np.nansum(Acell[JB, IB][maskB] * Rsq[JB, IB][maskB]) / nB)
      rmseA_hi = np.sqrt(np.nansum(Acell[JA, IA][maskA] * Rsq[JA, IA][maskA]) / nA)
      biasB_hi = np.nansum(Acell[JB, IB][maskB]* (HInep[JB,IB][maskB] - HIpms[JB,IB][maskB])) / nB
      biasA_hi = np.nansum(Acell[JA, IA][maskA]* (HInep[JA,IA][maskA] - HIpms[JA,IA][maskA])) / nA

    else:
      # No area weigthing:
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

def plot_ax0(ax0, time_m, data0, data_pms, sttl, lcntrl):
  lnw0 = 2
  for ik in range(Nexpts):
    clr1 = ECOLR[ik]
    if ENMBS[ik] == 1 and lcntrl:
      lnw = lnw0
      lstl = '-'
      marker='o'
      mksz=4
    else:
      lnw = lnw0
      lstl = '-'
      marker=None
      mksz=None
    ax0.plot(time_m, data0[:,ik], 
             linestyle=lstl, linewidth=lnw, marker=marker, markersize=mksz, color=clr1)

  if data_pms is not None:
    ax0.plot(time_m, data_pms, linestyle='-', linewidth=lnw0, color=pms_clr)

  ax0.set_xticks(time_m)
  ax0.grid('on')
  ax0.set_xlim([time_m[0]-0.25, time_m[-1]+0.25])
  ax0.set_title(sttl)

  return ax0

def plot_ice_stat(fgnmb, iarea_nep, iarea_pms, ivol_nep, ivol_pms, 
                  rmse_ai, rmse_hi, bias_ai, bias_hi, subregn, lcntrl=False):
  """
    lcntrl - show control run with a thicker line
  """

  time_m = np.arange(1,13)

  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()

  # Ice area:
  ax1 = plt.axes([0.06, 0.72, 0.4, 0.25])
  lnw0 = 2
  sttl = f'IceArea x10^3 km2, {subregn}, irlx tests'
  ax1 = plot_ax0(ax1, time_m, iarea_nep, iarea_pms, sttl, lcntrl)

  # Ice volume:
  ax2 = plt.axes([0.53, 0.72, 0.4, 0.25])
  sttl2 = f'IceVol km3, {subregn}'
  ax2 = plot_ax0(ax2, time_m, ivol_nep, ivol_pms, sttl2, lcntrl)

  # Rmse ice conc
  ax3 = plt.axes([0.06, 0.41, 0.4, 0.25])
  sttl3 = f'RMSE_ai, {subregn}'
  ax3 = plot_ax0(ax3, time_m, rmse_ai, None, sttl3, lcntrl)

  # Rmse ice thkn
  ax4 = plt.axes([0.53, 0.41, 0.4, 0.25])
  sttl4 = f'RMSE_hi (m), {subregn}'
  ax4 = plot_ax0(ax4, time_m, rmse_hi, None, sttl4, lcntrl)

  # Bias iconc
  ax5 = plt.axes([0.06, 0.1, 0.4, 0.25])
  sttl5 = f'bias_ai (m), {subregn}'
  ax5 = plot_ax0(ax5, time_m, bias_ai, None, sttl5, lcntrl)  

  # Bias ice thickness    
  ax6 = plt.axes([0.53, 0.1, 0.4, 0.25])
  sttl6 = f'bias_hi (m), {subregn}'
  ax6 = plot_ax0(ax6, time_m, bias_hi, None, sttl6, lcntrl)

  # Legend:
  ax7 = plt.axes([0.05,0.02,0.8,0.08])
  x0 = 0.02
  y0 = 0.32
  dy = 0.1
  xl = 0.03
  nlns = 2
  dx = 0.15
  for ik in range(Nexpts):
    clr1 = ECOLR[ik]
    if ENMBS[ik] == 1 and lcntrl:
      lnw = lnw0
      lstl = '-'
      marker='o'
      mksz=4
    else:
      lnw = lnw0
      lstl = '-'
      marker=None
      mksz=None

    expt_nmb = ENMBS[ik]
    icol = ik // 2
    xS = x0 + icol*dx
    yS = y0 - dy*(ik%nlns)
    xE = xS + xl
    xT = xE + 0.2*xl
    expt_name = msisrlx.irlx_tests_name(expt_nmb, regn)
    ax7.plot([xS,xE],[yS, yS],
             linestyle=lstl, linewidth=lnw, marker=marker, markersize=mksz, color=clr1)
    ax7.text(xT,yS,expt_name, va='center')
  ik += 1
  icol = ik // 2
  xS = x0 + icol*dx
  xE = xS + xl
  xT = xE + 0.2*xl
  yS = y0 - dy*(ik%nlns)
  ax7.plot([xS,xE],[yS, yS],'-', linewidth=lnw0, color=pms_clr)
  ax7.text(xT,yS,f'PIOMAS',va='center')
  ax7.set_xlim([0.,0.8])
  ax7.set_ylim([0.15,0.5])
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
  subregn   = 'BerS'

# Plot time series and monthly stat for RMSE and bias, Bering Sea:
print("Plotting ...")
fgnmb=1
time_yrs = DV[:,0]+(DV[:,1]-1)/12
fig1,ax11,ax12,ax13,ax14,ax15, ax16, ax17 = plot_ice_stat(fgnmb, IAREA_Ber_nep, IAREA_Ber_pms, \
                                         IVOL_Ber_nep, IVOL_Ber_pms, \
                                         RMSE_Ber_ai, RMSE_Ber_hi, BIAS_Ber_ai, \
                                         BIAS_Ber_hi, 'BerS', lcntrl=True)

fig2,ax21,ax22,ax23,ax24,ax25, ax26, ax27 = plot_ice_stat(2, IAREA_Arc_nep, IAREA_Arc_pms, \
                                         IVOL_Arc_nep, IVOL_Arc_pms, \
                                         RMSE_Arc_ai, RMSE_Arc_hi, BIAS_Arc_ai, \
                                         BIAS_Arc_hi, 'ArcOc', lcntrl=True)

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
 


# Estimate S change due to ice vol change:
Ivol_expt1 = IVOL_Arc_nep[:,0]
Ivol_expt2 = IVOL_Arc_nep[:,1]   
Sice = 3.4
Soc = 30.
dltIvol = Ivol_expt1-Ivol_expt2 # km3
Vfw = dltIvol*(1.-Sice/Soc)
dltZ = 50*1e-3  # ocean layer, km
Voc = np.sum(Aarc*dltZ) # ocean vol, km3

dltS = -Soc*Vfw/Voc


