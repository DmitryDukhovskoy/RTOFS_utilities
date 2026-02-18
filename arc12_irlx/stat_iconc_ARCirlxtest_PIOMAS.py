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
parser.add_argument(
    "--enmb",
    help="List of experiment numbers (e.g., 1 3 5 32 33)",
    type=int,
    nargs="+",             # <-- allows one or more integers and will generate a list
    required=True
)
args = parser.parse_args()

# Test runs were performed for only 1 year
YRS = args.yr if args.yr else 1995
interp = args.intrp if args.intrp else 0  
ENMBS = args.enmb if args.enmb else None

# Number of test runs:
Nexpts = len(ENMBS)

interp_mnthly = interp > 0  # for more accurate comparison, do time interpolation of PIOMAS 
                      # to get mnthly mean values, similar to how it is done in SIS2
                      # when deriving iconc ithkn for day=d0 from PIOMAS target fields
#use_mnth = not (args.fday and args.fday > 0)
regn = 'ARC'
lat_rlx = 60.  # compute statistics inside relax zone north of lat_rlx, make it 0 to ignore
use_mnth = True
plt_rgn = False # Show Arc and Ber regions

if lat_rlx > 1.e-16:
  print(f' ===  NOTE: statistics computed inside rlx zone north of {lat_rlx:.1f}N  ===')


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

#EXPTS = [1,2,3,4,5,21]   # 21 - 1hr rlx, no ridging
#EXPTS = [1,2,3,4,5]
#Nexpts = len(EXPTS)
RMSE_Arc_ai = np.zeros((12,Nexpts))
RMSE_Arc_hi = np.zeros((12,Nexpts))
BIAS_Arc_ai = np.zeros((12,Nexpts))
BIAS_Arc_hi = np.zeros((12,Nexpts))
IAREA_arc = np.zeros((12,Nexpts))
IVOL_arc  = np.zeros((12,Nexpts))
IAREA_pms = np.zeros((12))
IVOL_pms  = np.zeros((12))
TM = []

for kexpt in range(Nexpts):
  expt_nmb = ENMBS[kexpt]

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
      CIpms,_ = msisrlx.read_relax_piomas(dnmb0, pthrlx, 'iconc', dfpiomas=dfpiomas) 
      HIpms,_ = msisrlx.read_relax_piomas(dnmb0, pthrlx, 'ithkn', dfpiomas=dfpiomas) 

    if lat_rlx > 1.e-16:
      CIarc[hlat<lat_rlx] = np.nan
      HIarc[hlat<lat_rlx] = np.nan
      CIpms[hlat<lat_rlx] = np.nan
      HIpms[hlat<lat_rlx] = np.nan

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
    VolIce_arc  = HIarc*CIarc*Acell*1e-6   # x10^3 km3, HI - m, Acell - km
    VolA_arc = np.nansum(VolIce_arc[JA,IA])

    VolIce_pms  = HIpms*CIpms*Acell*1e-6   # x10^3 km3, HI - m, Acell - km
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

IAREA_arc = IAREA_arc*1e-6  # 1e6 km2
IAREA_pms = IAREA_pms*1e-6  # 1e6 km2

TM = np.array(TM)
DV = mtime.datevec1D(TM, fHR=False)
DV = np.array(DV).transpose()

#ECOLR = [[0.,0.4,0.9],
#         [0.9,0.5,0],
#         [0.,0.9,0.7],
#         [1.,0.9,0],
#         [0.8,0.,0.5],
#         [0.7, 1, 0.2]]
#
ECOLR = msisrlx.irlx_tests_colors()

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
    lcntrl - show control run with a special line
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

  btx = 'stat_iconc_ARCirlxtest_PIOMAS.py'
  bottom_text(btx, pos=[0.02,0.01])

  return fig1, ax1, ax2, ax3, ax4, ax5, ax6, ax7

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

if lat_rlx > 1.e-16:
  regn_name = f'ARC rlx{lat_rlx:.1f}N'
else:
  regn_name = 'ARC10k norlx'

fig1,ax11,ax12,ax13,ax14,ax15, ax16, ax17 = plot_ice_stat(1, IAREA_arc, IAREA_pms, \
                                    IVOL_arc, IVOL_pms, \
                                    RMSE_Arc_ai, RMSE_Arc_hi, BIAS_Arc_ai, 
                                    BIAS_Arc_hi, regn_name, lcntrl=True)

  

