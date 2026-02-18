"""
  Calculate trends in winter (May) and summer (Sept) SSS and SST
  inside the relaxation zone - Arctic NEP10k only
  to check for the drifts due to ice relaxation

  calc over the time that = the shortest time series

  from PHYS and BGC + IRLX h/casts - both with GLORYS nudging

  from PHYS no GLORYS but there is IRLX

  NOAA NWS EMC Dmitry Dukhovskoy 
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib
import xarray
import matplotlib.colors as colors
from yaml import safe_load
import argparse
from scipy.stats import linregress

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
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hausdorff')

import mod_time as mtime

parser = argparse.ArgumentParser()
parser.add_argument("--yrs", help="year to start calc trend, default=1993", type=int)
parser.add_argument("--yre", help="year to end, default=2019", type=int)
args = parser.parse_args()

YRS = args.yrs if args.yrs else 1993
YRE = args.yre if args.yre else 2019


EXPTS = ['NEPbgc_nudged_hindcast02',
         'NEPphys_nudged_hindcast',
         'NEPphys_nonudg_irlx_hcast']
nexp = len(EXPTS)

EXPTS_INFO=['NEPbgc GLORYS IRLX12hrs',
            'NEPphys GLORYS no IRLX',
            'NEPphys no GLORYS, IRLX24hrs']


def get_timeyr(TM):
  TM = np.array(TM)
  DV = mtime.datevec1D(TM, fHR=False)
  DV = np.array(DV).transpose()
  time_yrs = DV[:,0]+(DV[:,1]-1)/12

  return time_yrs

regn_name = 'NEP10k Arctic Chukchi'
#expt_nameB = 'NEPbgc_nudged_hindcast02'
#expt_nameP = 'NEPphys_nudged_hindcast'

xtck = [x for x in range(YRS,YRE+1)]

def read_output(flout):
  pthoutp = '/work/Dmitry.Dukhovskoy/anls_output/NEPbgc_hindcast02'
  dflout = os.path.join(pthoutp,flout)
  print(f'Reading  {dflout}')
  DTM = np.load(dflout)
  SFLD = DTM['SFLD']
  TFLD = DTM['TFLD']
  TM   = DTM['TM']

  return TM, SFLD, TFLD

def derive_timeser(F1d, YRS, YRE, MM, TM):
  DV = mtime.datevec1D(TM)
  YEARS = DV[0]
  MONTHS = DV[1]
  # Check that requested years are there:
  assert YRS >= YEARS[0], f'Requested year {YRS} is outside saved time {YEARS[0]}'
  assert YRE <= YEARS[-1], f'Requested year {YRE} is outside saved time {YEARS[-1]}'
  iS = np.where(YEARS == YRS)[0][0]
  iE = np.where(YEARS == YRE)[0][-1]

  # Subset data for year range and filter by months:
  TSsub = F1d[iS:iE+1]
  MMsub = MONTHS[iS:iE+1]
  TS = TSsub[MMsub == MM]

  return TS

def lregr_stat(TSER):
  nrec = len(TSER)
  XX = [x for x in range(1,nrec+1)]

  # Compute trend:
  LSF = linregress(XX,TSER)
  alf0 = LSF.intercept
  alf1 = LSF.slope
  p_val = LSF.pvalue

  # Results
  print(f"Slope: {LSF.slope}")
  print(f"Intercept: {LSF.intercept}")
  print(f"R-squared: {LSF.rvalue**2}")
  print(f"P-value: {LSF.pvalue}")
  #print(f"Standard error (slope): {LSF.stderr}")

  STAT = np.zeros((3))
  STAT[0] = alf0
  STAT[1] = alf1
  STAT[2] = p_val

  return STAT
  

CLR = [[0.,0.3,1],
       [0.9,0.3,0],
       [0.,0.9,0.2],
       [1.,0.9,0],
       [0.8,0.,0.5],
       [0.7, 1, 0.2]]


nyrs = YRE-YRS+1
MW = 5
MS = 9
STAT_SW = np.zeros((3,nexp))  # SSS winter stat: alf0, alf1, p-value
STAT_SS = np.zeros((3,nexp))  # SSS summer
STAT_TW = np.zeros((3,nexp))
STAT_TS = np.zeros((3,nexp))
SSTW = np.zeros((nyrs,nexp))
SSTS = np.zeros((nyrs,nexp))
SSSW = np.zeros((nyrs,nexp))
SSSS = np.zeros((nyrs,nexp))
for iexp in range(nexp):
  expt_name = EXPTS[iexp]
  floutP = f'{expt_name}_surfTS_NEParct.npz'
  TM, SFLD, TFLD = read_output(floutP)
  #tyrs = get_timeyr(TM)

  SSS_wint = derive_timeser(SFLD, YRS, YRE, MW, TM)
  SSS_summ = derive_timeser(SFLD, YRS, YRE, MS, TM)
  SST_wint = derive_timeser(TFLD, YRS, YRE, MW, TM)
  SST_summ = derive_timeser(TFLD, YRS, YRE, MS, TM)

  STAT = lregr_stat(SSS_wint)
  STAT_SW[:,iexp] = STAT
  STAT = lregr_stat(SSS_summ)
  STAT_SS[:,iexp] = STAT
  STAT = lregr_stat(SST_wint)
  STAT_TW[:,iexp] = STAT
  STAT = lregr_stat(SST_summ)
  STAT_TS[:,iexp] = STAT

  SSTW[:,iexp] = SST_wint
  SSTS[:,iexp] = SST_summ
  SSSW[:,iexp] = SSS_wint
  SSSS[:,iexp] = SSS_summ

def plot_trends(ax1,sttl,STAT, FLDS, sinfo=''):
  """
    Plot least sq. fit to the SST/SSS data
  """
  tyrs = [x for x in range(YRS,YRE+1)]
  hndls = []
  for iexp in range(nexp):
    expt_info = EXPTS_INFO[iexp]
    expt_name = EXPTS[iexp]
    clr1 = CLR[iexp]
    ax1.plot(tyrs, FLDS[:,iexp],'.', markersize=12, color=clr1)
    # Regr line:
    alf0 = STAT[0,iexp]
    alf1 = STAT[1,iexp]
    pval = STAT[2,iexp]
    lfit = alf0 + alf1*XX

    ln1, = ax1.plot(tyrs, lfit, '-', linewidth=2, color=clr1, label=expt_info)
    hndls.append(ln1)

    sinfo += f'{expt_name}:\n'
    str   = f'alf0={alf0:.4f}\n' +\
            f'alf1={alf1:.4f}\n' +\
            f'pval={pval:.6f}\n'
    sinfo += str

    ax1.set_xticks(tyrs)
    ax1.set_xlim([YRS-0.1,YRE+0.1])
    ax1.grid('on')
    ax1.set_title(sttl)

  return hndls, sinfo


expt_info = EXPTS_INFO[iexp]
tyrs = [x for x in range(YRS,YRE+1)]
nrec = len(tyrs)
XX = np.arange(1,nrec+1)

plt.ion()

# Plot SSS winter/ summer
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.7, 0.8, 0.25])
sttl = f'hcasts SSS wint' 
hndls, sinfo_sw = plot_trends(ax1,sttl,STAT_SW,SSSW, sinfo='SSS wint\n') 
  
ax2 = plt.axes([0.1, 0.35, 0.8, 0.25])
sttl = f'hcasts SSS summ' 
_, sinfo_ss = plot_trends(ax2,sttl,STAT_SS,SSSS, sinfo='SSS summ\n') 


# Legend
ax2 = plt.axes([0.63, 0.05, 0.35, 0.08])
ax2.legend(handles=hndls, loc='lower right')
ax2.axis('off')

# Info:
ax3 = plt.axes([0.05, 0.03, 0.45,0.25])
ax3.text(0,0,sinfo_sw)
ax3.text(0.6,0,sinfo_ss)
ax3.axis('off')

btx = 'calc_trends_surfTS.py'
bottom_text(btx, pos=[0.02,0.02])

# Plot SST winter/ summer
fig2 = plt.figure(2,figsize=(9,8))
plt.clf()
ax21 = plt.axes([0.1, 0.7, 0.8, 0.25])
sttl = f'hcasts SST wint' 
hndls, sinfo_tw = plot_trends(ax21,sttl,STAT_TW,SSTW, sinfo='SST wint\n') 
  
ax22 = plt.axes([0.1, 0.35, 0.8, 0.25])
sttl = f'hcasts SST summ' 
_, sinfo_ts = plot_trends(ax22,sttl,STAT_TS,SSTS, sinfo='SST summ\n') 


# Legend
ax23 = plt.axes([0.63, 0.05, 0.35, 0.08])
ax23.legend(handles=hndls, loc='lower right')
ax23.axis('off')

# Info:
ax24 = plt.axes([0.05, 0.03, 0.45,0.25])
ax24.text(0,0,sinfo_tw)
ax24.text(0.6,0,sinfo_ts)
ax24.axis('off')


