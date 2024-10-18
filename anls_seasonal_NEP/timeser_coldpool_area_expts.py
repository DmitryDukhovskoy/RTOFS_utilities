"""
  Plot time series of the cold pool area for all expts

  see: coldpool_area.py
  in the Bering Sea
  bottom water < 2C
  deep water (> 250 m) is excluded
  following:
  On the variability of the Bering Sea Cold Pool and implications 
       for the biophysical environment
  2022
 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8979450/
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

nens = 2
YRS  = 1993
MMS  = 4
DDS  = 1
RUNS = [f'NEPphys_frcst_climOB_{YRS}-{MMS:02d}-e{nens:02d}',
        'NEP_physics_GOFS-IC',
        'NEPphys_LZRESCALE_climOB_1993_04-e02']
nruns = len(RUNS)
expt  = "seasonal_fcst"

CLRS = np.array([[1., 0.5, 0],
                [0., 0.4, 0.9],
                [0., 0.9, 0.2]])

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)


# Plot time series:
plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])

LNS  = []
Tmin = 1.e10
Tmax = -1.
pthdump = pthseas["MOM6_NEP"][expt]["pthdump"] 
for irun in range(nruns):
  runname = RUNS[irun]
  fldump  = f'{runname}_coldpool_area.pkl'
  drfldump = os.path.join(pthdump, fldump)
  print(f'Loading <--- {drfldump}')

  with open(drfldump, 'rb') as fid:
    TM, CPA = pickle.load(fid)


  cff   = 1.e-5
  CPsc  = CPA*cff
  Tmin  = min([Tmin,TM[0]])
  Tmax  = max([Tmax,TM[-1]]) 
  clr = CLRS[irun,:]
  ln1, = ax1.plot(TM, CPsc, color=clr, linewidth=2, label=runname)

  LNS.append(ln1)


# Show Month day1:
dv1   = mtime.datevec(Tmin)
dnmb1 = mtime.datenum([dv1[0],dv1[1],1])
dv2   = mtime.datevec(Tmax + 31)
dnmb2 = mtime.datenum([dv2[0],dv2[1],1])
TMp   = np.arange(dnmb1,dnmb2+1)
DVp   = mtime.datevec2D(TMp)
Ip = np.where(DVp[:,2] == 1)
Tticks  = TMp[Ip]
DVticks = DVp[Ip,:3].squeeze()
ntcks   = len(Tticks)
tck_lbls = []
for ik in range(ntcks):
  ff = f'{DVticks[ik,1]}/{DVticks[ik,2]}'
  tck_lbls.append(ff)
  

ax1.set_xticks(Tticks)
ax1.grid(True)
ax1.set_xticklabels(tck_lbls)
ax1.set_xlim([Tmin, Tmax])
ax1.set_ylabel(f'km2 x {cff}')
ax1.set_xlabel(f'Months starting {DVp[0,0]}/{DVp[0,1]}/{DVp[0,2]}')
ax1.set_title(f'NEP experiments, cold pool area km2 x {cff}')

ax3 = plt.axes([0.1, 0.2, 0.6, 0.22])
lgd = plt.legend(handles=LNS, loc='upper right')
ax3.axis('off')

btx = 'timeser_coldpool_area_expts.py'
bottom_text(btx, pos=[0.2, 0.3])


