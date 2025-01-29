"""
  Plot monhtly mean area cold pool on the Bering Sea shelf
  bottom water < 2C for T classes
  derived in:
   seasonal f/casts:  derive_coldpool_area_Tclass_seasfcst.py
   glorys: derive_coldpool_area_Tclass_glorys.py
   gofs3.1: derive_coldpool_area_Tclass_gofs31.py

  Plot bottom T and cold pool in the Bering Sea by seasons:
  bottom water < 2C
  e.g. On the variability of the Bering Sea Cold Pool and implications 
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
from matplotlib.patches import Polygon

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

#data_plot = 'seasonal_fcst'  # from seasonal f/casts output
#data_plot = 'glorys'   # from GLORYS data
data_plot = 'gofs31'   # from GOFS3.1 

# Initial date
# Look at ens run #1 - the only ens. that has 5-day av. output fields
expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
varnm    = 'salin'  # temp (potential) / salin
#dnmbS    = mtime.datenum([2015,1,1])
# Averaging time period:
MMS   = 1    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YAVRG = [x for x in range(2005,2015)]

nensR    = 1
expt_nmb = 2   # 2 - seas f/casts with dailyOB


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthanls = pthseas['MOM6_NEP'][expt]['pthanls'].format(expt_nmb=expt_nmb)
match data_plot:
  case 'seasonal_fcst':
    expt_name = f'NEPphys_frcst_dailyOB-expt{expt_nmb:02d}'
    dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_fcst_{YAVRG[0]}-{YAVRG[-1]}.pkl')
  case 'glorys':
    expt_name = 'GLORYS12v1' 
    dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_glorys_{YAVRG[0]}-{YAVRG[-1]}.pkl')
  case 'gofs31':
    expt_name = 'GOFS3.1-53.X'
    dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_gofs31_{YAVRG[0]}-{YAVRG[-1]}.pkl')
 
run_info = f'{expt_name} init MM={MMS} e{nensR:02d}, conservT bottom: {min(YAVRG)}-{max(YAVRG)}' 

print(f'Plotting {varnm} {expt_name} ')
print(f'{run_info}')


print(f'Loading monthly coldpool area --> {dflout}')
with open(dflout, 'rb') as fid:
  CPA, TCLASS = pickle.load(fid)
nTC    = len(TCLASS)

cff = 1.e-5
CPA = CPA*cff # km2 x 1e-5

# Get median, percentiles:
lprc   = 10.
CPmed  = np.median(CPA, axis=0)
CPlprc = np.percentile(CPA, lprc, axis=0)
CPuprc = np.percentile(CPA, (100-lprc), axis=0)

CPT   = np.sum(CPA, axis=2)
CPTm  = np.median(CPT, axis=0)
CPTl  = np.percentile(CPT, lprc, axis=0)
CPTu  = np.percentile(CPT, (100-lprc), axis=0)

plt.ion()

CLRS = np.array([[0., 0.2, 0.9],
                 [0., 0.8, 1],
                 [0.7, 0., 1],
                 [0.9, 0.4, 0],
                 [0.5, 0.3, 0]])
CPRC = np.array([[0.9,  0.93, 1],
                 [0.9,  0.99, 1],
                 [0.97, 0.9, 1],
                 [0.95, 0.93, 0.9]])

#Xtime = [x for x in range(1,13)]
Xtime = np.arange(1,13)

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
LNS = []
for icl in range(nTC-1):
  clr0 = CLRS[icl, :]
  CPa  = CPmed[:, icl]
  t1 = TCLASS[icl]
  t2 = TCLASS[icl+1]
  tline = f'{t1:.1f} < t < {t2:.1f}'
  ln1, = ax1.plot(Xtime, CPa, linewidth=2, color=clr0, label=tline)
  LNS.append(ln1)
#
# Total area for all classes:
clr0 = [0,0,0]
ln1, = ax1.plot(Xtime, CPTm, linewidth=2, color=clr0, label='Total')
LNS.append(ln1)

# Show percentiles:
for icl in range(nTC-1):
  CPlprc_cold = CPlprc[:,icl]
  CPuprc_cold = CPuprc[:,icl]
  verts = [*zip(Xtime,CPuprc_cold), *zip(np.flip(Xtime),np.flip(CPlprc_cold))]
  fclr = CPRC[icl,:]
  poly = Polygon(verts, facecolor=fclr)
  ax1.add_patch(poly)

# Percentiles for total CP area:
Tverts = [*zip(Xtime,CPTu), *zip(np.flip(Xtime),np.flip(CPTl))]
Tpoly = Polygon(Tverts, facecolor=(0.95,0.95,0.95))
ax1.add_patch(Tpoly)

ax1.set_xticks(Xtime)
ax1.set_yticks([x for x in range(0,10)])
ax1.grid(True)
ax1.set_xlim([0.9, 12.1])
ax1.set_ylim([-0.2, 9.5])
#ax1.set_xticklabels(tck_lbls)
ax1.set_ylabel(f'km2 x {cff}')
ax1.set_xlabel(f'Months')
ax1.set_title(run_info)

ax3 = plt.axes([0.1, 0.2, 0.6, 0.22])
lgd = plt.legend(handles=LNS, loc='upper right')
ax3.axis('off')

btx = 'coldpool_area_Tclass_month.py'
#btx = 'coldpool_area_Tclass_seasfcst.py'
bottom_text(btx, pos=[0.1, 0.1])



