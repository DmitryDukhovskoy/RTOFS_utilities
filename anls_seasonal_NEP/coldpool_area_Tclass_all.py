"""
  Monthly cold pool areas from MOM6-SIS2, GLORYS, GOFS3.1
  Plot monhtly mean area cold pool on the Bering Sea shelf
  bottom water < 2C for T classes
  for seasonal_forecasts, GOFS, GLORYS 
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
from scipy.stats import f_oneway

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


DATA_PLOT = ['seasonal_fcst', 'glorys', 'gofs31']
tplot = 3 # 1= -2:-1, 2 = -1:0, 3= 0:+1, 4= +1:+2, 5= -2:+2 (all classes) 

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

def reading_data(dflout):
  print(f'Loading monthly coldpool area --> {dflout}')
  with open(dflout, 'rb') as fid:
    CPA, TCLASS = pickle.load(fid)

  return CPA, TCLASS

def get_tclass(CPA, tplot):
  kdm, jdm, idm = CPA.shape
  if tplot > idm:
    grp = np.sum(CPA, axis=2)
  else:
    grp = CPA[:,:,tplot-1].squeeze()

  return grp

pthanls = pthseas['MOM6_NEP'][expt]['pthanls'].format(expt_nmb=expt_nmb)
ics = -1
for dpl in DATA_PLOT:
  ics += 1
  match dpl:
    case 'seasonal_fcst':
      expt_name = f'NEPphys_frcst_dailyOB-expt{expt_nmb:02d}'
      dflout = os.path.join(pthanls,f'Bering_coldpoolarea_fcst_{YAVRG[0]}-{YAVRG[-1]}.pkl')
      CPA, TCLASS = reading_data(dflout)
      grp1 = get_tclass(CPA, tplot) 
    case 'glorys':
      expt_name = 'GLORYS12v1' 
      dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_glorys_{YAVRG[0]}-{YAVRG[-1]}.pkl')
      CPA, _ = reading_data(dflout)
      grp2 = get_tclass(CPA, tplot) 
    case 'gofs31':
      expt_name = 'GOFS3.1-53.X'
      dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_gofs31_{YAVRG[0]}-{YAVRG[-1]}.pkl')
      CPA, _  = reading_data(dflout)
      grp3 = get_tclass(CPA, tplot) 

  cff = 1.e-5
  CPA = CPA*cff # km2 x 1e-5

  # Get median, percentiles:
  lprc   = 25.
  if tplot == 5: 
    CPT = np.sum(CPA, axis=2)
  else:
    CPT = CPA[:,:,tplot-1].squeeze()  

  CPmed  = np.median(CPT, axis=0)
  CPmean = np.mean(CPT, axis=0)
  CPlprc = np.percentile(CPT, lprc, axis=0)
  CPuprc = np.percentile(CPT, (100-lprc), axis=0)
  CPmin  = np.min(CPT, axis=0)
  CPmax  = np.max(CPT, axis=0)

  if ics == 0:
    AA  = np.zeros((len(DATA_PLOT), 6, len(CPmed)))

#    STAT = np.array([CPmed, CPlprc, CPuprc, CPmin, CPmax])
  AA[ics,:,:] = np.array([CPmed, CPlprc, CPuprc, CPmin, CPmax, CPmean])

nTC    = len(TCLASS)

# Perform 1-way ANOVA to test if the means btw seas. f/casts and GLORYS/GOFS are the same
# Note however that normality assumption is not held for all months/groups
# ANOVA hypothesis: 
# H0 (null hypothesis): μ1 = μ2 = μ3 = … = μk (It implies that the means of all the population are equal)
# H1 (alternative): It states that there will be at least one population mean that differs from the rest
anova12_pval = np.zeros((12))  # anova between seasonal f/cast and glorys
anova13_pval = np.zeros((12))  # -"-  -"- f/cast and gofs31
anova23_pval = np.zeros((12))  # -"-  -"- glorys/gofs

for imo in range(12):
  a1 = grp1[:,imo]
  a2 = grp2[:,imo]
  a3 = grp3[:,imo]
  anova12_pval[imo] = f_oneway(a1,a2)[1]
  anova13_pval[imo] = f_oneway(a1,a3)[1]
  anova23_pval[imo] = f_oneway(a2,a3)[1]


CLRS = np.array([[0., 0.4, 0.9],
                [0.8, 0., 1],
                [1.0, 0.5, 0]])

CPRC = np.array([[0.9,  0.94, 1],
                 [0.98,  0.9, 1],
                 [1, 0.95, 0.9]])

Xtime_tck = [x for x in range(1,13)]
Xtime = np.arange(1,13) + 0.5

xlbl = 'Months'
ylbl = f'km2 x {cff}' 
if tplot<5:
  tclass = f'{TCLASS[tplot-1]:.1f} : {TCLASS[tplot]:.1f} C'
else:
  tclass = 'Total'

sttl = f'Coldpool area, Tclass={tclass} {min(YAVRG)}-{max(YAVRG)}' 


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])


# Show percentiles for MOM6 only:
#for igrp in range(len(DATA_PLOT)):
for igrp in range(0,1):
#  CPlprc = AA[igrp,3,:].squeeze()  # min
#  CPuprc = AA[igrp,4,:].squeeze()  # max
  CPlprc = AA[igrp,1,:].squeeze()  # 25 prc
  CPuprc = AA[igrp,2,:].squeeze()  # 75 prc
  verts = [*zip(Xtime,CPuprc), *zip(np.flip(Xtime),np.flip(CPlprc))]
  fclr = CPRC[igrp,:]
  poly = Polygon(verts, facecolor=fclr)
  ax1.add_patch(poly)

LNS = []
for igrp in range(len(DATA_PLOT)):
  clr0   = CLRS[igrp, :]
  CPmn   = AA[igrp,5,:].squeeze()
  tline  = DATA_PLOT[igrp]
  ln1, = ax1.plot(Xtime, CPmn, linewidth=2, color=clr0, label=tline)
  LNS.append(ln1)
  ax1.plot(Xtime, CPmn, 'o', markerfacecolor=clr0, mec='none')
#
ax1.set_xticks(Xtime_tck)
#ax1.set_yticks([x for x in range(0,10)])
ax1.grid(True)
ax1.set_xlim([1, 13])
#ax1.set_ylim([-0.2, 9.5])
#ax1.set_xticklabels(tck_lbls)
ax1.set_ylabel(f'km2 x {cff}')
ax1.set_xlabel(f'Months')
ax1.set_title(sttl)

# Plot ANOVA tests for simulation groups:
alf = 0.05
dltY = 0.015
dx = 0.05
ax3 = plt.axes([0.1, 0.28, 0.8, 0.12])
clrok = [0,1,0]
icc = -1
for igr in ['12','13','23']:
  icc += 1
  match igr:
    case '12':
      fanova = anova12_pval
      clr1 = CLRS[0]
      clr2 = CLRS[1]
    case '13':
      fanova = anova13_pval
      clr1 = CLRS[0]
      clr2 = CLRS[2]
    case '23':
      fanova = anova23_pval
      clr1 = CLRS[1]
      clr2 = CLRS[2]

  y0 = icc+0.5
  for imo in range(12):
    x0 = Xtime_tck[imo]
    mu_same = fanova[imo] > alf  # if true: accept null hypothesis, means are the same

#    ax3.plot(x0, y0, marker='s', ms=10, markerfacecolor=clr1, mec='none')
#    ax3.plot(x0+dx, y0, marker='s', ms=10, markerfacecolor=clr2, mec='none')
    if mu_same:
#      ax3.plot(x0+2.5*dx, y0, marker='o', ms=10, markerfacecolor=clrok, mec='none')
      ax3.plot(x0+0.5, y0, marker='o', ms=10, markerfacecolor=clrok, mec='none')
    else:
      ax3.plot(x0+0.5, y0, marker='x', ms=10, color=[1,0,0], markeredgewidth=2)

#import matplotlib.transforms
Yticklab=['fcst-glor','fcst-gofs','glor-gofs']
ax3.set_xticks(Xtime_tck)
ax3.set_xlim([Xtime_tck[0],Xtime_tck[-1]+1])
ax3.grid('on')
ax3.set_xlabel('Months')
ax3.set_ylim([0,3])
ax3.set_title(f'ANOVA tests, mean1=mean2, conf={alf:.2f}')
ax3.set_yticks([1,2,3])
ax3.set_yticklabels(Yticklab, va='top')
  
ax2 = plt.axes([0.5, 0.1, 0.4, 0.2])
lgd = plt.legend(handles=LNS, loc='lower right')
ax2.axis('off')

btx = 'coldpool_area_Tclass_all.py'
bottom_text(btx, pos=[0.1, 0.1])





