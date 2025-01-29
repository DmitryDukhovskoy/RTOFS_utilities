"""
  Boxplots
  Plot seasonal mean area cold pool on the Bering Sea shelf
  bottom water < 2C for T classes
  for seasonal_forecasts, GOFS, GLORYS 
  and compare to WOA23

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
import scipy.stats as stats

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
tplot = 4 # 1= -2:-1, 2 = -1:0, 3= 0:+1, 4= +1:+2, 5= -2:+2 (all classes) 

# Initial date
# Look at ens run #1 - the only ens. that has 5-day av. output fields
expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
YAVRG = [x for x in range(2005,2015)]

expt_nmb = 2   # 2 - seas f/casts with dailyOB


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

def reading_data(dflout):
  print(f'Loading monthly coldpool area --> {dflout}')
  with open(dflout, 'rb') as fid:
    CPA, TCLASS = pickle.load(fid)

  return CPA, TCLASS

def average_seasons(CPA):
  """
  Average monthly data by seasons: JFM, MAM, ...
  """
  kdm,jdm,idm = CPA.shape
  CPAs = np.zeros((kdm,4,idm))  # years, seasons, T classes
  for iseas in range(4):
    i1 = iseas*3
    i2 = i1+3
    CPAs[:,iseas,:] = np.mean(CPA[:,i1:i2,:], axis=1)

  return CPAs

def get_tclass(CPA, tplot):
  kdm, jdm, idm = CPA.shape
  if tplot > idm:
    grp = np.sum(CPA, axis=2)
  else:
    grp = CPA[:,:,tplot-1].squeeze()

  return grp

cff = 1.e-5
pthanls = pthseas['MOM6_NEP'][expt]['pthanls'].format(expt_nmb=expt_nmb)
ics = -1
for dpl in DATA_PLOT:
  ics += 1
  match dpl:
    case 'seasonal_fcst':
      expt_name = f'NEPphys_frcst_dailyOB-expt{expt_nmb:02d}'
      dflout = os.path.join(pthanls,f'Bering_coldpoolarea_fcst_{YAVRG[0]}-{YAVRG[-1]}.pkl')
      CPA, TCLASS = reading_data(dflout)
      CPAs = average_seasons(CPA)
      grp1 = get_tclass(CPAs, tplot)*cff 
    case 'glorys':
      expt_name = 'GLORYS12v1' 
      dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_glorys_{YAVRG[0]}-{YAVRG[-1]}.pkl')
      CPA, _ = reading_data(dflout)
      CPAs = average_seasons(CPA)
      grp2 = get_tclass(CPAs, tplot)*cff
    case 'gofs31':
      expt_name = 'GOFS3.1-53.X'
      dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_gofs31_{YAVRG[0]}-{YAVRG[-1]}.pkl')
      CPA, _  = reading_data(dflout)
      CPAs = average_seasons(CPA)
      grp3 = get_tclass(CPAs, tplot)*cff

  CPAs = CPAs*cff # km2 x 1e-5
#  CPAs = CPA*cff

  # Get median, percentiles:
  lprc   = 25.
  if tplot == 5: 
    CPT = np.sum(CPAs, axis=2)
  else:
    CPT = CPAs[:,:,tplot-1].squeeze()  

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

# Use WOA as a banchmark for comparison
# means
expt_name = 'WOA23'
dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_WOA23_{YAVRG[0]}-{YAVRG[-1]}.pkl')
CPAs, _  = reading_data(dflout)
if tplot > CPAs.shape[1]:
  woa = np.sum(CPAs, axis=1)*cff
else:
  woa = CPAs[:, tplot-1].squeeze()*cff
print(f'nTC = {nTC}')

# 1-sample t-test for testing H0: mean = mean_woa
ttest_fcst = np.zeros((4))
ttest_glor = np.zeros((4))
ttest_gofs = np.zeros((4))

for iseas in range(4):
  a1 = grp1[:,iseas]
  a2 = grp2[:,iseas]
  a3 = grp3[:,iseas]
  _, ttest_fcst[iseas] = stats.ttest_1samp(a=a1, popmean=woa[iseas])
  _, ttest_glor[iseas] = stats.ttest_1samp(a=a2, popmean=woa[iseas])
  _, ttest_gofs[iseas] = stats.ttest_1samp(a=a3, popmean=woa[iseas])
 

CLRS = np.array([[0., 0.4, 0.9],
                [0.8, 0., 1],
                [1.0, 0.5, 0]])

#Xtime = [x for x in range(1,13)]
Xtime = np.arange(1,5)

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])

# Plot box plots
xlbl = 'Months'
ylbl = f'km2 x {cff}'
if tplot<5:
  tclass = f'{TCLASS[tplot-1]:.1f} : {TCLASS[tplot]:.1f} C'
else:
  tclass = 'Total'

sttl = f'Coldpool area, Tclass={tclass} {min(YAVRG)}-{max(YAVRG)}'

btx = 'coldpool_area_Tclass_boxplotsWOA23.py'

ax1 = mutil.plot_boxplot_v2(ax1, AA, sttl=sttl, CLRS=CLRS, XX=Xtime, \
                    xlbl=xlbl, ylbl=ylbl, lgnd_names=DATA_PLOT, btx=btx)

# Show WOA mean:
clrwoa = [1,0,0]
for iseas in range(4):
  x1 = iseas+1
  x2 = x1+1
  yy = woa[iseas]
  ax1.plot([x1,x2],[yy,yy],'-',linewidth=2, color=clrwoa)

# T-tests results:
alf = 0.05
dltY = 0.015
dx = 0.05
ax3 = plt.axes([0.1, 0.28, 0.8, 0.12])
clrok = [0,1,0]
icc = -1
for igr in ['fcst','glor','gofs']:
  icc += 1
  match igr:
    case 'fcst':
      ttest = ttest_fcst
    case 'glor':
      ttest = ttest_glor
    case 'gofs':
      ttest = ttest_gofs 

  y0 = icc+0.5
  for imo in range(len(Xtime)):
    x0 = Xtime[imo]
    mu_same = ttest[imo] > alf  # if true: accept null hypothesis, mean = WOA mean

    if mu_same:
      ax3.plot(x0+0.5, y0, marker='o', ms=10, markerfacecolor=clrok, mec='none')
    else:
      ax3.plot(x0+0.5, y0, marker='x', ms=10, color=[1,0,0], markeredgewidth=2)

#import matplotlib.transforms
Yticklab = ['fcst','glor','gofs']
Xticklab = ['JFM','AMJ','JAS','OND']
ax3.set_xticks(Xtime)
ax3.set_xlim([Xtime[0],Xtime[-1]+1])
ax3.set_xticklabels(Xticklab, ha='left')
ax3.grid('on')
ax3.set_xlabel('Seasons')
ax3.set_ylim([0,3])
ax3.set_title(f'ANOVA tests, mean1=mean2, conf={alf:.2f}')
ax3.set_yticks([1,2,3])
ax3.set_yticklabels(Yticklab, va='top')





