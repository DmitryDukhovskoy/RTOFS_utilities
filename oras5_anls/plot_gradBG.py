"""
  Plot ssh gradient that determines the BG intensity
  following Prosh & Johnson
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
from copy import copy
import matplotlib.colors as colors
from yaml import safe_load
import argparse
import pandas as pd
import time 

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
import mod_interp1D as mint1d
import mod_oras as moras
importlib.reload(moras)


parser = argparse.ArgumentParser()
parser.add_argument("--yrs", help="start year", type=int)
parser.add_argument("--yre", help="end year", type=int)
args = parser.parse_args()

YRS = 1988
YRE = 2024
if args.yrs:
  YRS = args.yrs
  YRE = YRS
if args.yre:
  YRE = args.yre


pthoras = '/work/Dmitry.Dukhovskoy/data/ORAS5'
pthout = '/work/Dmitry.Dukhovskoy/anls_output/oras5/BG_anls'
fout = f'BG_gradH_{YRS}-{YRE}.nc'
dfout = os.path.join(pthout,fout)
print(f'Loading {dfout}')

cff = 1.e7
cfar = 1.e-6
ds_bg = xarray.open_dataset(dfout)
TMd  = ds_bg['time'].values
HMAX = ds_bg['ssh_max'].values
GHMN = ds_bg['gradh_mean'].values*cff  # [m/m]
BGAR = ds_bg['BG_area'].values*cfar    # km2

# Delete unfinished years, if any:
ibad = np.where(TMd < 0)[0]
if len(ibad) > 0:
  ibad0 = ibad[0]
  ibad0 = (ibad0 // 12) * 12  # Keep only completed years
  yrs_finished = ibad0 // 12 - 1
  YRE = YRS + yrs_finished
  TMd = TMd[:ibad0]
  HMAX = HMAX[:ibad0]
  GHMN = GHMN[:ibad0]
  BGAR = BGAR[:ibad0]


dnmb_ref = mtime.datenum([1900,1,1])
TM = TMd+dnmb_ref
DV = mtime.datevec1D(TM)

years = DV[0]
months = DV[1]


# Annual means:
def compute_annual_mean(A1d):
  nrec = len(A1d)
  nyrs = nrec // 12
  A2d = np.reshape(A1d,(12,nyrs), order='F')
  Amn = np.nanmean(A2d, axis=0)

  return Amn

GHMN_ann = compute_annual_mean(GHMN)
BGAR_ann = compute_annual_mean(BGAR)
HMAX_ann = compute_annual_mean(HMAX)
yrs_plt = np.arange(YRS,YRE+2)
nyrs = YRE-YRS+1

clrmp = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBuGn_r')
rmin = -0.5
rmax = 0.5

time_days = TMd-TMd[0]+15
# Find years:
nrec = len(TM)
time_years = np.zeros((nrec))
for ii in range(nrec):
  d0 = TM[ii]
  dv0 = mtime.datevec(d0)
  yr0 = dv0[0]
  _,jd0 = mtime.dnmb2jday(d0)
  nd = mtime.month_days(2,yr0)
  if nd == 29:
    ndays_yr = 366
  else:
    ndays_yr = 365
  time_years[ii] = yr0 + (jd0-1)/ndays_yr

xtks = [x for x in range(years[0],years[-1]+1)]
xtk_labels = [str(int(t)) if i % 2 == 0 else '' for i, t in enumerate(xtks)]

def plot_ts_annmean(ax0, time_years, yrs_plt, TS, TS_ann, xtks, xtk_labels, sttl='time series', \
                    clrts=[0,0.4,0.8], clrgh=[0.8,0.5,0]):
  ax0.plot(time_years, TS, color=clrts)
  # Plot annual means:
  for ii, dmm in enumerate(yrs_plt[:-1]):
    y1 = yrs_plt[ii]
    y2 = yrs_plt[ii+1]
    val = TS_ann[ii]
    ax0.plot([y1,y2],[val,val],'-',linewidth=1.5,color=clrgh)

  ax0.grid('on')
  ax0.set_xticks(xtks)
  ax0.set_xticklabels(xtk_labels)
  ax0.set_xlim([YRS,YRE+1])
  ax0.set_title(sttl)

  return ax0  

print(f'Plotting BG for {YRS} - {YRE}')
plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.55, 0.85, 0.4])
sttl1 = f'gradH m/m * {1/cff:.1e}'
ax1 = plot_ts_annmean(ax1, time_years, yrs_plt, GHMN, GHMN_ann,  xtks, xtk_labels, sttl=sttl1, \
                      clrts=[0.,0.5,0.9], clrgh = [0.9,0.4,0])

sttl2 = f'BG Area km2 * {1/cfar:.1e}'
ax2 = plt.axes([0.1,0.08, 0.85, 0.4])
ax2 = plot_ts_annmean(ax2, time_years, yrs_plt, BGAR, BGAR_ann, xtks, xtk_labels, sttl=sttl2, \
                      clrts=[0.,0.8,0.3], clrgh = [0.6,0.0,0.8])

btx = 'plot_gradBG.py'
bottom_text(btx, pos=[0.02, 0.02]) 

fig2 = plt.figure(2,figsize=(9,8))
plt.clf()
ax3 = plt.axes([0.1, 0.55, 0.85, 0.4])
sttl3 = 'max SSH, m'
ax3 = plot_ts_annmean(ax3, time_years, yrs_plt, HMAX, HMAX_ann, xtks, xtk_labels, sttl=sttl3, \
                      clrts=[1.,0.5,0.], clrgh = [0.,0.5,0.9])

btx = 'plot_gradBG.py'
bottom_text(btx)



