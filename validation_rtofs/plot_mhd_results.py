"""
  Plot MHD statistics, results from:
  gulfstream_stat.py
  and calc_saturation_mhd.py
"""
import os
import numpy as np
import sys
import importlib
import matplotlib
import matplotlib.pyplot as plt
import datetime
import pickle

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hausdorff')

import mod_read_hycom as mhycom
import mod_misc1 as mmisc
importlib.reload(mmisc)
import mod_hausdorff_distance as mmhd
import mod_utils_fig as mufig
import mod_utils as mutil
import mod_rtofs as mrtofs
import mod_gulfstream as mgulf
importlib.reload(mgulf)
import mod_time as mtime
importlib.reload(mtime)
from mod_utils_fig import bottom_text

pthsave = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/data_outp/validation_rtofs/'

expt      = 'paraD5'
fl_prd    = 'MHDproduct.pkl'
fmhd_prd  = pthsave + fl_prd
fl_prl    = 'MHD' + expt + '.pkl'
fmhd_prl  = pthsave + fl_prl
fl_satr   = 'MHD_saturation_value.pkl'
fmhd_satr = pthsave + fl_satr

print('Loading ' + fmhd_prd)
with open(fmhd_prd,'rb') as fid:
  [MHDprd,TMprd] = pickle.load(fid)

print('Loading ' + fmhd_prl)
with open(fmhd_prl,'rb') as fid:
  [MHDprl,TMprl] = pickle.load(fid)

print('Loading ' + fmhd_satr)
with open(fmhd_satr,'rb') as fid:
  [MHDsatr,TMsatr] = pickle.load(fid)

# # of days wrt to production time
dTprd  = TMprd-TMprd[0]
dTprl  = TMprl-TMprd[0]
dTsatr = TMsatr-TMprd[0] 
day1str = mtime.datestr(TMprd[0], show_hr=False)
rdateS  = mtime.datestr(TMprd[0], show_hr=False)
rdateE  = mtime.datestr(TMprd[-1], show_hr=False)

# Problem with NAVO GStr extending into the Gulf at days 71, 72:
MHDprd = np.where(MHDprd>180.,np.nan,MHDprd)
MHDprl = np.where(MHDprl>180.,np.nan,MHDprl)

# Saturation Value:
satVal  = 0.95*np.mean(MHDsatr) # Definition of saturation value as 95% 
#satVal  = np.median(MHDsatr)
#satVal1 = np.percentile(MHDsatr,25)
#satVal2 = np.percentile(MHDsatr,75)

md1  = np.nanpercentile(MHDprl,50.)
iqr2 = np.nanpercentile(MHDprl,90.)
iqr1 = np.nanpercentile(MHDprl,10.)

mdP1  = np.nanpercentile(MHDprd,50.)
iqrP2 = np.nanpercentile(MHDprd,90.)
iqrP1 = np.nanpercentile(MHDprd,10.)

mqmin = round(min([iqr1,iqrP1]),2)
mqmax = round(max([iqr2,iqrP2]),2)


clr1 = [0., 0.4, 0.9]
clr2 = [0., 0.8, 0.2]
clr3 = [1., 0., 0.]


plt.ion()

fig1 = plt.figure(1,figsize=(9,8), constrained_layout=False)
fig1.clf()

ax1 = fig1.add_axes([0.1, 0.5, 0.8, 0.4])
ln1, = ax1.plot(dTprl, MHDprl, color=clr1, linewidth=2, label="paraD5")
ln2, = ax1.plot(dTprd, MHDprd, color=clr2, linewidth=2, label="Porduct")
ln3, = ax1.plot([dTsatr[0], dTsatr[-1]], [satVal, satVal], color=clr3, linewidth=2, label="SaturVal")


xt1 = 0
xt2 = np.max(dTprd)

ax1.set_xlim([xt1,xt2])
ax1.set_xticks(np.arange(0,xt2,5))
ax1.grid(True)
ax1.set_xlabel('Days since ' + day1str)
ax1.set_title('{0} & Production vs NAVO, MHD (km), Gulf Stream, {1}-{2}'.\
                format(expt,rdateS,rdateE))

# Plot boxplot:
dxx  = 0.1
xx0  = 1.
ax2 = fig1.add_axes([0.1, 0.08, 0.22, 0.35])
ax2.plot(xx0,md1,'.',markersize=20, color=clr1)
ax2.plot([xx0,xx0], [iqr1,iqr2], '-', linewidth=2, color=clr1)
ax2.plot([xx0-dxx,xx0+dxx], [iqr1,iqr1], '-', linewidth=2, color=clr1)
ax2.plot([xx0-dxx,xx0+dxx], [iqr2,iqr2], '-', linewidth=2, color=clr1)

xx0  = 2.
ax2.plot(xx0,mdP1,'.',markersize=20, color=clr2)
ax2.plot([xx0,xx0], [iqrP1,iqrP2], '-', linewidth=2, color=clr2)
ax2.plot([xx0-dxx,xx0+dxx], [iqrP1,iqrP1], '-', linewidth=2, color=clr2)
ax2.plot([xx0-dxx,xx0+dxx], [iqrP2,iqrP2], '-', linewidth=2, color=clr2)

ax2.set_xlim([0.5, 2.5])
ax2.set_ylim([mqmin-5,mqmax+5])
ax2.set_xticks([1,2])
ax2.grid(True)
ax2.set_xticklabels(['paraD5','Product'])

ax2.set_title('Inter-decile range, MHD(km)')


ax3 = plt.axes([0.7, 0.2, 0.25, 0.25])
lgd = plt.legend(handles=[ln1, ln2, ln3], loc='upper left')
ax3.axis('off')


btx = 'plot_mhd_results.py'
bottom_text(btx, pos=[0.02,0.01])







