"""
  Plot all models
  T/S Profile from model simulations 
  subsampled at locations of WOD observations

  First, run extract_TSprofWOD_{mom6, gofs, rtofs}.py

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from netCDF4 import Dataset as ncFile

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
import mod_mom6_valid as mom6vld
import mod_misc1 as mmisc
import mod_time as mtime
import mod_WODdata as mwod
importlib.reload(mwod)
importlib.reload(mom6vld)

from matplotlib.patches import Polygon

# WOD data
pthwod = '/scratch1/NCEPDEV/stmp4/Dmitry.Dukhovskoy/WOD_profiles/'
pthoutp= '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/MOM6_CICE6/ts_prof/'

modrun = 'rtofs'
YR1    = 2018
YR2    = 2022
mo1    = 1
mo2    = 12
#regn   = 'AmundsAO'
#regn   = 'NansenAO'
#regn   = 'MakarovAO'
regn  = 'CanadaAO'

REGNS = mom6vld.ts_prof_regions()
x0 = REGNS[regn]["lon0"]
y0 = REGNS[regn]["lat0"] 
dx = REGNS[regn]["dlon"]
dy = REGNS[regn]["dlat"]

print(f'Plotting WOD T/S for selected UID {regn} {YR1} {YR2} mo={mo1}-{mo2}\n')

flts_wod   = f'WODTS_{regn}_{YR1}-{YR2}.pkl'
flts_rtofs = f'rtofsTS_{regn}_{YR1}-{YR2}.pkl'
flts_gofs  = f'gofsTS_{regn}_{YR1}-{YR2}.pkl'
flts_mom   = f'mom6TS_{regn}_{YR1}-{YR2}.pkl'
dflts_wod  = os.path.join(pthoutp, flts_wod)
dflts_rtofs = os.path.join(pthoutp, flts_rtofs)
dflts_gofs  = os.path.join(pthoutp, flts_gofs)
dflts_mom = os.path.join(pthoutp, flts_mom)

# Load T/S profiles
prct = 5
plow = prct
pup  = 100-plow
ZZi, Tmn_wod, Tpl_wod, Tpu_wod, Smn_wod, Spl_wod, Spu_wod = \
           mom6vld.read_prof_stat(dflts_wod, plow=prct)
_, Tmn_mom, _, _, Smn_mom, _, _     = mom6vld.read_prof_stat(dflts_mom, plow=prct)
_, Tmn_rtofs, _, _, Smn_rtofs, _, _ = mom6vld.read_prof_stat(dflts_rtofs, plow=prct)
_, Tmn_gofs, _, _, Smn_gofs, _, _   = mom6vld.read_prof_stat(dflts_gofs, plow=prct)

izm = min(np.where((ZZi < -100.) & (np.isnan(Tmn_wod)))[0])
#zlim = np.floor(ZZi[izm])
zlim = -900.

# Patch the 5-95 percentile region:
tlow = np.flipud(Tpl_wod[np.where(~np.isnan(Tpl_wod))])
zlow = np.flipud(ZZi[np.where(~np.isnan(Tpl_wod))])
tup  = Tpu_wod[np.where(~np.isnan(Tpu_wod))]
zup  = ZZi[np.where(~np.isnan(Tpu_wod))]

slow = np.flipud(Spl_wod[np.where(~np.isnan(Spl_wod))])
sup  = Spu_wod[np.where(~np.isnan(Spu_wod))] 
tverts = [*zip(tlow,zlow), *zip(tup,zup)]
sverts = [*zip(slow,zlow), *zip(sup,zup)]

btx  = 'plot_TSprof_allModvsWOD.py'

plt.ion()
clrW_mn  = [0., 0.0, 0.]
clrM_mn  = [0.0, 0.4, 0.9]
clrR_mn  = [0.2, 0.9, 0.]
clrG_mn  = [1., 0.7, 0.]
clrW_pct = [.6, 0.6, 0.6]


fig1 = plt.figure(1,figsize=(9,9))
plt.clf()

sttl = f'Models vs WOD T {YR1}-{YR2} {regn}'

ax1   = plt.axes([0.08, 0.15, 0.4, 0.8])
ln1,  = ax1.plot(Tmn_wod, ZZi, '-', linewidth=2, color=clrW_mn, label="Obs")
ln2,  = ax1.plot(Tpl_wod,ZZi,'-', color=clrW_pct, label=f"{plow}-{pup}%")
tpoly = Polygon(tverts, facecolor='0.6', edgecolor='none')
ax1.add_patch(tpoly)
#ax1.plot(Tpu_wod,ZZi,'-', color=clrW_pct)

ln3, = ax1.plot(Tmn_mom, ZZi, '-', linewidth=2, color=clrM_mn, label="MOM6")
ln4, = ax1.plot(Tmn_gofs, ZZi, '-', linewidth=2, color=clrG_mn, label="GOFS")
ln5, = ax1.plot(Tmn_rtofs, ZZi, '-', linewidth=2, color=clrR_mn, label="RTOFS")

ax1.set_ylim([zlim, 0])
ax1.grid('on')
ax1.set_title(sttl)

sttl = f'Models vs WOD S {YR1}-{YR2} {regn}'
ax2 = plt.axes([0.58, 0.15, 0.4, 0.8])
ax2.plot(Smn_wod, ZZi, '-', linewidth=2, color=clrW_mn) 
#ax2.plot(Spl_wod, ZZi, '-', color=clrW_pct)
#ax2.plot(Spu_wod, ZZi, '-', color=clrW_pct)
spoly = Polygon(sverts, facecolor='0.6', edgecolor='none')
ax2.add_patch(spoly)

ax2.plot(Smn_mom, ZZi, '-', linewidth=2, color=clrM_mn)
ax2.plot(Smn_gofs, ZZi, '-', linewidth=2, color=clrG_mn)
ax2.plot(Smn_rtofs, ZZi, '-', linewidth=2, color=clrR_mn)

ax2.set_ylim([zlim, 0])
ax2.grid('on')
ax2.set_title(sttl)

ax3 = plt.axes([0.6, 0.01, 0.18, 0.12])
lgd = plt.legend(handles=[ln1,ln2], loc='lower left')
ax3.axis('off')

ax4 = plt.axes([0.8, 0.01, 0.18, 0.12])
lgd = plt.legend(handles=[ln3,ln4,ln5], loc='lower left')
ax4.axis('off')

bottom_text(btx, pos=[0.01, 0.01], fsz=8)




