"""
  Plot
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

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text
import mod_mom6_valid as mom6vld
import mod_misc1 as mmisc
import mod_time as mtime
import mod_WODdata as mwod
importlib.reload(mwod)
importlib.reload(mom6vld)


# WOD data
pthwod = '/scratch1/NCEPDEV/stmp4/Dmitry.Dukhovskoy/WOD_profiles/'
pthoutp= '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/MOM6_CICE6/ts_prof/'

modrun = 'rtofs'
YR1    = 2018
YR2    = 2022
mo1    = 1
mo2    = 12
regn   = 'AmundsAO'
#regn   = 'NansenAO'
#regn   = 'MakarovAO'
#regn  = 'CanadaAO'

REGNS = mom6vld.ts_prof_regions()
x0 = REGNS[regn]["lon0"]
y0 = REGNS[regn]["lat0"] 
dx = REGNS[regn]["dlon"]
dy = REGNS[regn]["dlat"]

print(f'Plotting WOD T/S for selected UID {regn} {YR1} {YR2} mo={mo1}-{mo2}\n')

#flts_out = f'WODTS_{regn}_{YR1}-{YR2}.pkl'
flts_out = f'{modrun}TS_{regn}_{YR1}-{YR2}.pkl'
dfltsout = os.path.join(pthoutp, flts_out)

# Load T/S profiles
print(f'Loading {dfltsout}')
with open(dfltsout, 'rb') as fid:
  [SPROF, TPROF] = pickle.load(fid)

ZZi   = TPROF.zz
TT    = TPROF.prof1d
SS    = SPROF.prof1d

# Average and find percentiles
plow = 5.
pup  = 100-plow
Tmn  = np.nanmean(TT, axis=0)
Tpl  = np.nanpercentile(TT, plow, axis=0)
Tpu  = np.nanpercentile(TT, pup, axis=0)

Smn  = np.nanmean(SS, axis=0)
Spl  = np.nanpercentile(SS, plow, axis=0)
Spu  = np.nanpercentile(SS, pup, axis=0)


btx  = 'plot_modelTSprof_vsWOD.py'

plt.ion()
clr_mn  = [0, 0.3, 0.7]
clr_pct = [.7, 0.7, 0.7]


fig1 = plt.figure(1,figsize=(9,9))
plt.clf()

sttl = f'{modrun} T {YR1}-{YR2} {regn}'
izm = min(np.where((ZZi < -100.) & (np.isnan(Tmn)))[0])
zlim = np.floor(ZZi[izm])
ax1 = plt.axes([0.08, 0.1, 0.4, 0.8])
ln1, = ax1.plot(Tmn, ZZi, '-', linewidth=2, color=clr_mn, label="mean")
ln2, = ax1.plot(Tpl,ZZi,'-', color=clr_pct, label=f"{plow}-{pup}%")
ax1.plot(Tpu,ZZi,'-', color=clr_pct)

ax1.set_ylim([zlim, 0])
ax1.grid('on')
ax1.set_title(sttl)

lgd = plt.legend(handles=[ln1,ln2], loc='lower right')


sttl = f'{modrun} S {YR1}-{YR2} {regn}'
ax2 = plt.axes([0.58, 0.1, 0.4, 0.8])
ax2.plot(Smn, ZZi, '-', linewidth=2, color=clr_mn) 
ax2.plot(Spl, ZZi,'-', color=clr_pct)
ax2.plot(Spu, ZZi,'-', color=clr_pct)

ax2.set_ylim([zlim, 0])
ax2.grid('on')
ax2.set_title(sttl)

bottom_text(btx, pos=[0.02, 0.03])




