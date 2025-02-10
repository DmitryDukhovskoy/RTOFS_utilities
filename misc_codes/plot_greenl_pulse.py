"""
  Plot Greenland water response to a pulse of freshwater
  from my paper "Timescales of greenland freshwater anomalies ..."
"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
import random


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

random.seed

Time = np.arange(100)
# Bump function:
mua = 20
sgma = 2
muA = 2000   # mean FW pulse amplitude
sgmA = 10

k = 1/15
VFW = []
icc = 0
nruns = 200
for ii in range(nruns):
  aa =  random.normalvariate(mua, sgma)
  dltT = random.uniform(1,10)  # the duration of the freshwater pulse
  bb = aa+dltT

  A = random.normalvariate(muA, sgmA)
  if muA < 50:
    muA = 50

# Heaviside step fn:
  uta = np.where(Time < aa, 0, 1)
  utb = np.where(Time < bb, 0, 1)
  
  vfw = A/k*(uta*(1.-np.exp(-k*(Time-aa))) - utb*(1.-np.exp(-k*(Time-bb))))

  vfw = np.expand_dims(vfw, axis=0)
  if icc == 0:
    VFW = vfw.copy()
  else:
    VFW = np.append(VFW, vfw, axis=0)

  icc += 1


vfw_mn = np.nanmean(VFW, axis=0)

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])

for ii in range(nruns):
  ax1.plot(Time, VFW[ii,:])

ax1.plot(Time, vfw_mn, '-', linewidth=3, color=[0,0,0])

ax1.set_xticks([x for x in range(0,100,10)])
ax1.grid('on')

ax1.set_xlabel('Time, yrs')
ax1.set_ylabel('Freshwater volume, km3')

ax1.set_title('Freshwater volume in the SPNA related to FW pulses')

from mod_utils_fig import bottom_text
btx = 'plot_greenl_pulse.py'
bottom_text(btx) 













