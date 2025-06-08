import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
#import pickle
import matplotlib.pyplot as plt
from yaml import safe_load

import mod_utils_ob as mutob
importlib.reload(mutob)


PPTHN = []
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')
sys.path.append('./seasonal-workflow')
import mod_time as mtime
import mod_mom6 as mmom6
from mod_utils_fig import bottom_text


# Plot different relaxation time scales
# For sea ice, hours
t=np.arange(100)
TAU = [1/1., 1/2., 1/4., 1/12., 1/120.]

y_init=1.

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))

fig1.clf()
ax1  = plt.axes([0.1, 0.3, 0.8, 0.6])
LNS = []
for n in range(len(TAU)):
  tau = TAU[n]
  dltY = y_init*np.exp(-t*tau)
  ln1, = ax1.plot(t, dltY, linewidth=2, label=f'1/tau={1/tau:.1f} hrs')
  LNS.append(ln1)

sttl = f'Relaxation for different e-folding time scales'
ax1.grid('on')
ax1.set_title(sttl)
ax1.set_xlabel('Time, hours')
ax1.set_xlim([0,np.max(t)])
ax1.xaxis.set_ticks(list(np.arange(0,np.max(t),10)))

ax1.plot([0,np.max(t)],[np.exp(-1),np.exp(-1)],'--',color=[0.6,0.6,0.6])

ax3 = plt.axes([0.55, 0.04, 0.4, 0.2])
lgd = plt.legend(handles=LNS, loc='upper right')
ax3.axis('off')

btx='diagr_relax_time_scale.py'
bottom_text(btx)



