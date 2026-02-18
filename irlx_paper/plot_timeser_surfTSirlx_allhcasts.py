"""
  Plot timeseries of monthly surface T/S
  inside the relaxation zone
  to check for the drifts due to ice relaxation

  from PHYS and BGC + IRLX h/casts - both with GLORYS nudging
  extracted in timeser_surfTSicerlx_hind.py

  from PHYS no GLORYS but there is IRLX
  extracted in plot_timeser_surfTSirlx_phys_nonudg.py

  NOAA NWS EMC Dmitry Dukhovskoy
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib
import xarray
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
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hausdorff')

import mod_time as mtime

EXPTS = ['NEPbgc_nudged_hindcast02',
         'NEPphys_nudged_hindcast',
         'NEPphys_nonudg_irlx_hcast']
nexp = len(EXPTS)

EXPTS_INFO=['NEPbgc GLORYS IRLX12hrs',
            'NEPphys GLORYS no IRLX',
            'NEPphys no GLORYS, IRLX24hrs']


def get_timeyr(TM):
  TM = np.array(TM)
  DV = mtime.datevec1D(TM, fHR=False)
  DV = np.array(DV).transpose()
  time_yrs = DV[:,0]+(DV[:,1]-1)/12

  return time_yrs

#rlxt = 12
regn_name = 'NEP10k Arctic Chukchi'
#expt_nameB = 'NEPbgc_nudged_hindcast02'
#expt_nameP = 'NEPphys_nudged_hindcast'

YRS=1993
YRE=2024
xtck = [x for x in range(YRS,YRE+1)]

def read_output(flout):
  pthoutp = '/work/Dmitry.Dukhovskoy/anls_output/NEPbgc_hindcast02'
  dflout = os.path.join(pthoutp,flout)
  print(f'Reading  {dflout}')
  DTM = np.load(dflout)
  SFLD = DTM['SFLD']
  TFLD = DTM['TFLD']
  TM   = DTM['TM']

  return TM, SFLD, TFLD

CLR = [[0.,0.2,.8],
       [0.9,0.3,0.0],
       [0.,0.8,0.3],
       [0.,0.9,0.2],
       [0.8,0.,0.5],
       [0.7, 1, 0.2]]


plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.6, 0.8, 0.35])
ax2 = plt.axes([0.1, 0.15, 0.8, 0.35])

hndls = []
for iexp in range(nexp):
  expt_name = EXPTS[iexp]
  floutP = f'{expt_name}_surfTS_NEParct.npz'
  TM, SFLD, TFLD = read_output(floutP)
  tyrs = get_timeyr(TM)
  clr1 = CLR[iexp]

  expt_info = EXPTS_INFO[iexp]
  if iexp==0:
    lwd=3.5
  elif iexp==1:
    lwd=2.5
  else:
    lwd=1.5

  ln1, = ax1.plot(tyrs, SFLD, '-', linewidth=lwd, color=clr1, label=expt_info)
  hndls.append(ln1)
  # Plot T
  ax2.plot(tyrs, TFLD, '-', linewidth=lwd, color=clr1)


ax1.set_xticks(xtck)
ax1.set_xlim([YRS,YRE+1])
ax1.grid('on')
sttl = f'Hcasts SSS spat.avrg, {regn_name}'
ax1.set_title(sttl)

ax2.set_xticks(xtck)
ax2.set_xlim([YRS,YRE+1])
ax2.grid('on')
sttl = f'Hcasts SST spat.avrg, {regn_name}'
ax2.set_title(sttl)

# Legend
ax2 = plt.axes([0.63, 0.01, 0.35, 0.08])
ax2.legend(handles=hndls, loc='lower right')
ax2.axis('off')



