"""
  Statistics of Vol Fluxes from 0.08 expt
  From several years of Vol Fluxes

  fluxes are extracted in 
  /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/anls_fluxes_straits/
  extr_TSVdaily_straits08.m and extr_TSVdaily_straits04.m
  binary for python saved in 

  /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/anls_Fram2/save_VolFWTflux008_112.m
  save_VolFWTflux004.m
  
"""
import os
import numpy as np 
from copy import copy
import importlib
import matplotlib.pyplot as plt
import sys

sys.path.append('/home/ddmitry/codes/MyPython/hycom_utils')
sys.path.append('/home/ddmitry/codes/MyPython/draw_map')
sys.path.append('/home/ddmitry/codes/MyPython')

from mod_utils_fig import bottom_text
plt.ion()  # enables interactive mode

import mod_polynom as mpol
importlib.reload(mpol)

FlxNm    = 'VolFlux'  #  FWFlux or VolFlux
strnm    = 'FramStr'
pthout08 = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/'
#pthout04 = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.04/data_straits/'
#pthout04 = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/023/strait_fluxes/'
expt8 = 112 # 123 - new HYCOM-CICE5, or 112 - old run

# Read section grid, bottom depth
# see: anls_Fram2/save_sect_grid.m
fflux = '{0}hycom008_{1:03d}_{2}_sect_grid.dat'.format(pthout08,expt8,strnm)
print('Reading section grid: ' + fflux)
Hbtm, Long, Latd = mpol.read_sect_grid(fflux)


VFLX  = []
FWFLX = []
yr1   = 2005
yr2   = 2019
for yr in range(yr1,yr2+1):
  dL, Vflx, Tflx, FWFlx = mpol.read_vflux(yr,expt8,strnm,pthout08)
  Vflx  = np.expand_dims(Vflx, axis=0)
  FWFlx = np.expand_dims(FWFlx, axis=0)

  if len(VFLX) == 0:
    VFLX  = Vflx
    FWFLX = FWFlx
  else:
    VFLX  = np.append(VFLX, Vflx, axis=0)
    FWFLX = np.append(FWFLX, FWFlx, axis=0)


# Mean flux:
alf    = 10.
mVF    = np.nanmean(VFLX, axis=0)
pLVF   = np.percentile(VFLX, alf, axis=0)
pUVF   = np.percentile(VFLX, (100.-alf), axis=0)
mFW    = np.nanmean(FWFLX, axis=0)
pLFW   = np.percentile(FWFLX, alf, axis=0)
pUFW   = np.percentile(FWFLX, (100.-alf), axis=0)
#
# Total transport:
VFtot = np.nansum(VFlx, axis=1)
mnVF  = np.mean(VFtot)
stdVF = np.std(VFtot)


#  Plot final result:
clr1 = [0.9,0.5,0.]
clr2 = [0.,0.9,0.1]
clr3 = [0,0.2,0.8]
clr4 = [0.8,0,0.7]
clr5 = [0.,0.7,1.0]
clr6 = [1,0.95,0.2]


xl1 = XX8[0]
xl2 = np.ceil(XX8[-1])

plt.figure(1, figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.53, 0.8, 0.35])


ax1.plot(XX8, Flx8, '-',   color=clr1, linewidth=3, label='HYCOM08')
ax1.plot(XX4, Flx4, '-',   color=clr2, linewidth=3, label='HYCOM04')
ax1.set_xlim(xl1,xl2)

ax1.legend(loc='lower left')
ax1.grid(True)
#plt.grid()
if FlxNm=='VolFlux':
  sftot = 'Vol Flux 0.04 / 0.08 = {0:8.2f} / {1:8.2f} Sv'.format(Ftot4, Ftot8)
else:
  sftot = 'FW Flux 0.04 / 0.08 = {0:8.2f} / {1:8.2f}  mSv'.format(Ftot4*1e-3, Ftot8*1e-3)

ctl = '0.08-{2}, 0.04-{3:03d}, {0} {1} '.format(sftot,yr,expt8, expt4)
plt.title(ctl)

btx='plot_VolFlux.py'
bottom_text(btx)

