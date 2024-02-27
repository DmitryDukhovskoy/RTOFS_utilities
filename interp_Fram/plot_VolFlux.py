"""
  Plot Vol Fluxes from 0.08 and 0.04 simualtions
  fluxes are extracted in 
  /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/anls_fluxes_straits/
  extr_TSVdaily_straits08.m and extr_TSVdaily_straits04.m
  binary for python saved in 
  
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
pthout04 = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.04/data_straits/'
#pthout04 = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/023/strait_fluxes/'
expt8 = 112 # 123 - new HYCOM-CICE5, or 112 - old run
expt4 = 23


yr = 2017
dL8, Vflx8, Tflx8, FWFlx8 = mpol.read_vflux(yr,expt8,strnm,pthout08)

XX8 = np.cumsum(dL8*1.e-3) 
XX8 = XX8-dL8[0]*1.e-3

if FlxNm == 'VolFlux':
  Flx8 = Vflx8
else:
  Flx8 = FWFlx8

Ftot8 = np.nansum(Flx8)


dL4, Vflx4, Tflx4, FWFlx4 = mpol.read_vflux(yr,expt4,strnm,pthout04, res='004')

XX4 = np.cumsum(dL4*1.e-3)
XX4 = XX4-dL4[0]*1.e-3

if FlxNm == 'VolFlux':
  Flx4 = Vflx4
else:
  Flx4 = FWFlx4

Ftot4 = np.nansum(Flx4)


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

