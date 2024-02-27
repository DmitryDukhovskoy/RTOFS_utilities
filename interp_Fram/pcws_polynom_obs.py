"""
 Find polynomial for observed depth-integrated transport
 picewise polynomial: try linear, quadr, cubic
 then Hermite polynomial
 interpolation given N mooring locations
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

plt.close('all')
rVFlx = False   # rVFlx - read/plot Vol Flux, otherwise - FW flux

from mod_utils_fig import bottom_text
plt.ion()  # enables interactive mode

import mod_polynom as mpol
importlib.reload(mpol)

yr = 2019
strnm = 'FramStr'
pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/'

dL, Vflx, Tflx, FWflx = mpol.read_vflux(yr,strnm,pthout)

if rVFlx:
  Flx = Vflx
  FlxNm = 'VolFlux'
else:
  Flx = FWflx
  FlxNm = 'FWFlux'

# Construct distance array in 100 km
XX = np.cumsum(dL*1.e-3) 
XX = XX-dL[0]*1.e-3

# Order of Polynom
Np = 1; 

# Equidistant points:
Xeqd = XX

# Start with regularly distributed points:
nP = XX.shape[0]
ERR1inf  = [] # 1st dgr polynom
ERR12nrm = []
ERR2inf  = [] # 2-nd dgr polynom
ERR22nrm = []
ERRHinf  = []
ERRH2nrm = []
n1 = 4
n2 = 40
iplt = 40
if iplt>n2:
 iplt = n2

btx = 'pcws_polynom_obs.py'

clr1 = [0.9,0.4,0.]
clr2 = [0.,0.8,0.2]
clr3 = [1,0.2,0.8]
for N in range(n1,n2+1):
  dx = (nP-1)/(N-1)
  iXp = np.round(np.arange(0,nP-0.1,dx)).astype(int)
  iXp = list(iXp)

  Yp = Flx[iXp]
  Xp = XX[iXp]


  Plagr1 = []  # estimated flux for equidistant polynomial points
  Plagr2 = []
  Phermt = []
  for ip in range(nP):
    yip1 = mpol.pcws_lagr1(Xp,Yp,XX[ip])
    yip2 = mpol.pcws_lagr2(Xp,Yp,XX[ip])
    yip3 = mpol.pcws_hermite(Xp,Yp,XX[ip])
    Plagr1.append((yip1))
    Plagr2.append((yip2))
    Phermt.append((yip3))

  Plagr1 = np.array(Plagr1)
  Plagr2 = np.array(Plagr2)
  Phermt = np.array(Phermt)

  f_plt = True
  if f_plt and N == iplt:
    fig1 = plt.figure(1)
    plt.clf()
    line1, = plt.plot(XX,Flx,label='HYCOM')
    plt.plot(Xp,Yp,'ro')
    line2, = plt.plot(XX,Plagr1,color=clr1,label='Lagr1') #interpolated
    line3, = plt.plot(XX,Plagr2,color=clr2,label='Lagr2')  # interpolated
    line4, = plt.plot(XX,Phermt,color=clr3,label='Hermt')
    plt.legend(handles=[line1,line2,line3,line4])
    plt.grid()
    ctl = '{0} {1}'.format(FlxNm,yr)
    plt.title(ctl)
    bottom_text(btx,pos=[0.05,0.02])

# Exact solution:
  Yex = Flx
  err1_inf = np.max(abs(Yex-Plagr1))/np.max(abs(Yex))
  err2_inf = np.max(abs(Yex-Plagr2))/np.max(abs(Yex))
  errH_inf = np.max(abs(Yex-Phermt))/np.max(abs(Yex))
  ERR1inf.append((err1_inf))
  ERR2inf.append((err2_inf))
  ERRHinf.append((errH_inf))
# 2-norm
  err12nrm = np.sqrt(np.dot((Yex-Plagr1),(Yex-Plagr1)))
  err22nrm = np.sqrt(np.dot((Yex-Plagr2),(Yex-Plagr2)))
  errH2nrm = np.sqrt(np.dot((Yex-Phermt),(Yex-Phermt)))
  ERR12nrm.append((err12nrm))
  ERR22nrm.append((err22nrm))
  ERRH2nrm.append((errH2nrm))

f_plterr = True
plt.ion()
if f_plterr:
  Nx = np.arange(n1,n2+1).astype(int)
  fig2 = plt.figure(2)
  plt.clf()
  line1, = plt.plot(Nx,ERR12nrm,'.-',color=clr1,label='Lagr1')
  line2, = plt.plot(Nx,ERR22nrm,'.-',color=clr2,label='Lagr2')
  line3, = plt.plot(Nx,ERRH2nrm,'.-',color=clr3,label='Hermt')
  plt.legend(handles=[line1,line2,line3])
  plt.xlabel('N polynom global nodes')
  plt.grid()
  ctl = 'Sqrt 2-norm Error {1}'.format(FlxNm,yr)
  plt.title(ctl)
  bottom_text(btx,pos=[0.05,0.02])





