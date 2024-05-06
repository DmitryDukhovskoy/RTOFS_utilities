import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
#import struct
import datetime
import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
import timeit
import yaml
from netCDF4 import Dataset as ncFile

#PPTHN = '/home/Dmitry.Dukhovskoy/python'
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

import mod_interp1D as minterp
importlib.reload(minterp)

def func(xx):
  yy = np.exp(-xx/1.5)*np.cos(xx**3)
  return yy

XX = np.arange(0,3.14,0.01)
#YY = np.exp(-XX/1.5)*np.cos(XX**3)
YY = func(XX)

Ip = list(k for k in range(100,290,18))
Xp = XX[Ip]
Yp = YY[Ip]

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
ax1.plot(XX,YY)
ax1.plot(Xp,Yp,'o')

XI  = np.arange(Xp[0],Xp[-1]-1.e-3,0.01)
nn  = len(XI)
Yi  = np.zeros((nn))
Y1i = np.zeros((nn))
for ik in range(nn):
  xx = XI[ik]
  Pn = minterp.pcws_lagr3(Xp,Yp,xx)
  Yi[ik] = Pn
  P1 = minterp.pcws_lagr1(Xp,Yp,xx)
  Y1i[ik] = P1

ax1.plot(XI,Yi,'.-')
ax1.plot(XI,Y1i,'.-')
ax1.set_title('Cubic and linear interpolation of f(x) for given nodes')

from mod_utils_fig import bottom_text
btx = 'test_interp.py'
bottom_text(btx)

# 2D field:
# Interpolate along columns:
nnd=4
X2 = np.column_stack((Xp[:nnd],Xp[1:nnd+1],Xp[2:nnd+2],Xp[3:nnd+3],Xp[4:nnd+4],\
                      Xp[5:nnd+5]))
a1   = np.linspace(1,4,4)
a2   = np.linspace(1,6,6)
#_,X2 = np.meshgrid(a2,a1)
Y2 = np.column_stack((Yp[:nnd],Yp[1:nnd+1],Yp[2:nnd+2],Yp[3:nnd+3],Yp[4:nnd+4],\
                      Yp[5:nnd+5]))

#xx = 2.75*np.ones((Y2.shape[1]))
xx = 0.57*(X2[-1,:] - X2[0,:])+X2[0,:]
P2d = minterp.lagr_polynom2D(X2, Y2, xx)
P1d = minterp.lagr_polynom(X2[:,-1],Y2[:,-1],xx[-1])
Pp3 = minterp.pcws_lagr3(Xp,Yp,xx[-1])

ax1.plot(xx, P2d, 'o')

#ax2 = plt.axes([0.1, 0.08, 0.8, 0.38])
#ax2.plot(X2,Y2,'.-')
#ax2.plot(xx,P2d,'o')







