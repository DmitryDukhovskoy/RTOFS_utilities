"""
 Find node locations for
 picewise continuous 1st degree polynomial 
 for observed depth-integrated transport
 interpolation given N mooring locations

 Start with equidistant N nodes then 
 adjust the node within the segment with largest L2 error
 i.e. 2-norm
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
#tplt = 5  

from mod_utils_fig import bottom_text
plt.ion()  # enables interactive mode

import mod_polynom as mpol
importlib.reload(mpol)

# min/max of polynom nodes:
n1 = 4
n2 = 18

#yr = 2019  # yr = 0 - some analytical functions
yr  = 0
strnm = 'FramStr'
pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/'

if yr == 0:
  nn = 211
  dL, Vflx, Tflx, FWflx = mpol.anlfn(nn)
else:
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

# Control - regularly distributed points:
nP = XX.shape[0]
ERR2nrm = []

btx = 'best_node_pcws_polynom1.py'

# Exact solution:
Yex = Flx
# 2-norm error
ERR2GLB      = []
ERR2GLB_EQ   = []
ERR2GLB_GE   = []
ERR2GLB_MMX  = []
ERR2GLB_MMXI = []

for N in range(n1,n2+1):
  print('Polynomials with {0} nodes'.format(N))
# Polynomial with equidistant nodes:
  PlagrE,XpE,iXpE = mpol.poly1eqd(XX,Flx,N)
  ErrSgmE, ErrGlbE, ERR2E = mpol.norm2err(XpE,Flx,iXpE,PlagrE)
  ERR2GLB_EQ.append(ErrGlbE)

# Polynom 1dgr with nodes adjusted using max error
  Plagr1,Xp,iXp = mpol.poly1maxerr(XX,Flx,N)     
  ErrSgm, ErrGlb, ERR2 = mpol.norm2err(Xp,Flx,iXp,Plagr1)
  ERR2GLB.append(ErrGlb)    
   
# Polynom 1dgr with nodes adjusted minimizing global error
#  selecting the location giving the smallest global error
# iteratevely
  P1GE,XpGE,iXpGE = mpol.poly1minErrGlb(XX,Flx,N, fguess='equidist')     
  ErrSgmGE, ErrGlbGE, ERR2GE = mpol.norm2err(XpGE,Flx,iXpGE,P1GE)
  ERR2GLB_GE.append(ErrGlbGE)    
   
# Minmax polynom - start with points at max errors and
# move them around to minimize GlbErr
  P1MMX,XpMMX,iXpMMX = mpol.poly1minmaxErrGlb(XX,Flx,N)
  ErrSgmMMX, ErrGlbMMX, ERR2MMX = mpol.norm2err(XpMMX,Flx,iXpMMX,P1MMX)
  ERR2GLB_MMX.append(ErrGlbMMX)

# Minmax polynom - start with points at max errors and
# iteratively find locations that give smallest global error
  P1MMXI,XpMMXI,iXpMMXI = mpol.poly1minErrGlb(XX,Flx,N, fguess='maxerr')
  ErrSgmMMXI, ErrGlbMMXI, ERR2MMXI = mpol.norm2err(XpMMXI,Flx,iXpMMXI,P1MMXI)
  ERR2GLB_MMXI.append(ErrGlbMMXI)



# Total Flux:
Ftot = np.nansum(Flx)

#  Plot final result:
clr1 = [0.9,0.5,0.]
clr2 = [0.,0.9,0.1]
clr3 = [0,0.2,0.8]
clr4 = [0.8,0,0.7]
clr5 = [0.,0.7,1.0]
clr6 = [1,0.95,0.2]

Yp   = Flx[iXp]
YpE  = Flx[iXpE]
YpGE = Flx[iXpGE]

xl1 = XX[0]
xl2 = np.ceil(XX[-1])

plt.figure(1, figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.53, 0.8, 0.35])
ax1.plot(XX, Flx, '-',   color=clr1, linewidth=3, label='HYCOM')
ax1.plot(XX, Plagr1,'-', color=clr2, label='Plnm1 maxErr')
ax1.plot(XX, PlagrE,'-', color=clr3, label='Plnm1 eqdst')
ax1.plot(XX, P1GE,'-',   color=clr4, label='Plnm1 minGErr')
ax1.plot(XX, P1MMX,'-',  color=clr5, label='Plnm1 minmax')
ax1.plot(XX, P1MMXI,'-', color=clr6, label='Plnm1 mmaxItr')
#ax1.plot(Xp,Yp,'g.')
#ax1.plot(XpE,YpE,'b.')
ax1.set_xlim(xl1,xl2)
#ax1.set_ylim(t1l1,t1l2)

ax1.legend(loc='upper right')
ax1.grid(True)
#plt.grid()
if rVFlx:
  sftot = 'Vol Flux={0:8.2f} Sv'.format(Ftot*1e-6)
else:
  sftot = 'FW Flux={0:8.2f} mSv'.format(Ftot*1e-3)

ctl = '{0} {1} Npnts={2}'.format(sftot,yr,N)
plt.title(ctl)

NP = np.arange(n1,n2+1)
ax2 = plt.axes([0.1, 0.1, 0.8, 0.35])
ax2.plot(NP, ERR2GLB,     color=clr2, label='Polyn maxErr')
ax2.plot(NP, ERR2GLB, '.', color=clr2)
ax2.plot(NP, ERR2GLB_EQ,  color=clr3, label='Polyn eqdst')
ax2.plot(NP, ERR2GLB_EQ, '.', color=clr3)
ax2.plot(NP, ERR2GLB_GE,  color=clr4, label='Polyn minGErr')
ax2.plot(NP, ERR2GLB_GE, '.', color=clr4)
ax2.plot(NP, ERR2GLB_MMX, color=clr5, label='Polyn minmax')
ax2.plot(NP, ERR2GLB_MMX, '.', color=clr5)
ax2.plot(NP, ERR2GLB_MMXI, color=clr6, label='Polyn mmaxItr')
ax2.plot(NP, ERR2GLB_MMXI, '.', color=clr6)
ax2.legend(loc='upper right')
ax2.set_xlim(n1-0.2,n2+0.2)
ax2.set_xticks(np.arange(n1,n2+1))
ax2.grid(True)
#plt.grid()
ctl = 'Err2 Global'
plt.title(ctl)

bottom_text(btx)






