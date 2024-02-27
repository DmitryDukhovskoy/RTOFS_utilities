"""
 Find best locations for moorings in the Fram strait
 for annual FW flux
 i.e. find best node locations for
 picewise continuous 1st degree polynomial 
 for observed depth-integrated transport
 interpolation given N mooring locations

 Use equidistant and minmax approaches 
"""
import os
import numpy as np
from copy import copy
import importlib
import matplotlib.pyplot as plt
import sys
import pickle

sys.path.append('/home/ddmitry/codes/MyPython/hycom_utils')
sys.path.append('/home/ddmitry/codes/MyPython/draw_map')
sys.path.append('/home/ddmitry/codes/MyPython')

plt.close('all')
rVFlx = True   # rVFlx - read/plot Vol Flux, otherwise - FW flux
#tplt = 5  

from mod_utils_fig import bottom_text
plt.ion()  # enables interactive mode

import mod_polynom as mpol
importlib.reload(mpol)

Mplnm = 'MMXI' # method to find best polynom nodes
strnm = 'FramStr'
pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/'

expt = 112
n1 = 6    # min # of nodes
n2 = 19   # max # of nodes

btx = 'obsloc_Fram_FWyr.py'


# Order of Polynom
Np = 1; 

PNOD = mpol.poly_stat()
for yr in range (2005,2020):

  dL, Vflx, Tflx, FWflx = mpol.read_vflux(yr, expt, strnm, pthout)

  if rVFlx:
    Flx = Vflx
    FlxNm = 'VolFlux'
  else:
    Flx = FWflx
    FlxNm = 'FWFlux'
# Construct distance array in 100 km
  XX = np.cumsum(dL*1.e-3) 
  XX = XX-dL[0]*1.e-3

# Equidistant points:
  Xeqd = XX

# Control - regularly distributed points:
  nP = XX.shape[0]
  ERR2nrm = []


# Exact solution:
  Yex = Flx

  for N in range(n1,n2+1):
    print('Polynomials with {0} nodes'.format(N))

    if Mplnm == 'EQDST':
  # Polynomial with equidistant nodes:
      PlagrE,XpE,iXpE = mpol.poly1eqd(XX,Flx,N)
      ErrSgm, ErrGlb, ERR2 = mpol.norm2err(XpE,Flx,iXpE,PlagrE)
      Xnds = XpE
      Inds = iXpE
    elif Mplnm == 'MXERR':
  # Polynom 1dgr with nodes adjusted using max error
      Plagr1,Xp,iXp = mpol.poly1maxerr(XX,Flx,N)     
      ErrSgm, ErrGlb, ERR2 = mpol.norm2err(Xp,Flx,iXp,Plagr1)
      Xnds = Xp
      Inds = iXp
    elif Mplnm == 'MINGLERR':
  # Polynom 1dgr with nodes adjusted minimizing global error
  #  selecting the location giving the smallest global error
      P1GE,XpGE,iXpGE = mpol.poly1minErrGlb(XX,Flx,N)     
      ErrSgm, ErrGlb, ERR2 = mpol.norm2err(XpGE,Flx,iXpGE,P1GE)
      Xnds = XpGE
      Inds = iXpGE
    elif Mplnm == 'MINMAX': 
  # Minmax polynom - start with points at max errors and
  # move them around to minimize GlbErr
      P1MMX,XpMMX,iXpMMX = mpol.poly1minmaxErrGlb(XX,Flx,N)
      ErrSgm, ErrGlb, ERR2 = mpol.norm2err(XpMMX,Flx,iXpMMX,P1MMX)
      Xnds = XpMMX
      Inds = iXpMMX
    elif Mplnm == 'MMXI':
# Minmax polynom - start with points at max errors and
# iteratively find locations that give smallest global error
      P1MMXI,XpMMXI,iXpMMXI = mpol.poly1minErrGlb(XX,Flx,N, fguess='maxerr')
      ErrSgm, ErrGlb, ERR2 = mpol.norm2err(XpMMXI,Flx,iXpMMXI,P1MMXI)
      Xnds = XpMMXI
      Inds = iXpMMXI

    PNOD.add_data(yr, N, Xnds, Inds, ErrGlb)


#
# Save polyn. nodes:
#foutp = pthout + 'polynodes_' + strnm + '_' + FlxNm + '_' + Mplnm + '.pkl'
foutp = ('{0}polynodes_{1}_{2}_{3}_hycom008-{4}.pkl'.format(pthout,strnm,FlxNm,Mplnm,expt))
print('Saving ---> ' + foutp)
with open(foutp,'wb') as fid:
  pickle.dump(PNOD, fid)

print(' All Done ')


