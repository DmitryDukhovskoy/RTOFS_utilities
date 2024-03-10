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
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
import yaml

PTHR = '/Users/ddmitry/python'
sys.path.append(PTHR+'/MyPython/hycom_utils')
sys.path.append(PTHR+'/MyPython/draw_map')
sys.path.append(PTHR+'/MyPython')
sys.path.append(PTHR+'/MyPython/mom6_utils')

import mod_time as mtime
from mod_utils_fig import bottom_text
import mod_polynom as mpol
importlib.reload(mpol)
import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)
import mod_misc1 as mmisc
import mod_mom6 as mom6util
import mod_colormaps as mcmp
import mod_read_hycom as mhycom
importlib.reload(mom6util)
importlib.reload(mmisc)
importlib.reload(mcmp)


plt.close('all')
rVFlx = True   # rVFlx - read/plot Vol Flux, otherwise - FW flux
#tplt = 5  

plt.ion()  # enables interactive mode


# min/max of polynom nodes:
n1 = 5
n2 = 18
YRM = 2017

nrun  = 'ARCc0.08'  # MOM6, RTOFS, GOFS3.1, ARCc0.08

#yr = 2019  # yr = 0 - some analytical functions
# Use volume flux estimates computed extrUVFlx_hycom.py 
# COAPS server
yr  = 0
sctnm = 'Fram79s2'
fld2d = 'Unrm'

mS=1
dS=1
if nrun == 'MOM6':
  expt = '003'
  YR   = YRM
elif nrun == 'RTOFS':
  expt = 'product' # 003 or product
  YR   = YRR
#mS=9
#dS=15
elif nrun == 'GOFS3.1':
  expt = '93.0'
  YR   = YRM
elif nrun == 'ARCc0.08':
  expt = '112'
  YR = YRM

dnmb1 = mtime.datenum([YR,mS,dS])
dnmb2 = mtime.datenum([YR,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)

pthoutp = '/Users/ddmitry/DATA/ARCc0.08/data_straits/'
floutp  = f"008arc-{expt}_{fld2d}xsct_{dv1[0]}" + \
        f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
pthgrid = '/Users/ddmitry/DATA/ARCc0.08/topo_grid_arc08/'
ftopo   = 'regional.depth'
fgrid   = 'regional.grid'
_, _, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

ffout = pthoutp + floutp
print('Loading ' + ffout)
  
with open(ffout, 'rb') as fid:
  F2D, UFLX = pickle.load(fid)


# 2D fields are at half-grid points
TM   = F2D.TM
Unrm = F2D.Fld2D  # 2D Flow: Time x depth x Width
II   = F2D.Iindx
JJ   = F2D.Jindx
XX   = F2D.LON
YY   = F2D.LAT
dL   = F2D.Lsgm
Hbtm = F2D.Hbtm
ZZi  = F2D.ZZi
# Depth-integrated: full grid
VFlx = UFLX.trnsp*1e-6  # 1D depth-integrated flow, Sv

#  dL, VFlx, Tflx, FWflx = mpol.read_vflux(yr,strnm,pthout)


if rVFlx:
  Flx = np.nanmean(VFlx, axis=0)  # time mean V flux
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






