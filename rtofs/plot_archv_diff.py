"""
  Plot diagnostic output fields from
  hycom_diff
  that produces difference between the
  archive inc (increments added) - initial (background)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
import importlib
import struct
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
#from mpl_toolkits.basemap import Basemap, shiftgrid

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
import mod_read_ncoda as rncoda

#rdate   = '20221117'
rdate   = '20220616'
pthbs   = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/'
pth1    = pthbs+'hycom/incup/'
pthbin  = pth1
pthbgr  = pthbs+'hycom/ncoda_archv_inc/'  # background field dir
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

fldplt = 'v-vel.'  # thknss temp salin u-vel. v-vel.
jday = rncoda.rdate2julian(rdate)
yr, mo, mday, hr = rncoda.parse_rdate(rdate)

# fldiff - Output file with differences btw b/ground and analysis
# flbgr  - background field from HYCOM f/cast
fldiff  = 'incupd.2022_{0:03d}_18'.format(jday-1)
flbgr   = 'archv.2022_{0:03d}_00'.format(jday) 

IDM = 4500
JDM = 3298
KDM = 41

ftopo = 'regional.depth'
fgrid = 'regional.grid'

import mod_read_hycom
#importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo, read_hycom, \
                           read_topo
import mod_utils as utls
#importlib.reload(utls)
from mod_read_hycom import read_2Dbinary
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text

btx = 'plot_archv_diff.py'

#  HH = read_topo(pthgrid,ftopo,nn,mm)
LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)

huge = 1.e20
rg   = 9806.

# Read layer depths from the background field
fina = pthbgr+flbgr+'.a'
finb = pthbgr+flbgr+'.b'
print('Reading background: '+fina)
fld  = 'thknss'
F,nn,mm,ll = read_hycom(fina,finb,fld,rLayer=1)
F[np.where(F>huge)] = np.nan
F = F/rg
F[np.where(F<0.001)] = 0.

dH = np.zeros((ll,mm,nn))
dH[0,:,:] = F
for kk in range(2,ll+1):
  F,nn,mm,lmm = read_hycom(fina,finb,fld,rLayer=kk)
  F = F/rg
  F[np.where(F>huge)] = np.nan
  F[np.where(F<0.001)] = 0.
  dH[kk-1,:,:] = F

ZZ, ZM = zz_zm_fromDP(dH, f_btm=False)
kdm = ZM.shape[0]
jdm = ZM.shape[1]
idm = ZM.shape[2]

# ===========================
# Read incup diff files
#
fina = pthbin+fldiff+'.a'
finb = pthbin+fldiff+'.b'
print('Processing '+fina)

fld = fldplt
A3d = np.array([])
for kk in range (1,KDM+1):
  F,n1,m1,l1 = read_hycom(fina,finb,fld,rLayer=kk)
  F[np.where(F>huge)] = np.nan
  if fld == "thknss":
    F = F/rg

  if A3d.size == 0:
    A3d = F.copy()
    A3d = np.expand_dims(A3d, axis=0)
  else:
    F = np.expand_dims(F, axis=0)  # Pa -> m
    A3d = np.append(A3d, F, axis=0)

if fldplt == "thknss":
  rmin = -50.
  rmax = 50.
  eps0 = 500.
  rtitle = 'dPanls-dPbgr, m'
elif fldplt == "temp":
  rmin = -1.
  rmax = 1.
  eps0 = 3.
  rtitle = 'Tanls-Tbgr, C'
elif fldplt == "salin":
  rmin = -0.5
  rmax = 0.5
  eps0 = 1.
  rtitle = 'Sanls-Sbgr, psu'
elif fldplt == "u-vel.":
  rmin = -0.25
  rmax = 0.25
  eps0 = 0.5
  rtitle = 'Uanls-Ubgr, m/s'
elif fldplt == "v-vel.":
  rmin = -0.25
  rmax = 0.25
  eps0 = 0.5
  rtitle = 'Vanls-Vbgr, m/s'

# Check max/min
# Plot maps of max/min incup diff
import plot_sect as psct
importlib.reload(psct)

f_pltmax = True
if f_pltmax:
  import mod_utils as mutil
  importlib.reload(mutil) 

  clrmax = psct.clrmp_WhYlOrRd()
  clrmin = psct.clrmp_BlGrWh()
  
  maxA = mutil.find_max3d(A3d,dH)
  minA = mutil.find_min3d(A3d,dH)

  a1,a2 = np.where(maxA==np.nanmax(maxA))
  jmax = a1[0]
  imax = a2[0]
  a1,a2 = np.where(minA==np.nanmin(minA))
  jmin = a1[0]
  imin = a2[0]

  ctl = 'incup {0} {1} eps={2:4.1f} rdate={3} max({5},{6})={4:5.1f}'.\
         format(rtitle,fldiff,eps0,rdate,np.nanmax(maxA),imax,jmax) 
  fgnmb = 20
  J,I = np.where(maxA >= eps0)
  psct.plot_2Dmap(maxA,HH,fgnmb=fgnmb,clrmp=clrmax,rmin=0.,rmax=rmax,\
               ctl=ctl,btx=btx,Ipnt=I,Jpnt=J)


  ctl = 'incup {0} {1} eps={2:4.1f} rdate={3} min({5},{6})={4:5.1f}'.\
         format(rtitle,fldiff,eps0,rdate,np.nanmin(minA),imin,jmin)
  fgnmb = 21
  J,I = np.where(minA <= -eps0)
  psct.plot_2Dmap(minA,HH,fgnmb=fgnmb,clrmp=clrmin,rmin=rmin,rmax=0.,\
               ctl=ctl,btx=btx,Ipnt=I,Jpnt=J)

#KK,JJ,II = np.where(A3d > eps0)
# ==========================
# Plot selected W-E profiles  
# Plot vertical section of background and increm rho
# and difference btw the two
SCTP = ["GoM1"]
#SCTP = ["Test2"]
#SCTP = ["SeaJpn"]
nsct = len(SCTP)

plt.ion()
nfg = 0
#f_plt = True
#if f_plt:
for ii in range(nsct):
  nfg += 1
  xsct = SCTP[ii]
  SCT = psct.Jsections()
  aa = SCT.get(xsct)
  i1 = aa[0]
  i2 = aa[1]
  j1 = aa[2]


# find index to print out layers:
# ipp,jpp - global indices
# IP0, JP0 - local wrt to section, for plotting
# JP0 = 0 and not needed for E-W section case
  ip0 = -1
  ipp,jpp,IP0,JP0 = psct.find_indx_lonlat(-170.42,-14.25,LON,LAT,xsct=xsct)
  if xsct == "GoM1":
    ip0 = 72  # GoM corresponds the itest, jtest in ncoda_archv_inc
  elif xsct == "SeaJpn":
    ip0 = 56
  elif xsct == "Test1":
    ip0 = 102
  elif xsct == "Test2":
    ip0 = 154

  clrmp = copy(plt.cm.BrBG_r)
  stl = ('incup: {3}, {0} {1}, rdate={2}'.\
         format(xsct,pthbin[-13:]+fldiff,rdate,rtitle))
  Hb,XX,Asct,ZZs = psct.plot_rho_Jsct(xsct, A3d, ZZ, HH, LON, LAT,\
                    fgnmb=nfg, stl=stl, sct_show=True, btx=btx, \
                    rmin=rmin, rmax=rmax, clrmp=clrmp,jlrs=ip0)

# Diagnostics for section:
  amax = np.nanmax(Asct)
  amin = np.nanmin(Asct)
  stxt = 'min/max = {0:10.4f}/{1:10.4f}\n'.format(amin,amax) 
  stxt = stxt+'section i={0}:{1}, j={2}'.format(i1,i2,j1)
  ax3 = plt.axes([0.1, 0.17, 0.3, 0.05])
  ax3.text(0.1,0.1,stxt)
  ax3.axis('off')

#  print_1D(dZZs[:,ip0])
  if ip0 > 0:
    psct.print_2D(ZZs[:,ip0],Asct[:,ip0],kend=kdm)
 
  f_3d = False
  if f_3d:
#    imm = 1454
#    jmm = 1324
#    imm = 1372
#    jmm = 1336
    imm = 2954
    jmm = 782
    dz  = dH[:,jmm,imm]
    aa  = A3d[:,jmm,imm]
    psct.print_3D(dz,aa,np.cumsum(dz),wd=14,prc=6)


