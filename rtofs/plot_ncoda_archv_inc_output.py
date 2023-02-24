"""
  Plot output fields from
  ncoda_archv_inc.f --> archv_1_inc.${archday}.[a,b]
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

sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
import mod_read_ncoda as rncoda

#rdate   = '20221117'
rdate   = '20220616'
pthbs   = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/'
pth1    = pthbs+'hycom/ncoda_archv_inc/'
pthbin  = pth1
pthbgr  = pthbs+'hycom/ncoda_archv_inc/'  # background field dir
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
pthrsrt = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/rtofs_expts/'+\
          'rtofs_para8.1/navy_restart/'

fldplt = 'salin'  # temp salin u-vel. v-vel.
jday = rncoda.rdate2julian(rdate)
yr, mo, mday, hr = rncoda.parse_rdate(rdate)

# Select - updated (True) or background (False) fields:
f_updated = True
# flupdt - updated fields=background+increments output from ncoda_archv_inc
# flbgr  - background field from HYCOM f/cast
flupdt = 'archv_1_inc.{0}_{1:03d}_00'.format(yr,jday)
flbgr   = 'archv.{0}_{1:03d}_00'.format(yr,jday) 

IDM = 4500
JDM = 3298
KDM = 41

get_topo = True
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

btx = 'plot_ncoda_archv_inc_output.py'

#  HH = read_topo(pthgrid,ftopo,nn,mm)
LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)

huge = 1.e20
rg   = 9806.

# Read layer depths
if f_updated:
# Updated fields with increments:
  fina = pthbin+flupdt+'.a'
  finb = pthbin+flupdt+'.b'
else:
# Background fields:
  fina = pthbin+flbgr+'.a'
  finb = pthbin+flbgr+'.b'


print('Reading layer thickness: '+fina)
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
    F = np.expand_dims(F, axis=0)
    A3d = np.append(A3d, F, axis=0)

# Check max/min
#IJ = np.where(A3d >= 2000.)
amin = np.nanmin(A3d)
amax = np.nanmax(A3d)
print('Fld={0} min/max = {1:8.4f}/{2:8.4f}'.format(fldplt,amin,amax))

# ==========================
# Plot selected W-E profiles  
# Plot vertical section of background and increm rho
# and difference btw the two
SCTP = ["GoM1"]
#SCTP = ["GoM1","SeaJpn"]
nsct = len(SCTP)

plt.ion()
import plot_sect as psct
importlib.reload(psct)

clrmp = copy(plt.cm.BrBG_r)
if fldplt == "temp":
  clrmp = copy(plt.cm.rainbow)
#  rmin = -2.
#  rmax = 24.
elif fldplt == "salin":
  clrmp = copy(plt.cm.gist_ncar_r)
#  clrmp = copy(plt.cm.rainbow)
#  rmin = 34.5
#  rmax = 36.5
elif fldplt == 'u-vel.' or fldplt == 'v-vel.':
  clrmp = copy(plt.cm.seismic)

clrmp.set_bad(color=[0.3,0.3,0.3])

jlrs = -1
nfg = 0
for ii in range(nsct):
  nfg += 1
  xsct = SCTP[ii]
  SCT = psct.Jsections()
  aa = SCT.get(xsct)
  i1 = aa[0]
  i2 = aa[1]
  j1 = aa[2]

  dmm = A3d[:,j1,i1:i2].flatten()
  dmm = dmm[~np.isnan(dmm)]

  if fldplt == 'u-vel.' or fldplt == 'v-vel.':
    rmin, rmax = psct.minmax_clrmap(dmm,pmin=5.,pmax=95.,cpnt=0.01,fsym=True)
  else:
    rmin, rmax = psct.minmax_clrmap(dmm,pmin=10.,pmax=90.,cpnt=0.01)

  if xsct == "GoM1":
    jp0 = 72  # GoM corresponds the itest, jtest in ncoda_archv_inc
  elif xsct == "SeaJpn":
    jp0 = 55

  if f_updated:
    cfls = 'Updated fld'
  else:
    cfls = 'B/ground fld'

  stl = ('{3}, {0} {4} {1}, rdate={2}'.\
         format(xsct,fina[-35:],rdate,fldplt,cfls))
  Hb, Xsct, Zsct, Asct = psct.plot_A3d_Jsct(xsct, A3d, ZZ, HH, LON, LAT,\
                           fgnmb=nfg, stl=stl, sct_show=True, btx=btx, \
                           rmin=rmin, rmax=rmax, clrmp=clrmp, jlrs=jp0)

  if jp0>0:
    psct.print_2D(Zsct[:,jp0],Asct[:,jp0],kend=KDM)

  f_3d = False
  if f_3d:
    imm = 1454
    jmm = 1324
#    imm = 1372
#    jmm = 1336
#    imm = 2954
#    jmm = 782
    dz  = dH[:,jmm,imm]
    aa  = A3d[:,jmm,imm]
    psct.print_3D(dz,aa,np.cumsum(dz),wd=14,prc=6)



