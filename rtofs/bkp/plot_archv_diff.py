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

sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython')

rdate   = '20221117'
pthbs   = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/'
pth1    = pthbs+'hycom/incup/'
pthbin  = pth1
pthbgr  = pthbs+'hycom/ncoda_archv_inc/'  # background field dir
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

fldiff  = 'incupd.2022_320_18'
flbgr   = 'archv.2022_321_00' 


IDM = 4500
JDM = 3298
KDM = 41

# Plot vertical section of background and increm rho
# and difference btw the two
SCTP = ["GoM1"]
nsct = len(SCTP)


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

btx = 'plot_archv_diff.py'

get_topo = True
if get_topo:
#  HH = read_topo(pthgrid,ftopo,nn,mm)
  LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)
  get_topo = False  

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


# Read diff
fina = pthbin+fldiff+'.a'
finb = pthbin+fldiff+'.b'
print('Processing '+fina)

fld = 'thknss'
dP_diff = np.array([])
for kk in range (1,KDM+1):
  F,n1,m1,l1 = read_hycom(fina,finb,fld,rLayer=kk)
  F[np.where(F>huge)] = np.nan
  if dP_diff.size == 0:
    dP_diff = F.copy()
    dP_diff = np.expand_dims(dP_diff, axis=0)
  else:
    F = np.expand_dims(F, axis=0)/rg  # Pa -> m
    dP_diff = np.append(dP_diff, F, axis=0)


# ==========================
# Plot selected W-E profiles  
plt.ion()
import plot_sect as psct
importlib.reload(psct)

nfg = 0
for ii in range(nsct):
  nfg += 1
  xsct = SCTP[ii]
  SCT = psct.Jsections()

  clrmp = copy(plt.cm.BrBG_r)
  stl = ('hycom_diff incr: dP_anls-dP_bckgr (m), {0} {1}, rdate={2}'.\
         format(xsct,pthbin[-13:]+fldiff,rdate))
  Hb, XX = psct.plot_rho_Jsct(xsct, dP_diff, ZZ, HH, LON, LAT,\
              fgnmb=nfg, stl=stl, sct_show=True, btx=btx, \
              rmin=-50., rmax=50., clrmp=clrmp)







