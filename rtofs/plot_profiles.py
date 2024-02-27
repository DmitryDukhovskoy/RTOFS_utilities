"""
  Plot T/S/U profiles  
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
#importlib.reload(rncoda)

# incup files = n-24 files
# bkgrd files = n-00 restart files (from previous f/cast)
date_type = 'incup'  # incup - incr update, bkgrd, anls, NRL_rstrt
fldplt    = 'temp'  # temp, salin, u-vel., v-vel. restart: saln
hr0       = 0

rdate0  = '20220616'  # analysis date - when nowcast is done
rdate0  = '20230125'
rdate   = '20230117'
#pth1    = '/scratch2/NCEPDEV/marine/Dan.Iredell/para8d/rtofs.'
#pth1    = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/rtofs_expts/'+\
#           'rtofs_para8.1b/rtofs.'
pth1    = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/'+\
           'rtofs_para7b/hycom/rtofs.'
pth1    = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/wcoss2.paraA/hycom/rtofs.'

date_outp = rncoda.date_type_output(rdate0,"incup")
date_fcst = rncoda.date_type_output(rdate0,"fcast")
date_dir  = date_fcst[0:8]

if pth1[-15:] == 'ncoda_archv_inc':
  pthbin = pth1+'/'
else:
  pthbin  = pth1+rdate+'/'

# Specify location:
#xloc = 175.77
#yloc = -1.34
xloc = 149.6
yloc = 29.0
iloc = []
jloc = []


pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

jday = rncoda.rdate2julian(rdate)
yr, mo, mday, hr = rncoda.parse_rdate(rdate)
hr   = hr0

# Background field:
#flnm = 'archv.{0}_{1:03d}_00'.format(yr,jday)   # background field
# Increment:
#flnm = 'archv_1_inc.2022_{0:03d}_00'.format(jday)  # added increments  ncoda_archv_inc
# incrementally updated fields:
flnm = 'archv.{0}_{1:03d}_{2:02d}'.format(yr,jday,hr)
flnm = 'rtofs_glo.t00z.n-24.archv'  # incup file
flnm = 'rtofs_glo.t00z.n-24.restart'  # NRL restart

fina = pthbin+flnm+'.a'
finb = pthbin+flnm+'.b'

IDM = 4500
JDM = 3298
KDM = 41

f_restart = flnm[0:7]=='restart'
if not f_restart:
  f_restart = flnm[-7:]=='restart'

get_topo = True
ftopo = 'regional.depth'
fgrid = 'regional.grid'


import mod_read_hycom
#importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo, read_hycom, \
                           read_topo, read_hycom_restart
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text

print('Processing '+fina)

def print_1D(A,wd=8,prc=2):
  ndim1 = A.shape[0]
  for k in range (ndim1):
    print('{0}: {1:{width}.{precis}f}'.format(k+1,A[k],width=wd,precis=prc))


huge = 1.e20
rg   = 9806.
fld  = 'thknss'
if f_restart:
  fld = 'dp'
  F,nn,mm,ll = read_hycom_restart(fina,finb,fld,IDM,JDM,rLayer=1)
else:
  F,nn,mm,ll = read_hycom(fina,finb,fld,rLayer=1)

F[np.where(F>huge)] = np.nan
F = F/rg
F[np.where(F<0.001)] = 0.

# For example 30S, 30W is  i=3199 j=1112
dH = np.zeros((ll,mm,nn))
dH[0,:,:] = F
for kk in range(2,ll+1):
  if f_restart:
    fld = 'dp'
    F,nn,mm,llm = read_hycom_restart(fina,finb,fld,IDM,JDM,rLayer=kk)
  else:
    F,nn,mm,lmm = read_hycom(fina,finb,fld,rLayer=kk)

  F = F/rg
  F[np.where(F>huge)] = np.nan
  F[np.where(F<0.001)] = 0.
  dH[kk-1,:,:] = F

# ================================
ZZ, ZM = zz_zm_fromDP(dH)
kdm = ZM.shape[0]
jdm = ZM.shape[1]
idm = ZM.shape[2]

if get_topo:
#  HH = read_topo(pthgrid,ftopo,nn,mm)
  LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)
  get_topo = False  


dZZ = np.abs(np.diff(ZZ, axis=0))
btx = 'plot_profiles.py'


fld = fldplt
print('Reading '+fina)

fld = fldplt
if f_restart and fldplt == 'salin':
  fld = 'saln'
  fldplt = 'saln'

A3d = np.array([])
for kk in range (1,KDM+1):
  if f_restart:
    F,nn,mm,llm = read_hycom_restart(fina,finb,fld,IDM,JDM,rLayer=kk)
  else:
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

#import mod_utils as mutil
#importlib.reload(mutil)

# ==========================
# Plot selected W-E profiles  
import plot_sect as psct
importlib.reload(psct)
i0, j0, imm, jmm = psct.find_indx_lonlat(xloc,yloc,LON,LAT)
lon0 = LON[j0,i0]
lat0 = LAT[j0,i0]

#  drfnm = "rtofs."+rdate+"/"+flnm
drfnm = pthbin[-16:]+flnm
stl = ('{0}  {9} \n{1} \n{2}/{3}/{4} \n i={5} j={6}, lon={7:5.2f}, lat={8:5.2f}'.\
        format(fldplt,fina[-35:],rdate[0:4],rdate[4:6],rdate[6:8],\
        i0,j0,lon0,lat0,date_type))
#psct.plot_1prof_map(A3d[:,j0,i0],ZM[:,j0,i0],HH,ctl=stl,zlim=-800.,i1=i0,j1=j0)
psct.plot_Nprof_map(A3d,ZM,i0,j0,HH,dI=0,ctl=stl,zlim=-800.)
bottom_text(btx)

#  print_1D(dZZs[:,jp0])
#psct.print_2D(ZM[:,j0,i0],A3d[:,j0,i0],kend=kdm)







