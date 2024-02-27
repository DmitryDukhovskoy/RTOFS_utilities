"""
  Plot section with vertical layers
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
import importlib
import struct
import pickle
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

rdate0  = '20230407'  # fcast date used in rtofs.YYYMMDD
# rdate = date of the f/cast
# incup 
#toutp = "incup"   # type of archv: incup - from 6hr incr update, fcast
hr0   = -24   # use -24 - for ncup, n-24 ... n00, hindcast, f01, ..., f12, .. f/cast
             # all hours are wrt to f/cast date rdate0

pth1    = '/scratch2/NCEPDEV/marine/Dan.Iredell/wcoss.paraB/rtofs.'
#pth1    = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/rtofs_expts/'+\
#           'rtofs_para8.1b/rtofs.'
#pth1    = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/'+\
#           'rtofs_para7b/hycom/rtofs.'
#pth1    = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/wcoss2.paraX/rtofs.'

#date_outp = rncoda.date_type_output(rdate0,toutp)
#date_fcst = rncoda.date_type_output(rdate0,"fcast")
#date_dir  = date_fcst[0:8]


if pth1[-15:] == 'ncoda_archv_inc':
  pthbin = pth1+'/'
else:
  pthbin  = pth1+rdate0+'/'

rdate = rdate0


pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

#flnm = 'rtofs_glo.t00z.n-24.archv'    # 6hr ingested NCODA incremental updates
#flnm = 'restart_r2022111600_930'      # NRL restart with no NCODA 
#flnm = 'rtofs_glo.t00z.n00.archv'    # 24-hr hindcast or background for next day
#if rdate == '20220616':
#  flnm = 'archv.2022061600'
jday = rncoda.rdate2julian(rdate)
yr, mo, mday, hr = rncoda.parse_rdate(rdate)
hr   = hr0

# Background field:
#flnm = 'archv.{0}_{1:03d}_00'.format(yr,jday)   # background field
# Increment:
#flnm = 'archv_1_inc.2022_{0:03d}_00'.format(jday)  # added increments  ncoda_archv_inc
# incrementally updated fields:
#flnm = 'archv.{0}_{1:03d}_{2:02d}'.format(yr,jday,hr) # not renamed
if hr0 < 0:
  sfx = 'n-{0:02d}'.format(abs(hr0))
elif hr0 > 0 and hr0 <= 24:
  sfx = 'n{0:02d}'.format(hr0)
else:
  sfx = 'f{0:02d}'.format(hr0)


flnm = 'rtofs_glo.t00z.{0}.archv'.format(sfx)   # renamed

fina = pthbin+flnm+'.a'
finb = pthbin+flnm+'.b'

IDM = 4500
JDM = 3298

f_restart = flnm[0:7]=='restart'

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

# Read 1 layer to obtain dimensions:
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
def find_indx_lonlat(x0,y0):
  """
  Find closest grid point to lon/lat coordinate
  """
  if x0 > 180.:
    x0 = x0-360.
  
  dmm = np.sqrt((LON-x0)**2+(LAT-y0)**2)
  jj0, ii0 = np.where(dmm == np.min(dmm))

  return ii0[0], jj0[0]

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
import mod_utils as mutil
importlib.reload(mutil)

grdZ = mutil.anls_lrthkn(dZZ,lrintrf=False)
btx = 'plot_xsct_layres.py'


# ==========================
# Plot selected W-E profiles  
import plot_sect as psct
importlib.reload(psct)

itst = 2412
jtst = 1798
f_plt = True
if f_plt:
#  xsct = "SeaJpn"
#  xsct = "Test1"
#  xsct = "Test2"
  xsct = "GoM1"
#  xsct = "NAtl1"
# find index to print out layers:
  ipp,jpp,IP0,JP0 = psct.find_indx_lonlat(-37.681,-47.96,LON,LAT,xsct=xsct)

  if xsct == "GoM1":
    jp0 = 72  # GoM corresponds the itest, jtest in ncoda_archv_inc
  elif xsct == "SeaJpn":
    jp0 = 56
  elif xsct == "Test1":
    jp0 = 102
  elif xsct == "NAtl1":
    jp0 = 50
  elif xsct == "Test2":
    jp0 = 60

#  drfnm = "rtofs."+rdate+"/"+flnm
  drfnm = pthbin[-16:]+flnm
  ZZs, Hb, XX, dZZs = psct.plot_Jsection(xsct,ZZ,HH,LON,LAT,drfnm,sct_show=True,\
                      fgnmb=5,dcntr=1, zstart=-1000., jlrs=jp0, btx=btx)
#  print_1D(dZZs[:,jp0])
  psct.print_2D(ZZs[:,jp0],dZZs[:,jp0],kend=kdm)
 






