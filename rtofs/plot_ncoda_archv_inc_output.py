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
import importlib
import struct
#import netCDF4
#from netCDF4 import Dataset as ncFile
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

# Analysis date, i.e. rtofs f/cast date - 1 day
expt    = 'paraD'
#rdate   = '20220616'
rdate   = '20230418'
pthbs   = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/'
pth1    = pthbs+'hycom/ncoda_archv_inc/'
pthbin  = pth1
pthbgr  = pthbs+'hycom/ncoda_archv_inc/'  # background field dir
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

fldplt = 'salin'  # temp salin u-vel. v-vel.

import mod_time as mtime
#jday = rncoda.rdate2julian(rdate)
#yr, mo, mday, hr = rncoda.parse_rdate(rdate)
jday             = mtime.rdate2jday(rdate)
yr, mo, mday, hr = mtime.parse_rdate(rdate)

# Select - updated (True) or background (False) fields:
f_updated = True
# flupdt - updated fields=background+increments output from ncoda_archv_inc
# flbgr  - background field from HYCOM f/cast
flupdt = 'archv_1_inc.{0}_{1:03d}_00'.format(yr,jday)
flbgr   = 'archv.{0}_{1:03d}_00'.format(yr,jday) 


get_topo = True
ftopo = 'regional.depth'
fgrid = 'regional.grid'

import mod_read_hycom as rdhycom
#importlib.reload(rdhycom)
import mod_utils as utls
#importlib.reload(utls)
from mod_utils_fig import bottom_text

btx = 'plot_ncoda_archv_inc_output.py'

#  HH = read_topo(pthgrid,ftopo,nn,mm)
LON, LAT, HH = rdhycom.read_grid_topo(pthgrid,ftopo,fgrid)

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

IDM, JDM, KDM = rdhycom.hycom_dim(fina,finb)

print('Reading layer thickness: '+fina)
fld  = 'thknss'
F,nn,mm,ll = rdhycom.read_hycom(fina,finb,fld,rLayer=1)
F[np.where(F>huge)] = np.nan
F = F/rg
F[np.where(F<0.001)] = 0.

dH = np.zeros((ll,mm,nn))
dH[0,:,:] = F
for kk in range(2,ll+1):
  F,nn,mm,lmm = rdhycom.read_hycom(fina,finb,fld,rLayer=kk)
  F = F/rg
  F[np.where(F>huge)] = np.nan
  F[np.where(F<0.001)] = 0.
  dH[kk-1,:,:] = F

ZZ, ZM = rdhycom.zz_zm_fromDP(dH, f_btm=False)
kdm = ZM.shape[0]
jdm = ZM.shape[1]
idm = ZM.shape[2]

print('Processing '+fina)

fld = fldplt
A3d = np.array([])
for kk in range (1,KDM+1):
  F,n1,m1,l1 = rdhycom.read_hycom(fina,finb,fld,rLayer=kk)
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
# Plot vertical sections
# Plot vertical section of background and increm rho
# and difference btw the two
#SCTP = ["GoM1"]
#SCTP = ["GoM1","SeaJpn"]
SCTP = ["NPac5"]  # Sulu Sea
nsct = len(SCTP)
import mod_utils as mutil
importlib.reload(mutil)
XSCT = mutil.rtofs_sections()
SCTnames = list(XSCT.keys())

plt.ion()
import plot_sect as psct
importlib.reload(psct)

clrmp = copy(plt.cm.BrBG_r)
if fldplt == "temp":
#  clrmp = copy(plt.cm.rainbow)
  cmpr = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  cmpr.set_bad(color=[1,1,1])
  rmin = -2.
  rmax = 28.
  fld1 = 'SctT'
elif fldplt == "salin":
#  clrmp = copy(plt.cm.gist_ncar_r)
  cmpr = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  cmpr.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 30.
  rmax = 40.
  fld1 = 'SctS'
elif fldplt == 'u-vel.' or fldplt == 'v-vel.':
#  rmin, rmax = psct.minmax_clrmap(dmm,pmin=5.,pmax=95.,cpnt=0.01,fsym=True)
  rmin = -0.1
  rmax = 0.1
  clrmp = copy(plt.cm.seismic)

clrmp.set_bad(color=[0.3,0.3,0.3])

jlrs = -1
nfg = 0
f_intract = True
for ii in range(nsct):
  nfg += 2

  sct_name = SCTP[ii]
  isc = SCTnames.index(sct_name)

# Specify location of the section:
# xsection in mod_utils <--- get di, dj, ix0, jx0 from XSCT
  ix1  = XSCT[sct_name]["x1"]
  ix2  = XSCT[sct_name]["x2"]
  jx0  = XSCT[sct_name]["y0"]
  jx1  = XSCT[sct_name]["y1"]
  jx2  = XSCT[sct_name]["y2"]
  ix0  = XSCT[sct_name]["x0"]
  Regn = XSCT[sct_name]["Reg"]
  Sname= XSCT[sct_name]["Name"]

  lat0 = LAT[jx0,ix0]
  lon0 = LON[jx0,ix0]
  lon1 = LON[jx0,ix1]
  lon2 = LON[jx0,ix2]
  lat1 = LAT[jx1,ix0]
  lat2 = LAT[jx2,ix0]
#
  if f_updated:
  # Updated fields with increments:
    fhcm = flupdt+'.a'
    cfls = 'Updated fld'
  else:
  # Background fields:
    fhcm = flbgr+'.a'
    cfls = 'B/ground fld'

  ss1 = '{0} {1} {3} date: {2} \n'.format(expt,fhcm,rdate,cfls)
  ss1 = ss1 + 'Fort i/i, j/j: {0}/{1}, {2}/{3}\n'.format(ix1+1,ix2+1,jx0+1,jx0+1)
  ss1 = ss1 + 'Lon, lat: {0:7.2f}/{1:7.2f}E, {2:7.2f}/{3:7.2f}N\n'.\
         format(lon1,lon2,lat0,lat0)
  ss1 = ss1 + 'Fort i0={0}, j0={1}, lon0={2:7.2f}E, lat0={3:7.2f}N'.\
         format(ix0+1,jx0+1, lon0, lat0)

  stl  = 'EW section {2}, {0} {1}'.format(Sname, sct_name, fld)
  print('Plotting '+stl)
  ZZs, Hb, XX, dZZs = psct.plot_EWsectTS(ix1, ix2, ix0, jx0, ZZ, A3d, HH, LON, LAT,\
                      cmpr, lrstart=42, stl=stl, fgnmb=nfg, btx=btx, \
                      sinfo=ss1, f_intract=f_intract)


  # S- N section:
  ss2 = '{0} {1} {3} date: {2} \n'.format(expt,fhcm,rdate,cfls)
  ss2 = ss2 + 'Fort i/i, j/j: {0}/{1}, {2}/{3}\n'.format(ix0+1,ix0+1,jx1+1,jx2+1)
  ss2 = ss2 + 'Lon, lat: {0:7.2f}/{1:7.2f}E, {2:7.2f}/{3:7.2f}N\n'.\
         format(lon0,lon0,lat1,lat2)
  ss2 = ss2 + 'Fort i0={0}, j0={1}, lon0={2:7.2f}E, lat0={3:7.2f}N'.\
         format(ix0+1,jx0+1, lon0, lat0)

  stl  = 'SN section {2}, {0} {1}'.format(Sname, sct_name, fld)
  print('Plotting '+stl)
  ZZn,Hbn,XXn,dZZn = psct.plot_SNsectTS(jx1, jx2, ix0, jx0, ZZ, A3d, HH, LON, LAT,\
                      cmpr, lrstart=42, stl=stl, fgnmb=nfg+1, btx=btx, \
                      sinfo=ss2, f_intract=f_intract)

  f_3d = False
  if f_3d:
    imm = 600
    jmm = 1606
#    imm = 1372
#    jmm = 1336
#    imm = 2954
#    jmm = 782
    dz  = dH[:,jmm,imm]
    aa  = A3d[:,jmm,imm]
    psct.print_3D(dz,aa,np.cumsum(dz),wd=14,prc=6)



