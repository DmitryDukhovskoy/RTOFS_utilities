# Validation: RTOFS ssh vs altimetry-derived SSH global gridded 
# CMEMS near-real time copernicus
# https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046/description
#
# This product is processed by the DUACS multimission altimeter 
# data processing system.
# DUACS  is the operationnal multimission production system of 
# altimeter data developed by CNES/CLS.
#
# downloaded from
# https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046/download?dataset=dataset-duacs-nrt-global-merged-allsat-phy-l4
# 

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text

import mod_utils as mutil
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_time as mtime

expt   = 'paraD5'
rdateS = '20230302'
rdateE = '20230430'
rdate0 = '20230314'
sfx    = 'n-24'

f_save = False  
pthsave = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/data_outp/ssh_intrp2cmems/'

import mod_read_hycom as mhycom
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
XH, YH, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
idm  = XH.shape[1]
jdm  = XH.shape[0]
ijdm = idm*jdm
ymin = np.nanmin(YH)

#
# Create land-ocean mask from SSH cmems
pthssh = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/OBS/ssh_altimetry/'
flnm   = 'ssh_duacs_nrt_global_20230301_20230430.nc'
flssh  = pthssh + flnm

nc    = ncFile(flssh,'r')
XA    = nc.variables['longitude'][:].squeeze().data
YA    = nc.variables['latitude'][:].squeeze().data
idma  = XA.shape[0]
jdma  = YA.shape[0]
ijdma = idma*jdma
SSHA  = nc['adt'][0,:,:]

Lmsk = np.where(np.isnan(SSHA),0,1)
Iall = np.where(Lmsk.flatten()>0)[0]
nall = Iall.shape[0]

# Output with HYCOM 4 indces for interpolation:
pthout = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/data_outp/'
floutp = pthout + 'duacs_nrt_global_025x025_hycom_indx.pkl'
with open(floutp,'rb') as fid:
  IHYCOM, JHYCOM = pickle.load(fid)

IHYCOM = IHYCOM.astype(np.int32)
JHYCOM = JHYCOM.astype(np.int32)

pthbin = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/' + expt +\
         '/rtofs.' + rdate0 + '/'
flarchv= 'rtofs_glo.t00z.' + sfx + '.archv' 
fina   = pthbin + flarchv + '.a'
finb   = pthbin + flarchv + '.b'

import mod_read_hycom as mhycom
importlib.reload(mhycom)
from mod_utils_fig import bottom_text
print('Processing '+fina)
IDM, JDM, KDM = mhycom.hycom_dim(fina,finb)
 
huge = 1.e20
rg   = 9806.
SSH,n1,m1,l1 = mhycom.read_hycom(fina,finb,'srfhgt')
SSH[SSH>huge] = np.nan
SSH = SSH/9.806

# --------------
#  Interpolate
# --------------
import mod_bilinear as mblnr

# Find basis functions for a reference rectangle:
phi1,phi2,phi3,phi4 = mblnr.basisFn_RectRef()


print('Starting interpolation HYCOM ssh --> 0.25 cmems grid')
SSHi = np.zeros((jdma,idma))*np.nan
kcc  = 0
tic  = time.time()
ticR = time.time() 
for ikk in range(nall):
  I1     = Iall[ikk]
  jj, ii = np.unravel_index(I1, Lmsk.shape) 

  kcc += 1
  if (kcc % 50000) == 0:
    toc  = time.time()
    print(' {0:5.2f}% done, {1:6.2f} min, lat={2:5.1f} ...'.\
          format(kcc/nall*100.,(toc-tic)/60.,YA[jj]))
#    print(' ')
    ticR = time.time()    

  iH = IHYCOM[jj,ii,:]
  jH = JHYCOM[jj,ii,:]
# 
# Should not happen but check:
  if np.min(iH < 0) or np.min(jH < 0):
#    print('Missing hycom indices i={0} j={1}'.format(ii,jj))
    continue
  
# Switch the order of vertices to match the reference square Sr
# Nodes should go c/clockwise starting lower-left node first

  eH = SSH[jH,iH]
  xH = XH[jH,iH]
  yH = YH[jH,iH]
  x0 = XA[ii]
  y0 = YA[jj]
# Check land in 4 grid points:
  In = np.where(np.isnan(eH))[0].shape[0]
  if In > 0 and In < 4:
    eH = np.where(np.isnan(eH), np.nanmean(eH), eH)
  elif In == 4:
# Land point in HYCOM
    continue

# Map X,Y ---> Xhat, Yhat on reference quadrialteral 
# i.e. map WOA grid coordinate to a reference quadrilateral 
# to do bilinear interpolation 
  xht, yht = mblnr.map_x2xhat(xH,yH,x0,y0)

# Perform interpolation on reference rectangle, that is 
# similar to interp on actual rectangle
  hintp = mblnr.bilin_interp(phi1,phi2,phi3,phi4,xht,yht,eH)

  SSHi[jj,ii] = hintp

if f_save:
  fsshout = 'sshintrp_' + expt + '_' + rdate0 + sfx + '.pkl'
  flsave  = pthsave + fsshout
  print('Saving  ' + flsave)
  with open(flsave,'wb') as fid:
    pickle.dump(SSHi,fid)
  

print('All done')



# Plot interpolated SSH
f_plotssh = True
if f_plotssh:
  plt.ion()
  fig1 = plt.figure(1,figsize=(8,7), constrained_layout=False)
  fig1.clf()

  ax1 = fig1.add_axes([0.1, 0.55, 0.8, 0.4]) 
  im1 = ax1.pcolormesh(SSH, vmin=-1, vmax=1)
  ax1.set_title(expt + ' RTOFS ssh ' + rdate0 + ' ' + flarchv)

  cax = fig1.add_axes([0.91, 0.55, 0.015, 0.4])
  clb = fig1.colorbar(im1, cax=cax, extend='both')
  cax.set_yticklabels(cax.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.tick_params(direction='in', length=5)

  ax2 = fig1.add_axes([0.1, 0.07, 0.8, 0.4]) 
  im2 = ax2.pcolormesh(SSHi, vmin=-1, vmax=1)
  ax2.set_title('RTOFS ssh interpolated 0.25 CMEMS ')

  cax2 = fig1.add_axes([0.91, 0.07, 0.015, 0.4])
  clb2 = fig1.colorbar(im2, cax=cax2, extend='both')
  cax2.set_yticklabels(cax2.get_yticks())
  ticklabs = clb2.ax.get_yticklabels()
  clb2.ax.set_yticklabels(ticklabs,fontsize=10)
  clb2.ax.tick_params(direction='in', length=5)

  btx = 'interp_hycom2cmems.py'
  bottom_text(btx, pos=[0.02,0.01])


# Check  mapping for individual grid box
f_plt = False
if f_plt:
  i0 = ii
  j0 = jj
  x0 = XA[i0]
  y0 = YA[j0]

  aa = IHYCOM[j0,i0,:]
  bb = JHYCOM[j0,i0,:]

  xh = XH[bb,aa]
  yh = YH[bb,aa]
  xh1 = xh[0]
  yh1 = yh[0]
  xh2 = xh[1]
  yh2 = yh[1]
  xh3 = xh[2]
  yh3 = yh[2]
  xh4 = xh[3]
  yh4 = yh[3]

# Coordinates in the reference square:
  xht0, yht0 = mblnr.map_x2xhat(xH,yH,x0,y0)
  xht1, yht1 = mblnr.map_x2xhat(xH,yH,xh1,yh1)
  xht2, yht2 = mblnr.map_x2xhat(xH,yH,xh2,yh2)
  xht3, yht3 = mblnr.map_x2xhat(xH,yH,xh3,yh3)
  xht4, yht4 = mblnr.map_x2xhat(xH,yH,xh4,yh4)

  plt.ion()
  fig2 = plt.figure(1,figsize=(7,7), constrained_layout=False)
  fig2.clf()
  ax1 = fig2.add_axes([0.1, 0.6, 0.35, 0.35])
  ax1.plot(x0,y0,'r*')
  ax1.plot(xh,yh,'-',color=[0.6,0.6,0.6])
  ax1.plot(xh1,yh1,'.',color=[0,0.4,1])  # blue 1st pnt
  ax1.plot(xh2,yh2,'.',color=[0,1,0.2])  # green 2nd
  ax1.plot(xh3,yh3,'.',color=[0.8,0.4,0])
  ax1.plot(xh4,yh4,'.',color=[0.9,0.0,0.8])
  ax1.set_title('Original grid')

  ax2 = fig2.add_axes([0.1, 0.1, 0.35, 0.35])
  ax2.plot(xht0,yht0,'r*')
  ax2.plot(xht1,yht1,'.',color=[0,0.4,1])  # blue 1st pnt
  ax2.plot(xht2,yht2,'.',color=[0,1,0.2])  # green 2nd
  ax2.plot(xht3,yht3,'.',color=[0.8,0.4,0])
  ax2.plot(xht4,yht4,'.',color=[0.9,0.0,0.8])
  ax2.set_title('Reference square')




