# Find HYCOM indeces of 4 points around CMEMS grid point
# for interpolation onto CMCES grid
#
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

f_new = False # True - start from 0, False - continue from last saved rec

pthssh = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/OBS/ssh_altimetry/'
flnm   = 'ssh_duacs_nrt_global_20230301_20230430.nc'
flssh  = pthssh + flnm

nc = ncFile(flssh,'r')
XA = nc.variables['longitude'][:].squeeze().data
YA = nc.variables['latitude'][:].squeeze().data
idma = XA.shape[0]
jdma = YA.shape[0]
ijdma = idma*jdma
# Only for ocean points:
SSH    = nc['adt'][0,:,:]

import mod_read_hycom as mhycom
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
XH, YH, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
idm  = XH.shape[1]
jdm  = XH.shape[0]
ijdm = idm*jdm
ymin = np.nanmin(YH) 

XH1d = XH.reshape(ijdm)
YH1d = YH.reshape(ijdm)

# Output save:
pthout = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/data_outp/'
floutp = pthout + 'duacs_nrt_global_025x025_hycom_indx.pkl'


import mod_misc1 as mmisc
IHYCOM = np.zeros((jdma,idma,4))-999
JHYCOM = np.zeros((jdma,idma,4))-999

if not f_new:
  print('Continue from saved output ' + floutp)

  with open(floutp,'rb') as fid:
    IHYCOM, JHYCOM = pickle.load(fid)
    

print('Finding closest HYCOM indices ')
icc  = 0
itot = 0
iskp = 0
rtimeS = time.time()
for iia in range(idma):
  for jja in range(jdma):
    x0 = XA[iia]
    y0 = YA[jja]
    itot += 1
    if y0 <= ymin:
      continue
    if np.isnan(SSH[jja,iia]):
      continue

    icc += 1
    if icc%1000 == 0:
      rtime0 = time.time()
      dtime  = (rtime0-rtimeS)/60.
      print('  {0:6.2f}% processed, skipped={2}, ellapsed {1:8.2f} min ...'.\
            format(itot/ijdma*100.,dtime,iskp))

    if not f_new:
      aa = IHYCOM[jja,iia,:]
      bb = JHYCOM[jja,iia,:]
      if (np.min(aa) > -1) and (np.min(bb) > -1):
        iskp += 1
        continue

    dst  = np.sqrt((x0-XH)**2+(y0-YH)**2)
    dmin = np.nanmin(dst)
# Check if simple way fails, use slower spherical dist:
    if dmin > 0.17:
      print('i={0} j={1}  min dist {2:10.4}dgr > 0.17, use spherical dist'.\
            format(iia,jja,dmin))
      dst = mmisc.dist_sphcrd(y0,x0,YH1d,XH1d)
      dst = dst.reshape(jdm,idm)
      dmin = np.nanmin(dst)
      print('  From dist_sphcrd: min dist={0:10.4f}m'.format(dmin))


    a1, a2 = np.where(dst == np.nanmin(dst))
    jj0 = a1[0]
    ii0 = a2[0]

    dhx = XA[iia]-XH[jj0,ii0]
    dhy = YA[jja]-YH[jj0,ii0]

    if abs(dhx) > 0.:
      ii1 = int(ii0+np.sign(dhx))
    else:
      ii1 = ii0-1

    if abs(dhy) > 0.:
      jj1 = int(jj0+np.sign(dhy))
    else:
      jj1 = jj0-1

    if ii1 < 0:
      ii1 = idm-1
    elif ii1 >= idm:
      ii1 = 0

    if jj1 >= jdm:
      jj1 = jdm-1

# Indices of 4 HYCOM grid points around cmems grid point:
# clockwise direction
    if ii1 > ii0 and jj1 > jj0:
      aa = np.array([ii0,ii0,ii1,ii1])
      bb = np.array([jj0,jj1,jj1,jj0])
    elif ii1 > ii0 and jj1 < jj0:
      aa = np.array([ii0,ii1,ii1,ii0])
      bb = np.array([jj0,jj0,jj1,jj1])
    elif ii1 < ii0 and jj1 < jj0:
      aa = np.array([ii0,ii0,ii1,ii1])
      bb = np.array([jj0,jj1,jj1,jj0])
    elif ii1 < ii0 and jj1 > jj0:
      aa = np.array([ii0,ii1,ii1,ii0])
      bb = np.array([jj0,jj0,jj1,jj1])
    else:
      aa = np.array([ii0,ii0,ii1,ii1])
      bb = np.array([jj0,jj1,jj1,jj0])
 
    IHYCOM[jja,iia,:] = aa
    JHYCOM[jja,iia,:] = bb

    if icc%1000 == 0:
      print('Saving  ' + floutp)
      with open(floutp,'wb') as fid:
        pickle.dump([IHYCOM,JHYCOM],fid)
      print('saved')

# Finished
print('Saving  ' + floutp)
with open(floutp,'wb') as fid:
  pickle.dump([IHYCOM,JHYCOM],fid)
print('All done ')


#
# Check
f_plt = False
if f_plt:
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

  plt.ion()
  fig1 = plt.figure(1,figsize=(7,7), constrained_layout=False)
  fig1.clf()
  ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
  ax1.plot(x0,y0,'r*')

  ax1.plot(xh,yh,'-',color=[0.6,0.6,0.6])
  ax1.plot(xh1,yh1,'.',color=[0,0.4,1])  # blue 1st pnt
  ax1.plot(xh2,yh2,'.',color=[0,1,0.2])  # green 2nd
  ax1.plot(xh3,yh3,'.',color=[0.8,0.4,0])
  ax1.plot(xh4,yh4,'.',color=[0.9,0.0,0.8])



    





