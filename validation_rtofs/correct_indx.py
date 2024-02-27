# Fixing bug in find_hycom2cmems_indx.py
# index at the boundary
#
# Find HYCOM indeces of 4 points around CMEMS grid point
# for interpolation onto CMCES grid
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
with open(floutp,'rb') as fid:
  IHYCOM, JHYCOM = pickle.load(fid)
    

print('Finding closest HYCOM indices ')
icc  = 0
iskp = 0
rtimeS = time.time()
for iia in range(idma):
  for jja in range(jdma):
    icc += 1
    if icc%1000 == 0:
      rtime0 = time.time()
      dtime  = (rtime0-rtimeS)/60.
      print('  {0:6.2f}% processed, skipped={2}, ellapsed {1:8.2f} min ...'.\
            format(icc/ijdma*100.,dtime,iskp))

    aa = IHYCOM[jja,iia,:]
    bb = JHYCOM[jja,iia,:]
    if (np.min(aa) < 0) and (np.min(bb) < 0):
      iskp += 1
      continue

    if np.max(aa) >= idm:
      aa = np.where(aa>=idm, idm-1, aa)

    if np.max(bb) >=jdm:
      bb = np.where(bb>=jdm, jdm-1, bb)

    IHYCOM[jja,iia,:] = aa
    JHYCOM[jja,iia,:] = bb
    


print('Saving  ' + floutp)
with open(floutp,'wb') as fid:
  pickle.dump([IHYCOM,JHYCOM],fid)



    





