"""
  MOM6 utilities
  reading grid
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
from netCDF4 import Dataset as ncFile

def read_mom6grid(fgrid, grdpnt='hpnt'):
  """
    Read MOM6 grid
    default - coordinates of h-points (scalar) are derived
    options: 'qpnt' - coordinated in the upper-right corner
                      q-points, vorticity points

   https://mom6.readthedocs.io/en/main/api/generated/pages/Discrete_Grids.html
 
   q(i-1,j)   ------------  v(i,j) --------------  q(i,j)
     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  

   u(i-1,j)               * h(i,j)                u(i,j)

     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  
   q(i-1,j-1) ------------  v(i,j-1) ------------  q(i,j-1)


  """

  print('Reading MOM6 grid ' + grdpnt)
  print('Gird: ' + fgrid)
  nc  = ncFile(fgrid,'r')
  XX  = nc.variables['x'][:].data
  YY  = nc.variables['y'][:].data
  mm  = XX.shape[0]
  nn  = XX.shape[1]
  jdm = int((mm-1)/2)
  idm = int((nn-1)/2)


  if grdpnt == 'hpnt':
    LON = XX[1:mm:2, 1:nn:2]
    LAT = YY[1:mm:2, 1:nn:2]
  elif grdpnt == 'qpnt':
    LON = XX[0:mm-1:2, 1:nn-1:2]
    LAT = YY[0:mm-1:2, 1:nn-1:2]

  LON = np.where(LON < -180., LON+360., LON)

  return LON, LAT

def read_mom6depth(ftopo):
  """
    Read MOM6 depths at h-pnts
  """
  print('Reading MOM6 depths ' + ftopo)
  nc  = ncFile(ftopo,'r')
  HH  = nc.variables['depth'][:].data

# Convert depth to negatives:
# Land > 0
  HH  = np.where(HH > 0, -HH, 100.)

  return HH

def read_mom6lmask(ftopo):
  """
    Read MOM6 land mask
    1 = ocean, 0 = land
  """
  print('Reading MOM6 land mask ' + ftopo)
  nc   = ncFile(ftopo,'r')
  Lmsk = nc.variables['wet'][:].data

  return Lmsk


def read_mom6(foutp, fld, rLayer=-1, fnan=True, finfo=True):
  """
    Read 2D or 3D field from MOM6 output
    Specify Lyaer #  to read for 3D fields (from 1, ..., kdm)
    otherwise all layers are read in 3D array
    Replace out of range data with nan's - default

  """
  print('Reading MOM6 {0}: {1}'.format(fld, foutp))

  huge = 1.e18
  nc   = ncFile(foutp, 'r')

# Check variable dimension
# Assumed: AA(Time, dim1, dim2, optional: dim3)
  ndim = nc.variables[fld].ndim - 1 # Time excluded
  
# Check that time dim=1:
  tdim = nc.variables['Time'][:].data.shape[0]
  if tdim != 1:
    raise Exception('in netcdf output  Time dimensions is not 1: {0}'.\
                    format(tdim))

  if ndim == 2:
    FF = nc.variables[fld][:].squeeze().data
  elif ndim == 3:
    if rLayer <= 0:
      FF  = nc.variables[fld][:].squeeze().data
    else:
      klr = rLayer-1
      FF  = nc.variables[fld][:,klr,:,:].squeeze().data
  else:
    print('Output field dim should be 2 or 3')
    raise Exception('Check {0} dimension = {1}'.format(fld,ndim))

  AA = np.where(FF >= huge, np.nan, FF)
  if finfo:
    print('{0} min/max: {1} / {2}'.format(fld,np.nanmin(AA),np.nanmax(AA)))

  if fnan:
    FF = AA

  return FF

  

def create_time_array(date1, date2, dday, date_mat=False):
  """
    Create time array for plotting fields
    if date_mat: date1 and date2 are in mat format
    otherwise: dates = [year, month, day]

  """
  import mod_time as mtime
  importlib.reload(mtime)

  if not date_mat:
#    yr1, mo1, mday1 = date1 
#    yr2, mo2, mday2 = date2  
#
    dnmb1 = mtime.datenum(date1)
    dnmb2 = mtime.datenum(date2)
  else:
    dnmb1 = date1
    dnmb2 = date2

  ldays = np.arange(dnmb1,dnmb2,dday)   
  nrec  = ldays.shape[0]
  TPLT  = np.zeros(nrec, dtype=[('dnmb', float),
                                ('date', int, (4,)),
                                ('yrday', float)])

  for irc in range(nrec):
    dnmb = ldays[irc]
    DV   = mtime.datevec(dnmb)
    _, jday = mtime.dnmb2jday(dnmb)

    TPLT['dnmb'][irc]   = dnmb
    TPLT['yrday'][irc]  = jday
    TPLT['date'][irc,:] = DV[0:4]

  print('Creating TPLT Time array, start: ' + \
        '{0}/{1}/{2} {3}hr, '.format(TPLT['date'][0,0], \
        TPLT['date'][0,1], TPLT['date'][0,2], TPLT['date'][0,3]) + \
        'end: {0}/{1}/{2} {3}hr'.\
      format(TPLT['date'][-1,0], TPLT['date'][-1,1], TPLT['date'][-1,2],\
      TPLT['date'][-1,3]))

  return TPLT  

 




