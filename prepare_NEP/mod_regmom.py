"""
  Utility for regional MOM6-SIS2 
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pdb
import importlib
#import struct
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
from mod_utils_fig import bottom_text

def find_gridpnts_box(x0, y0, LON, LAT, dhstep=0.5):
  """
    Given pnt (x0,y0) find 4 grid points enclosing the pnt
    on a grid XX, YY - coordinates

     * -----------------  *
     |                    |
     |                    |
     |    *(x0,y0)        |
     |                    |
     |                    |
     |                    |
     * -----------------  *

    Specify offset for i,j indices if a subset XX, YY is ued
    instead of whole domain X and Y (to speed up distance calculation)
    to get the global indices 

    dhstep - approximate max grid horizontal stepping (in units of XX/YY)
             making it too small may result in errors as
             grid points outside this range will be discarded in order to
             accelerate searching algorithm
  """
  import mod_misc1 as mmisc1
  import mod_bilinear as mblnr

  mm = LON.shape[1]
  nn = LON.shape[0]
# Subsample the region:
  dy = dhstep
  JJ,II = np.where((LAT > y0-dy) & (LAT < y0+dy))

  XX = LON[JJ,II]
  YY = LAT[JJ,II]

  DD   = mmisc1.dist_sphcrd(y0,x0,YY,XX)
  kmin = np.where(DD == np.min(DD))[0][0]
#  jmin, imin = np.unravel_index(kmin, LON.shape)
#  jmin, imin = np.where(DD == np.min(DD))
  jmin = JJ[kmin]
  imin = II[kmin]
  xmin = LON[jmin, imin]
  ymin = LAT[jmin, imin]

# Find grid points around x0,y0:
  jm1   = jmin-1
  jp1   = jmin+1
  if jm1 < 0:
    jm1 = 0
  if jp1 >= mm:
    jp1 = mm-1

# Assuming global grid:
  im1   = imin-1
  ip1   = imin+1
  if im1 < 0:
    im1 = nn-1
  if ip1 >= nn:
    ip1 = 0

  Xs = LON[jm1:jp1+1,im1:ip1+1].flatten()
  Ys = LAT[jm1:jp1+1,im1:ip1+1].flatten()

# Find enclosed box:
  xrf = x0-1.
  if xrf < -180.:
    xrf = xrf+360.
  yrf = y0-1.
  xd0 = mmisc1.dist_sphcrd(y0,xrf,y0,x0)
  yd0 = mmisc1.dist_sphcrd(yrf,x0,y0,x0)
  ixx = []
  jxx = []
  for ibx in range(0,4):
    di = int(np.floor((ibx/2%2))*2-1)
    dj = ((ibx+1)%2)*2-1
    if abs(di) > 1 or abs(dj) > 1:
      raise Exception(f"Searching grid box di/dj>1: ibx={ibx} di={di} dj={dj}")
    inxt = imin+di
    jnxt = jmin+dj 

    IV = [imin, imin, inxt, inxt]
    JV = [jmin, jnxt, jnxt, jmin]
    IV, JV = mblnr.sort_gridcell_indx(IV, JV)

# Convert box vertices into Cartesian coord wrt x0,y0
    XX  = LON[JV,IV]
    YY  = LAT[JV,IV]
#    XV  = np.zeros((4))
#    YV  = np.zeros((4))
#    for ipp in range(0,4):
#      dlx     = mmisc1.dist_sphcrd(YY[ipp],xrf,YY[ipp],XX[ipp]) - xd0
#      dly     = mmisc1.dist_sphcrd(yrf,XX[ipp],YY[ipp],XX[ipp]) - yd0
#      XV[ipp] = dlx
#      YV[ipp] = dly
#
    XV, YV, x0c, y0c = mblnr.lonlat2xy_wrtX0(XX, YY, x0, y0)
    INp     = mmisc1.inpolygon_1pnt(x0c, y0c, XV, YV)
    if INp:
      ixx = IV
      jxx = JV
      break

# Plot box:
# ax1.plot(XV,YV,'.-')    
# ax1.plot(0,0,'r*')    # pnt should be inside XV,YV
  if len(ixx) == 0:
    print(f"Searching for {x0},{y0} enclosing grid box")
    raise Exception("Failed no grid box for x0,y0 found")
 
  ixx = np.array([imin, imin, inxt, inxt]).astype(int)
  jxx = np.array([jmin, jnxt, jnxt, jmin]).astype(int)
  
  return ixx, jxx

def fill_land(aa1,aa2,aa3,aa4,HH,A3d,JJ,II,Jocn,Iocn):
  """
    Fill all nans in 1D arrays - land points
    at least 1 array should be ocean pnt
  """
  la1  = len(np.where(np.isnan(aa1))[0])
  la2  = len(np.where(np.isnan(aa2))[0])
  la3  = len(np.where(np.isnan(aa3))[0])
  la4  = len(np.where(np.isnan(aa4))[0])
  kdmh = len(aa1)

  LL  = np.array([la1,la2,la3,la4])
  iL  = np.where(LL == kdmh)[0]
  inL = np.where(LL < kdmh)[0] 

# This should not happen - no land points
# assumed at least 1 array is nans
  if len(iL) == 0:
    return aa1, aa2, aa3, aa4

  AP  = np.column_stack((aa1,aa2,aa3,aa4))
# All land - find closes ocean point
  if len(iL) == 4:
#    raise Exception("All 4 1D arrays - land points, check Land masks")
    for ild in range(4):
      il1 = II[ild]
      jl1 = JJ[ild]    
      D   = np.square((Iocn-float(il1))**2 + (Jocn-float(jl1))**2)
      ixx = np.argmin(D)
      ioc = Iocn[ixx]
      joc = Jocn[ixx]
      AP[:, ild] = A3d[:,joc,ioc]

  APm = np.nanmean(AP, axis=1) 
  AP[:,iL] = APm[:,None]
 
  return AP[:,0].squeeze(), AP[:,1].squeeze(),\
         AP[:,2].squeeze(), AP[:,3].squeeze() 

def check_bottom(AA):
  """
    Make sure there are no NaNs in the 1D vertical data
    AA is 1D array
    If HYCOM land mask does not match MOM's --> all nans in the 1D profile
    from HYCOM
  """
  if AA[0] == np.nan:
    print(f"1D profile: all values are nans, Land mask mismatch")
    return 

  izb = np.argwhere(AA == np.nan)
  if len(izb) == 0:
    return

  AA[izb] = AA[min(izb)-1]

  return AA
  



