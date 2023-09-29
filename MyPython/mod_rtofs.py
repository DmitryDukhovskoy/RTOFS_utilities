"""
  RTOFS utilities
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct


def find_indx_lonlat(x0, y0, X0, Y0):
  """
  Map lon/ lat ---> index (i,j)
  Find closest grid point (ii0,jj0) to lon/lat coordinate
  For W-E sections, provide xsct name 
  then indices are given relative to the section 1st index
  Output: ii0, jj0 - indices closest to x0,y0 (lon, lat)
          ip0, jp0 - indices relative to 1st pnt in the section
  """
  if x0 > 180.:
    x0 = x0-360.

  XX = X0.copy()
  YY = Y0.copy()

  dmm = np.sqrt((XX-x0)**2+(YY-y0)**2)
  jj0, ii0 = np.where(dmm == np.min(dmm)) # global indices
  jj0 = jj0[0]
  ii0 = ii0[0]

#
# Check for singularity along 180/-180 longitude
  IDM = XX.shape[1]
  if ii0 > 0 and ii0 < IDM:
    xm1 = XX[jj0,ii0-1]
    xp1 = XX[jj0,ii0+1]
    if abs(xm1-x0) > 180. or abs(xp1-x0)>180.:
      J, I = np.where(XX < 0.)
      XX[J,I] = XX[J,I]+360.

      if x0 < 0:
        x0 = x0+360.

    dmm = np.sqrt((XX-x0)**2+(YY-y0)**2)
    jj0, ii0 = np.where(dmm == np.min(dmm)) # global indices
#
  return ii0[0], jj0[0]


def interp_lonlat2indx(x0,y0,LON,LAT):
  """
    Map lon/ lat ---> index (i,j)
    Find "exact" index (float number) by interpolation, i.e.
    Map x0,y0 (lon,lat) into index space
    LON, LAT - 2D arrays of long and latitudes of RTOFS grid
  """
  import mod_bilinear as mblinr
#
# Find 4 nodes around x0,y0
# arranging counter-clockwise 
  i1, j1 = find_indx_lonlat(x0,y0,LON,LAT)
  x1 = LON[j1,i1]
  y1 = LAT[j1,i1]

  JDM = LON.shape[0]
  IDM = LON.shape[1]

  if x1 <= x0:
    iv1 = i1
    iv2 = i1+1
  else:
    iv1 = i1-1
    iv2 = i1

  if y1 <= y0:
    jv1 = j1
    jv2 = j1+1
  else:
    jv1 = j1-1
    jv2 = j1

  iht1 = iv1  # for interpolation keep original values
  iht2 = iv2
  jht1 = jv1
  jht2 = jv2
  if iv1 < 0:
    iv1 = IDM-1
  if jv1 < 0:     # impossible case
    jv1 = 0
  if iv2 > IDM-1:
    iv2 = 0
  if jv2 > JDM-1:
    jv2 = JDM-1
#
# Form 1D arrays of X and Y coordinates of vertices
# of a grid cell including the interp point
  XX = np.array([LON[jv1,iv1],LON[jv1,iv2],LON[jv2,iv2],LON[jv2,iv1]])
  YY = np.array([LAT[jv1,iv1],LAT[jv1,iv2],LAT[jv2,iv2],LAT[jv2,iv1]])

# Singularity:
  if (np.max(XX)-np.min(XX)) > 180.:
    II =  np.where(XX<0)
    XX[II] = XX[II]+360.
    if x0 < 0.:
      x0 = x0+360.

# Vector of indices:
  IX = np.array([iht1,iht2,iht2,iht1])
  JX = np.array([jht1,jht1,jht2,jht2])

# Map X,Y ---> Xhat, Yhat on reference quadrialteral 
# i.e. map WOA grid coordinate to a reference quadrilateral 
# to do bilinear interpolation 
  xht, yht = mblinr.map_x2xhat(XX,YY,x0,y0)

# Perform interpolation on reference rectangle, that is 
# similar to interp on actual rectangle
  Iint = mblinr.blinrint(xht,yht,IX)
  Jint = mblinr.blinrint(xht,yht,JX)

# Check negative values:
  if Iint < 0:
    Iint = IDM+Iint
  elif Iint > IDM-1:
    Iint = Iint-(IDM-1)

  return Iint, Jint


def interp_indx2lonlat(ii0,jj0,LON,LAT):
  """
    Map index (i,j) ---> lon/lat
    Given "exact" index (float number) by interpolation, i.e.
    LON, LAT - 2D arrays of long and latitudes of RTOFS grid
  """
  import mod_bilinear as mblinr
  JDM = LON.shape[0]
  IDM = LON.shape[1]

# Find 4 nodes around x0,y0
# arranging counter-clockwise 
  iv1 = int(np.floor(ii0))
  jv1 = int(np.floor(jj0))
  iv2 = iv1 + 1
  jv2 = jv1 + 1
 
  iht1 = iv1  # for interpolation keep original values
  iht2 = iv2
  jht1 = jv1
  jht2 = jv2
  if iv1 < 0:
    iv1 = IDM-1
  if jv1 < 0:     # impossible case
    jv1 = 0
  if iv2 > IDM-1:
    iv2 = 0
  if jv2 > JDM-1:
    jv2 = JDM-1
#
# Form 1D arrays of X and Y coordinates of vertices
# of a grid cell including the interp point
  XX = np.array([LON[jv1,iv1],LON[jv1,iv2],LON[jv2,iv2],LON[jv2,iv1]])
  YY = np.array([LAT[jv1,iv1],LAT[jv1,iv2],LAT[jv2,iv2],LAT[jv2,iv1]])

# Singularity:
  if (np.max(XX)-np.min(XX)) > 180.:
    II =  np.where(XX<0)
    XX[II] = XX[II]+360.
    if x0 < 0.:
      x0 = x0+360.

# Vector of indices:
  IX = np.array([iht1,iht2,iht2,iht1])
  JX = np.array([jht1,jht1,jht2,jht2])

# Map X,Y ---> Xhat, Yhat on reference quadrialteral 
# i.e. map WOA grid coordinate to a reference quadrilateral 
# to do bilinear interpolation 
  xht, yht = mblinr.map_x2xhat(IX,JX,ii0,jj0)

# Perform interpolation on reference quadrilateral, that is 
# similar to interp on actual quadrilateral
  lon0 = mblinr.blinrint(xht,yht,XX)
  lat0 = mblinr.blinrint(xht,yht,YY)

  return lon0, lat0


def find_max3d(A3d,dZZ):
  """
  Find max values over all layers for 3D array
  avoid zero thickness layers
  """
  kdm = A3d.shape[0]
  jdm = A3d.shape[1]
  idm = A3d.shape[2]

  maxA = np.zeros((jdm,idm))-2.e20
  for k in range(kdm):
    aa    = np.squeeze(A3d[k,:,:])
    dz    = np.squeeze(dZZ[k,:,:])
    dmm   = aa-maxA
    dmm   = np.where(dz<1.e-6,-1.e3,dmm)
    maxA  = np.where(dmm>0.,aa,maxA)
    print('Max value: k={0} max={1:10.4f}'.format(k,np.nanmax(aa)))

  return maxA

def find_min3d(A3d,dZZ):
  """
  Find min values over all layers for 3D array
  avoid zero thickness layers
  """
  kdm = A3d.shape[0]
  jdm = A3d.shape[1]
  idm = A3d.shape[2]

  minA = np.zeros((jdm,idm))+2.e20
  for k in range(kdm):
    aa    = np.squeeze(A3d[k,:,:])
    dz    = np.squeeze(dZZ[k,:,:])
    dmm   = minA-aa
    dmm   = np.where(dz<1.e-6,-1.e3,dmm)  # zero-thickness layers
    minA  = np.where(dmm>0.,aa,minA)
    print('Min value: k={0} min={1:10.4f}'.format(k,np.nanmin(aa)))

  return minA

def find_max3d(A3d,dZZ):
  """
  Find max values over all layers for 3D array
  avoid zero thickness layers
  dZZ is 3D layer thickness array
  """
  kdm = A3d.shape[0]
  jdm = A3d.shape[1]
  idm = A3d.shape[2]

  maxA = np.zeros((jdm,idm))-2.e20
  for k in range(kdm):
    aa    = np.squeeze(A3d[k,:,:])
    dz    = np.squeeze(dZZ[k,:,:])
    dmm   = aa-maxA
    dmm   = np.where(dz<1.e-6,-1.e3,dmm)
    maxA  = np.where(dmm>0.,aa,maxA)
    print('Max value: k={0} max={1:10.4f}'.format(k,np.nanmax(aa)))

  return maxA

def find_min3d(A3d,dZZ):
  """
  Find min values over all layers for 3D array A3d
  avoid zero thickness layers 
  dZZ is 3D layer thickness array
  """
  kdm = A3d.shape[0]
  jdm = A3d.shape[1]
  idm = A3d.shape[2]

  minA = np.zeros((jdm,idm))+2.e20
  for k in range(kdm):
    aa    = np.squeeze(A3d[k,:,:])
    dz    = np.squeeze(dZZ[k,:,:])
    dmm   = minA-aa
    dmm   = np.where(dz<1.e-6,-1.e3,dmm)  # zero-thickness layers
    minA  = np.where(dmm>0.,aa,minA)
    print('Min value: k={0} min={1:10.4f}'.format(k,np.nanmin(aa)))

  return minA







