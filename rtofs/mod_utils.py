"""
  Utilities 
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
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

sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython')
import mod_read_hycom
#importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo, read_hycom, read_topo
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text

def ncoda_depths(zneg=True):
  """
    NCODA interface and mid-point depths
    zi = from hycom/ncoda_archv_inc/zi.txt
  """
  ZI = np.array([      
      0.00,
      1.00,
      3.00,
      5.00,
     11.00,
     21.00,
     35.00,
     53.00,
     67.00,
     85.00,
     99.00,
    117.00,
    131.00,
    149.00,
    171.00,
    189.00,
    211.00,
    229.00,
    251.00,
    269.00,
    291.00,
    309.00,
    371.00,
    389.00,
    451.00,
    469.00,
    531.00,
    589.00,
    651.00,
    709.00,
    771.00,
    829.00,
    971.00,
   1029.00,
   1271.00,
   1329.00,
   1571.00,
   1629.00,
   1871.00,
   2129.00,
   2371.00,
   2629.00])

  kzi = ZI.shape[0]
#
# Follow ncoda_archv_inc.f to define the mid-depths points in NCODA
# Note that in the ncoda code zz is NCODA z-level mid-depths points
# and zi is array with z-level interface depths
  ZM = np.zeros((kzi-1))
  for kk in range(1,kzi-1):
    ZM[kk] = 0.5*(ZI[kk-1]+ZI[kk])

  kzm = ZM.shape[0]

  if zneg:
    ZI = -ZI
    ZM = - ZM

  return ZI, kzi, ZM, kzm

def find_indx_lonlat(x0,y0,X0,Y0,xsct="none"):
  """
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
# If section provided, indices for this section:
  ip0 = -1
  jp0 = -1

  if not xsct=="none":
    SCT = Jsections()
    aa  = SCT.get(xsct)
    i1  = aa[0]
    i2  = aa[1]
    j1  = aa[2]
    ip0 = ii0[0]-i1
    jp0 = jj0[0]-j1

    return ii0[0], jj0[0], ip0, jp0
  else:
    return ii0[0], jj0[0]

def interp_indx_lonlat(x0,y0,LON,LAT):
  """
    Map x0,y0 (lon,lat) into index space, i.e.
    find "exact" index (float number) by interpolation
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


def anls_lrthkn(ZZ,lrintrf=True):
  """
  Analyze layer interface depths for large jumps
  ZZ - default, interface depths lrintrf
       if not - layer thickness (m)
  """
  if not lrintrf:
    dZZ = ZZ.copy()
  else:
    dZZ = np.abs(np.diff(ZZ, axis=0))

  kdm = dZZ.shape[0]
  jdm = dZZ.shape[1]
  idm = dZZ.shape[2]
 
  grdZ = np.zeros((jdm,idm))
  for k in range(kdm):
    aa = np.squeeze(dZZ[k,:,:])
#    aa = np.where(np.isnan(aa),0,aa)
    dzdx = np.zeros((jdm,idm))
    dzdy = np.zeros((jdm,idm))
    for i in range(1,idm):
      dzdx[:,i] = np.abs(aa[:,i-1]-aa[:,i])
    for j in range(1,jdm):
      dzdy[j,:] = np.abs(aa[j-1,:]-aa[j,:])

    dmm   = dzdx-dzdy
    gradz = np.where(dmm>0,dzdx,dzdy)
    dmm   = gradz-grdZ
    grdZ  = np.where(dmm>0,gradz,grdZ)

    print(' Grad lr thkn: k={0} max={1:8.2f} m'.format(k,np.nanmax(grdZ)))
  
  return grdZ

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


def clrmp_lmask(nclrs=2):
  """
    Create colormap for land mask with 2 colors
  """

  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  clrs   = cm.get_cmap('GnBu_r',nclrs)
  newclr = clrs(range(nclrs))
  newclr[0,:] = [1, 1, 1, 1]
  newclr[1,:] = [0.3, 0.3, 0.3, 1]

  newcmp  = ListedColormap(newclr)

  return newcmp



