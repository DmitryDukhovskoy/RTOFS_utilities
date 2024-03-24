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

def find_gridpnts_box(x0, y0, RR, XX, YY, ioffset=0, joffset=0, rr0=-1e30):
  """
    Given pnt (x0,y0) find 4 grid points enclosing the pnt
    rr0 - some radius within which the 4 points are being searched
    this is typically ~ grid diagonal distance
    can specify a single number (then RR is ignored) or
    provied a 2D RR (say RR = diagonal of grid cells)
    shape RR = shape XX = shape YY

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

  """
  import mod_misc1 as mmisc1

  DD = mmisc1.dist_sphcrd(y0,x0,YY,XX)
  jmin, imin = np.where(DD == np.min(DD))
  jmin = jmin[0]
  imin = imin[0]

  if rr0 <= 0.: 
    rr0  = RR[jmin,imin]

# Find 2 closest grid points first:
  JM, IM = np.where(DD <= rr0)
  DR     = DD[JM,IM]
  isort  = np.argsort(DR)[0:2]
  ixx    = IM[isort]
  jxx    = JM[isort]

#  Find 2 other grid points to enclose the grid pnt x0,y0:
  jv  = jxx[0]
  iv  = ixx[0]
  A   = [XX[jv,iv],YY[jv,iv]]
  jv  = jxx[1]
  iv  = ixx[1]
  B   = [XX[jv,iv],YY[jv,iv]]
  C   = [x0,y0]
  Dr0 = mmisc1.orientation(A,B,C)
  if ixx[0] == ixx[1]:
# 2 pnts along Y axis have been selected find the other 2 + or - dx:
    jv  = jxx[0]
    iv  = ixx[0]-1
    C   = [XX[jv,iv],YY[jv,iv]]
    Dr  = mmisc1.orientation(A,B,C)
    if np.sign(Dr) == np.sign(Dr0):
      di = -1
    else:
      di = 1
    ixx = np.append(ixx,ixx[0]+di)
    jxx = np.append(jxx,jxx[0])
    ixx = np.append(ixx,ixx[1]+di)
    jxx = np.append(jxx,jxx[1])
  else:
# 2 pnts along X axis have been selected find the other 2 + or - dy:
    jv = jxx[0]-1
    iv = ixx[0]
    C   = [XX[jv,iv],YY[jv,iv]]
    Dr  = mmisc1.orientation(A,B,C)
    if np.sign(Dr) == np.sign(Dr0):
      di = -1
    else:
      di = 1
    ixx = np.append(ixx,ixx[0])
    jxx = np.append(jxx,jxx[0]+di)
    ixx = np.append(ixx,ixx[1])
    jxx = np.append(jxx,jxx[1]+di)

  if ioffset > 0:
    ixx = ixx+ioffset
  if joffset > 0:
    jxx = jxx+joffset

  ixx.astype(int)
  jxx.astype(int)

  return ixx, jxx

