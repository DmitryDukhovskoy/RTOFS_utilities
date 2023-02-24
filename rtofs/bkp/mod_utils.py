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


def find_indx_lonlat(x0,y0,LON,LAT):
  """
  Find closest grid point to lon/lat coordinate
  x0,y0 - lon/lat of a point
  """
  if x0 > 180.:
    x0 = x0-360.

  dmm = np.sqrt((LON-x0)**2+(LAT-y0)**2)
  jj0, ii0 = np.where(dmm == np.min(dmm))

  return ii0[0], jj0[0]



