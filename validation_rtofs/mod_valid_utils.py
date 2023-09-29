"""
  Functions for validation codes
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
import mod_read_hycom as mrhycom
#import mod_write_hycom as mwhycom
import mod_utils_fig as mufig
import mod_misc1 as mmisc
importlib.reload(mmisc)

def ssh_regions():
  """
    Regions for SSH statistics
    Indices are 0.25-degree CMEMS grid of SSH field
  """
  REGNS = {
    "GOM"     : {
      "xl1"   : 325,
      "xl2"   : 398,
      "yl1"   : 430,
      "yl2"   : 484,
      "Ip"    : [322,388,397,398,392,387,383,370,355,334,323,323],
      "Jp"    : [483,484,461,450,451,451,448,444,430,430,448,471]
    },
    "Carib"   : {
      "xl1"   : 360,
      "xl2"   : 500,
      "yl1"   : 384,
      "yl2"   : 460,
      "Ip"    : [355,370,383,387,392,398,404,416,434,453,473,475,435,424,413,405,401,395,378,367,355],
      "Jp"    : [430,444,448,451,451,450,448,441,436,432,423,395,402,402,388,396,397,392,406,418,426]
    },  
    "GulfStr" : {
      "xl1"   : 380,
      "xl2"   : 570,
      "yl1"   : 455,
      "yl2"   : 550,
      "Ip"    : [388,412,429,454,470,487,498,507,562,564,397],
      "Jp"    : [483,517,536,546,541,552,554,549,548,464,461]
    },
    "Kurosh"  : {
      "xl1"   : 1200,
      "xl2"   : 1339,
      "yl1"   : 420,
      "yl2"   : 630,
      "Ip"    : [1200,1200,1339,1339],
      "Jp"    : [420,625,625,420]
    },
    "Agulhas" : {
      "xl1"   : 734,
      "xl2"   : 942,
      "yl1"   : 169,
      "yl2"   : 317,
      "Ip"    : [734,734,942,942],
      "Jp"    : [169,317,317,169]
    },
    "SOcean"  : {
      "xl1"   : 0,
      "xl2"   : 1339,
      "yl1"   : 0,
      "yl2"   : 230,
      "Ip"    : [0,0,1339,1339],
      "Jp"    : [0,237,237,0]
    },
  }

  return REGNS

def read_ssh_cmems(day_read):
  """
    Read daily SSH CMEMS 
    day_read - matlab format of day #
  """
  import mod_time as mtime
  from netCDF4 import Dataset as ncFile 

  rdate_cmems = mtime.dnumb2rdate(day_read,ihours=False)

  DV = mtime.datevec(day_read)
  YR = DV[0]
  MM = DV[1]
  MD = DV[2]

  pthssh = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/OBS/ssh_altimetry/'
  if MM == 3 or MM == 4:
    flnm   = 'ssh_duacs_nrt_global_20230301_20230430.nc'
  elif MM == 5:
    flnm   = 'ssh_duacs_nrt_global_20230501_20230531.nc'

  flssh  = pthssh + flnm

  nc    = ncFile(flssh,'r')
  TMA   = nc.variables['time'][:].squeeze().data
  dnmbR = mtime.rdate2datenum('19500101')  # ref day 
  TMA   = TMA + dnmbR
  TMA   = TMA.astype(np.int32)

#  print('TMA= {0}'.format(TMA[0]))
#  print('Date {0}/{1}/{2}'.format(YR,MM,MD))
  try:
    iTime = np.where(TMA == day_read)[0][0]
  except:
    print('Day {0} not found in SSH file {1}'.format(day_read,flssh))
    raise Exception('Day not found')

  dnmbA = TMA[iTime]
  SSHA = nc['adt'][iTime,:,:]

  return SSHA



