"""
  Create ocean mask netCDF file
  ocean fraction at T-cell centers  ( h-points)

"""
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
import netCDF4
from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')


pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/GDAS/global-workflow/fix/mom6/008/'
flout   = 'ocean_mask008.nc'

flmsk_out = pthout + flout
print('Creating ' + flmsk_out)

Not finished - not needed, found ocean mask file for 0.08



