"""
  Convert unformatted binary files with ocean increments from RTOFS-DA
  into netcdf format required by MOM6 for increments 
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
from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime

dnmb_in  = mtime.datenum([2022,2,15])
dnmb_out = mtime.datenum([2021,7,2])
pthincr  = '/scratch1/NCEPDEV/stmp2/Zulema.Garraffo/FV3_RT/expt_11.9/' + \
          'ncoda/hycom_var/restart/'

dvin       = mtime.datevec(dnmb_in)
_, jday_in = mtime.dnmb2jday(dnmb_in)
rdayin     = f"{dvin[0]}{dvin[1]:02d}{dvin[2]:02d}00"
sfxin      = f"{dvin[3]:04d}"
dvout      = mtime.datevec(dnmb_out)
_, jday_out= mtime.dnmb2jday(dnmb_out)
rdayout    = f"{dvout[0]}{dvout[1]:02d}{dvout[2]:02d}00"
sfxout     = f"{dvout[3]:04d}"

flS = 'salint_pre_1o4500x3297_' + f"{rdayin}_{sfxin}_analinc"
flT = 'seatmp_pre_1o4500x3297_' + f"{rdayin}_{sfxin}_analinc"
flU = 'uucurr_pre_1o4500x3297_' + f"{rdayin}_{sfxin}_analinc"
flV = 'vvcurr_pre_1o4500x3297_' + f"{rdayin}_{sfxin}_analinc"



