"""
  Read CFSR data
  6-hr data
"""
import os
import sys
import numpy as np
import importlib
import datetime
import time
import netCDF4
from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

atm_frc  = 'CFSR'
YR       = 2021
MM       = 9
DD       = 1
HR       = 0
datenm   = '{0}{1:02d}{2:02d}{3:02d}'.format(YR,MM,DD,HR)
pthout   = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/CFSR/{0}{1:02d}/'.\
           format(YR,MM)
flout    = 'cfsr.{0}.nc'.format(datenm)
fldir    = pthout + flout

nc = netCDF4.Dataset(fldir, 'r')

# Check any missing dates:
import mod_time as mtime

time_in = nc['time'][:].data[0]

# Convert to days YYYYMMDD and timeofday = 0, 21600, ... seconds
Ndays     = time_in/86400.
day_ref   = mtime.datenum([1970,1,1])
dv        = mtime.datevec(day_ref)
Time_dnmb = day_ref + Ndays
DV        = mtime.datevec(Time_dnmb)
SEC       = DV[3]*3600

print('Time from {0}: {1}/{2:02d}/{3:02d} {4:02d}hr {5:02d} min'.\
      format(flout, DV[0], DV[1], DV[2], DV[3], DV[4]))
