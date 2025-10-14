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
MM       = 1
DD       = 25
HR       = 0
SC       = HR*3600
datenm   = '{0}{1:02d}{2:02d}{3:05d}'.format(YR,MM,DD,SC)
pthrst   = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6_003/RESTART/'
flname   = 'DATM_CFSR.cpl.r.{0}-{1:02d}-{2:02d}-{3:05d}.nc'.format(YR,MM,DD,SC)
fdrstrt  = pthrst + flname

nc = netCDF4.Dataset(fdrstrt, 'r')

# Check any missing dates:
import mod_time as mtime

time_in = nc['time'][:].data[0]
day_ref = mtime.datenum([2020,12,30]) # end of previous cycle
dv        = mtime.datevec(day_ref)
Time_dnmb = day_ref + time_in
DV        = mtime.datevec(Time_dnmb)

print('Restart Time from {0}: {1}/{2:02d}/{3:02d} {4:02d}hr {5:02d} min'.\
      format(flname, DV[0], DV[1], DV[2], DV[3], DV[4]))
