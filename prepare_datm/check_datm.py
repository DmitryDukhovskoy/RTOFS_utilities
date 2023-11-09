"""
  Read DATM_????.datm.r.* files -  time stamps for data atm. forcing
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
date     = '2021-01-25-00000'
#date     = '2020-12-04-00000'
pthout   = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6_003/'
flout    = 'DATM_{0}.datm.r.{1}.nc'.format(atm_frc,date)
# File with merged atm forcing fields to get time stamps from:

fldatm_out = pthout + flout
nc = netCDF4.Dataset(fldatm_out, 'r')

date_yr = nc['date'][:].data
hr_day  = nc['timeofday'][:].data 

print('data shape: ')
print(date_yr.shape)
date_yr = date_yr.squeeze()

print('hourofday shape: ')
print(hr_day.shape)
hr_day = hr_day.squeeze()

# How many records per day:
II = np.where(hr_day == 0)[0]
dI = np.diff(II)
nrc_day = dI[0]
print('Number of records per 1 day: {0}'.format(nrc_day))
if nrc_day != np.min(dI) or nrc_day != np.max(dI):
  print('Check records, some are missing: min/max nrcords={0}/{1}'.\
        format(np.min(dI), np.max(dI)))

# Check any missing dates:
import mod_time as mtime
d1 = date_yr[0]
d2 = date_yr[-1]

def datestr2dnmb(datestr):
  """
  Convert YYYYMMDD integer date string to date number 
  """
  yr1   = int(np.floor(datestr/10000))
  mm1   = int(np.floor((datestr-yr1*10000)/100))
  dd1   = datestr - (yr1*10000+mm1*100)
  dnmb1 = mtime.datenum([yr1,mm1,dd1])
  return dnmb1

dnmb1 = int(datestr2dnmb(d1))
dnmb2 = int(datestr2dnmb(d2))

nrec      = len(date_yr)
nrec_true = (dnmb2 - dnmb1 + 1)*nrc_day 
print('Missing records: {0}'.format(nrec_true-nrec))

icc = -1
for dnmb in range(dnmb1,dnmb2+1):
  for ihr in range(4):
    icc += 1
    dstr0 = date_yr[icc]
    dnmb0 = datestr2dnmb(dstr0)
    dv0   = mtime.datevec(dnmb0)
    dv    = mtime.datevec(dnmb)
    if dnmb0 != dnmb:
      print('Date expected: {0}/{1}/{2} dayrec={3}'.\
      format(dv[0],dv[1],dv[2],ihr))
      print('Date in file: {0}/{1}/{2} dayrec={3}'.\
      format(dv0[0],dv0[1],dv0[2],ihr))
#      raise Exception ('Date mismatch, missing data')
      icc -= 1


