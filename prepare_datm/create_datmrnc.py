"""
  For continued set of cycles, new year of atm data
  need to change DATM_*.datm.r.*nc file
  which is a list of dataes/times for merged atm. file
  created from GFS, CFSR or other atm. files 

  The time stamp file is created based on the datm input file
  with merged h-hr atm. data
  e.g.: cfsr_202101_202201.nc

  Check created datm file:
  check_datm.py
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

#pthinp   = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6_003/'
#flinp    = 'DATM_CFSR.datm.r.2020-12-30-00000.nc'
pthout   = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6_003/'
flout    = 'DATM_CFSR.datm.r.2021-01-25-00000.nc'
# File with merged atm forcing fields to get time stamps from:
pthdata  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/CFSR/'
fldata   = 'cfsr_202101_202201.nc'
char_len = 50

#fldatm_in  = pthinp + flinp
fldatm_out = pthout + flout
flatm_flds = pthdata + fldata

print('Creating ' + fldatm_out)
print('using ' + flatm_flds)

# Create new restart from template for writing CICE4 fields
import mod_datm_utils as mdatm
importlib.reload(mdatm)
#mdatm.datm_newfile(fldatm_in, fldatm_out)
try:
  os.remove(fldatm_out)
except:
  pass

import mod_time as mtime
# Derive time stamps:
# Convert to days YYYYMMDD and timeofday = 0, 21600, ... seconds
atmnc     = ncFile(flatm_flds, 'r')
Time      = atmnc['time'][:].data
Ndays     = Time/86400.
day_ref   = mtime.datenum([1970,1,1])
dv        = mtime.datevec(day_ref) 
Time_dnmb = day_ref + Ndays
DV        = mtime.datevec2D(Time_dnmb)
SEC       = DV[:,3]*3600
SEC       = SEC.astype(int)

# form date values:
nrec = DV.shape[0]
date_rec = np.zeros((1,1,nrec), dtype=int)
for ii in range(nrec):
  date_rec[0,0,ii] = int(DV[ii,0]*10000 + DV[ii,1]*100 + DV[ii,2])
 
str_info = 'DATM_INPUT/' + fldata
frmt = 'S{0}'.format(char_len) 
for ii in range(len(str_info), char_len):
  str_info = str_info + ' '

strval = np.array(['DATM_INPUT/' + fldata],frmt)
ncstr  = np.squeeze(netCDF4.stringtochar(strval))

print('Writing data to ' + fldatm_out)
# Create new nc file:
nc = netCDF4.Dataset(fldatm_out, 'w')

# Glbal attributes:
nc.title       = 'Atmospheric forcing time stamps'
nc.source      = '/home/Dmitry.Dukhovskoy/python/prepare_datm/create_datmrnc.py'
nc.institution = 'NOAA NWS EMC'
nc.history     = '{0} creation of DATM time stamps  netcdf file'.format(
                 datetime.datetime.now().strftime("%Y-%m-%d"))
nc.version     = 'nuopc_data_models_v0'

# Dimensions:
strlen   = nc.createDimension('strlen', char_len)
ntime    = nc.createDimension('nt', nrec)
nfiles   = nc.createDimension('nfiles', 1)
nstreams = nc.createDimension('nstreams', 1)

# Variables:
lb_var   = nc.createVariable('ymdLB', np.int32, ('nstreams'))
ub_var   = nc.createVariable('ymdUB', np.int32, ('nstreams'))
tlb_var  = nc.createVariable('todLB', np.int32, ('nstreams'))
tub_var  = nc.createVariable('todUB', np.int32, ('nstreams'))
nfl_var  = nc.createVariable('nfiles', np.int32, ('nstreams'))
ofst_var = nc.createVariable('offset', np.int32, ('nstreams'))
klvd_var = nc.createVariable('k_lvd', np.int32, ('nstreams'))
nlvd_var = nc.createVariable('n_lvd', np.int32, ('nstreams'))
kgvd_var = nc.createVariable('k_gvd', np.int32, ('nstreams'))
ngvd_var = nc.createVariable('n_gvd', np.int32, ('nstreams'))
nt_var   = nc.createVariable('nt', np.int32, ('nstreams','nfiles'))
hvdat_var= nc.createVariable('haveData', np.int32, ('nstreams','nfiles'))
flnm_var = nc.createVariable('filename', 'S1', ('nstreams','nfiles','strlen'))
date_var = nc.createVariable('date', np.int32, ('nstreams', 'nfiles', 'nt'))
tday_var = nc.createVariable('timeofday', np.int32, ('nstreams', 'nfiles', 'nt'))

# Load values
nfl_var[:]   = 1
ofst_var[:]  = 0
klvd_var[:]  = 1
nlvd_var[:]  = 1
kgvd_var[:]  = 1
ngvd_var[:]  = 1464
nt_var[:]    = nrec
hvdat_var[:] = 1
flnm_var[:]  = netCDF4.stringtochar(ncstr)
date_var[:]  = date_rec
tday_var[:]  = SEC

nc.close() 

print('All done')

