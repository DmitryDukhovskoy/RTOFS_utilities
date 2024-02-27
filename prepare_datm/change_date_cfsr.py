"""
  If some 6-hr file is missing, need to create a fake file
  by copying previous or next file 
  need to change date/time field in the netCDF file
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
import netCDF4
from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

# Time, date to use for new file:
YR       = 2021
MM       = 9
DD       = 1
HR       = 0
date_stmp= '{0}{1:02d}{2:02d}{3:02d}'.format(YR,MM,DD,HR)

pthdata  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/CFSR/{0}{1:02d}/'.\
           format(YR,MM)

ftmplt   = 'cfsr.2021090106.nc'   # file used as a template to create missing file
fnew     = 'cfsr.{0}.nc'.format(date_stmp)
fl_tmplt = pthdata + ftmplt
fl_new   = pthdata + fnew


# Create new restart from template for writing CICE4 fields
import mod_datm_utils as mdatm
importlib.reload(mdatm)
mdatm.datm_newfile(fl_tmplt, fl_new)

import mod_time as mtime
dnmb     = mtime.datenum([YR,MM,DD,HR])
ref_dnmb = mtime.datenum([1970,1,1])
Ndays    = dnmb-ref_dnmb
Nsecs    = Ndays*86400.

mdatm.modify_fld_nc(fl_new,'time',Nsecs)

print('All done')

