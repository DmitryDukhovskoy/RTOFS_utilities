"""
  Analize T/S profiles in the Arctic 
  rejected during QC

  Previous Profile Observations - downloaded from WOD18 website
  
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
from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text

import mod_utils as mutil
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_read_ncoda as rncoda
importlib.reload(rncoda)
#import mod_extrargo as exargo
#importlib.reload(exargo)
from mod_utils import tsarray

np.set_printoptions(precision=3)

# WOD data
pthwod = '/scratch1/NCEPDEV/stmp4/Dmitry.Dukhovskoy/WOD_profiles/'

rtofs_outp = 'incup'
#rtofs_outp = 'bkgrd'
rdate0     = '20230129'  # forecast date with Argo is assimilated

yr, mo, mday, hr = rncoda.parse_rdate(rdate0)
time_n00    = datetime.datetime(yr,mo,mday,hr,0,0)  # forecast time
# incup:
time_n24    = time_n00 - datetime.timedelta(days=1) # time of hindcast incup incorpor
hdate_n24   = time_n24.strftime('%Y%m%d%H')
rdate_n24   = rdate0

# background, analysis, and Argo obs for assimilation:
time_bkgrd  = time_n00 - datetime.timedelta(days=1)
hdate_bkgrd = time_bkgrd.strftime('%Y%m%d%H')
rdate_bkgrd = time_bkgrd.strftime('%Y%m%d')
rdate_anls  = rdate_bkgrd    # analysis time
yrA, moA, mdayA, hrA = rncoda.parse_rdate(rdate_anls)

# 12hr f/cast:
time_f12    = datetime.datetime(yr,mo,mday,hr,12,0)
hdate_f12   = time_f12.strftime('%Y%m%d%H')
rdate_f12   = rdate0

pthbin = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+\
          rdate0+'/ocnqc_logs/profile_qc/'
ffinal = pthbin+'TSargo_stat_'+rdate_anls+rtofs_outp+'.pkl'

print('Loading '+ffinal)
with open(ffinal,'rb') as fid:
  TSERR = pickle.load(fid)

# Find NCODA consecutive #s of rejected profiles
Srj  = TSERR.Sqcrj.astype(int)
Trj  = TSERR.Tqcrj.astype(int)
rnmb = TSERR.recn.astype(int)
Srj = Srj[np.where(Srj >= 0)]
Trj = Trj[np.where(Trj >= 0)]
ns  = Srj.shape[0]
nt  = Trj.shape[0]


# Indices of T or S rejected profile - total
TSrj = Srj.copy()
for ii in range(nt):
  smm = abs(Srj-Trj[ii])
  if min(smm) > 0:
    TSrj = np.append(TSrj,Trj[ii])

# Argo numbers:
Anmb = TSERR.numb
# Rejected Argo numbers:
Anmb_rj = Anmb[TSrj]
# Histogram of 
# Rejected by lat:
Lat      = TSERR.lat
Lon      = TSERR.lon
Lat_rjct = Lat[TSrj]
Lon_rjct = Lon[TSrj]
latArc   = 70.
Iarc     = np.where(Lat_rjct >= latArc)
LatArj   = Lat_rjct[Iarc]
LonArj   = Lon_rjct[Iarc]
AnmbArj  = Anmb_rj[Iarc]

# Search region:
dx  = 40.   # lon range
dy  = 2.    # lat range
x0  = -132.8
y0  = 78.33 
tm0 = datenum([yrA,moA,mdayA])
# Specify 1st - last month/days
# where to search for data
mo1 = 10
mo2 = 3 

import mod_WODdata as mwod
importlib.reload(mwod)

# Select lon, lat to search for WOD profiles
obtype = 'CTD'
furl  = '{0}{1}/info_{1}.nc'.format(pthwod,obtype)
LAT, LON, TM, UID = mwod.search_UID(furl,x0,y0,dx,dy,mnth1=mo1,mnth2=mo2)
DV = mmisc.datevec1D(TM)

obtype = 'DRB'
furl  = '{0}{1}/info_{1}.nc'.format(pthwod,obtype)
LAT, LON, TM, UID = mwod.search_UID(furl,x0,y0,dx,dy,mnth1=mo1,mnth2=mo2)
DV = mmisc.datevec1D(TM)

obtype = 'PFL'
furl  = '{0}{1}/info_{1}.nc'.format(pthwod,obtype)
LAT, LON, TM, UID = mwod.search_UID(furl,x0,y0,dx,dy,mnth1=mo1,mnth2=mo2)
DV = mmisc.datevec1D(TM)

f_plt = False
if f_plt:
  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()

  ax1 = plt.axes([0.1, 0.25, 0.8, 0.7])
 
  ax1.plot(LonArj,LatArj,'.')
  ax1.plot(x0,y0,'r*')
  ax1.plot([x0-dx,x0-dx],[y0-dy,y0+dy],'r-') 
  ax1.plot([x0+dx,x0+dx],[y0-dy,y0+dy],'r-') 
  ax1.plot([x0-dx,x0+dx],[y0-dy,y0-dy],'r-') 
  ax1.plot([x0-dx,x0+dx],[y0+dy,y0+dy],'r-') 





