"""
  Module for WOD data
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
from netCDF4 import Dataset as ncFile

import mod_misc1 as mmisc

def read_ncfld(furl, varnm):
  """
    Read netcdf file, field varnm
  """
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl,'r')
# lookup a variable
  dmm0 = nc.variables[varnm][:].squeeze()
  dmm = np.copy(dmm0)

  return dmm

def search_UID(furl, x0, y0, dx, dy,YR1=0, YR2=0, \
               mnth1=0, mnth2=0, mday1=0, mday2=0):
  """
   Search for UID - unique station number used at NESDIS
   for profiles within lon/lat range and time
   from month1 day1 to month2 day2
   in year YR1 - YR2
   if YR1=0, YR2=0 - all years
   if mnth1=0, mnth2=0 - all months 

   For all longitudes, select:
     x0 = 0., dx= 180.
  """
  LAT = read_ncfld(furl,'lat')
  LON = read_ncfld(furl,'lon')
  TMM = read_ncfld(furl,'time')
  UID = read_ncfld(furl,'cast')  # WOD unique identifier  

  II = np.where(LON > 180.)[0]
  LON[II] = LON[II]-360.
# Correct time: 
# time:units = "days since 1770-01-01 00:00:00" ;
  dnmb_ref = mmisc.datenum([1770,1,1])
  TM = TMM + dnmb_ref

  IX = np.arange(LON.shape[0])
# Find profiles 
# bounded by lat
  J = np.where( (LAT >= y0-dy) & (LAT <= y0+dy))[0]
  
  LAT = LAT[J]
  LON = LON[J]
  TM  = TM[J]
  IX  = IX[J]
  UID = UID[J]
  
# Bounded by Lon
  if x0 > 180.:
    x0 = x0-360.

  LON360     = LON.copy()
  if dx > 180.:
    dx = 180.

  if x0+dx > 180.:
    II         = np.where(LON < 0.)[0]
    LON360[II] = LON360[II]+360.
    J = np.where((LON >= x0-dx) & (LON2 <= x0+dx))[0]
  elif x0-dx < -180.:
    II         = np.where(LON >= 0.)[0]
    LON360[II] = LON360[II]-360.
    J = np.where((LON >= x0-dx) & (LON2 <= x0+dx))[0]    
  else:
    J = np.where((LON >= x0-dx) & (LON <= x0+dx))[0]

  LAT = LAT[J]
  LON = LON[J]
  TM  = TM[J]
  IX  = IX[J]
  UID = UID[J]

# Take all time: yr, months
  if (YR1==0 or YR2==0) & (mnth1==0 or mnth2==0) & \
     (mday1==0 or mday2==0): 
    return LAT, LON, TM, UID

# Time bounded
  DV = mmisc.datevec1D(TM)
  YY = DV[0]
  MM = DV[1]  # months
  DM = DV[2]  # month days

  if (YR1 > 0) | (YR2 > 0):
    J = np.where((YY >= YR1) & (YY <= YR2))
    MM = MM[J]
    DM = DM[J]
    LAT = LAT[J]
    LON = LON[J]
    TM  = TM[J]
    IX  = IX[J]
    UID = UID[J]

  if (mnth1 > 0) & (mnth2 > 0):
    if mnth1 > mnth2:
      J = np.where((MM >= mnth1) | (MM <= mnth2))
    else:
      J = np.where((MM >= mnth1) & (MM <= mnth2))

    MM = MM[J]
    DM = DM[J]
    LAT = LAT[J]
    LON = LON[J]
    TM  = TM[J]
    IX  = IX[J]
    UID = UID[J]

# Days limits are not finished:

  return LAT, LON, TM, UID

