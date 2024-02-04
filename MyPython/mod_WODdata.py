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

def read_ncfld(furl, varnm, array=False, finfo=True):
  """
    Read netcdf file, field varnm
  """
  if finfo:
    print("Reading {1} from {0}".format(furl,varnm))

  nc=ncFile(furl,'r')
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  try:
    nl = len(dmm0)
    dmm = np.copy(dmm0)
  except:
# scalalr
    dmm = float(dmm0)
    if int(dmm) == dmm:
      dmm = int(dmm)
    if array:
      dmm = np.array([dmm]) # avoid 0-length arrays 

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
    J = np.where((LON >= x0-dx) & (LON360 <= x0+dx))[0]
  elif x0-dx < -180.:
    II         = np.where(LON >= 0.)[0]
    LON360[II] = LON360[II]-360.
    J = np.where((LON >= x0-dx) & (LON360 <= x0+dx))[0]    
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

  if (YR1 > 0) & (YR2 > 0):
    print(f'WOD selecting: years = {YR1} - {YR2}')
    J = np.where((YY >= YR1) & (YY <= YR2))[0]
    print(f'Found {len(J)} records')
#    print(f'YY[1]={YY[J[0]]}')
    MM = MM[J]
    DM = DM[J]
    LAT = LAT[J]
    LON = LON[J]
    TM  = TM[J]
    IX  = IX[J]
    UID = UID[J]

    if len(J) == 0:
      return LAT, LON, TM, UID

  if (mnth1 > 0) & (mnth2 > 0):
    print(f'WOD selecting: months = {mnth1} - {mnth2}')
    if mnth1 > mnth2:
      J = np.where((MM >= mnth1) | (MM <= mnth2))[0]
    else:
      J = np.where((MM >= mnth1) & (MM <= mnth2))[0]
    print(f'Found {len(J)} records')

    MM = MM[J]
    DM = DM[J]
    LAT = LAT[J]
    LON = LON[J]
    TM  = TM[J]
    IX  = IX[J]
    UID = UID[J]

# Days limits are not finished:

  return LAT, LON, TM, UID

def read_WOA18_3D(finp, idm, jdm, kdm, f_info=True):
  """
    Read WOA18 climatology
    Fortran binary files
  """
  print('Reading WOA18 3D: ' + finp)
  ijdm = idm*jdm

  try:
    fga = open(finp, 'rb')
  except:
    print('Could not open '+finp)

  F3d = np.zeros((kdm,jdm,idm)) + 1.e30
  fga.seek(0)
  for kk in range(kdm):
#    nbts = np.fromfile(fga, dtype='>i4', count=1)[0]
    AA = np.fromfile(fga, dtype='>f4', count=ijdm)
    AA = np.reshape(AA,(jdm,idm), order='C')
    F3d[kk,:,:] = AA
    if f_info:
      print(f"Lr {kk}, min/max: {np.min(AA):.4f}/{np.max(AA):.4f}")
  
  return F3d

def woa_profile(finp, LON, LAT, idm, jdm, kdm, lon0, lat0, f_info=False):
  """
    Extract WOA profile for given location
  """
  F3d = read_WOA18_3D(finp, idm, jdm, kdm, f_info)

  if lon0 < 0.:
    II = np.where(LON>180.)[0]
    LON[II] = LON[II]-360.

  dst = np.square((LON-lon0)**2)
  I0  = np.where(dst == np.min(dst))[0][0]
  lng = LON[I0]
  dst = np.square((LAT-lat0)**2)
  J0  = np.where(dst == np.min(dst))[0][0]
  ltg = LAT[J0]
  dst = np.square((lng-lon0)**2+(ltg-lat0)**2)
  if dst > 0.25:
    print(' ERROR finding closest WOA18 pnt for {0:6.3f}N {0:6.3f}E'.\
           format(lat0, lon0))
    raise Exception('Check closest point in WOA18: {0:6.3f}N {0:6.3f}E'.\
           format(ltg, lng))

  Aprf = np.squeeze(F3d[:,J0,I0])

  return Aprf
  
def read_WOA18_grid(fgrid):
  """
    WOA18 Grid
  """
  print('Reading WOA18 grid: ' + fgrid)

  try: 
    fga = open(fgrid, 'rb')
  except:
    print('Could not open '+fgrid)

  fga.seek(0)
  nbts = np.fromfile(fga, dtype='>i4', count=1)[0]
  kdm  = np.fromfile(fga, dtype='>i4', count=1)[0]
  nbte = np.fromfile(fga, dtype='>i4', count=1)[0]

  nbts = np.fromfile(fga, dtype='>i4', count=1)[0]
  ZZ   = np.fromfile(fga, dtype='>f4', count=kdm)
  nbte = np.fromfile(fga, dtype='>i4', count=1)[0]

  nbts = np.fromfile(fga, dtype='>i4', count=1)[0]
  idm  = np.fromfile(fga, dtype='>i4', count=1)[0]
  nbte = np.fromfile(fga, dtype='>i4', count=1)[0]

  nbts = np.fromfile(fga, dtype='>i4', count=1)[0]
  LON  = np.fromfile(fga, dtype='>f4', count=idm)
  nbte = np.fromfile(fga, dtype='>i4', count=1)[0]

  nbts = np.fromfile(fga, dtype='>i4', count=1)[0]
  jdm  = np.fromfile(fga, dtype='>i4', count=1)[0]
  nbte = np.fromfile(fga, dtype='>i4', count=1)[0]

  nbts = np.fromfile(fga, dtype='>i4', count=1)[0]
  LAT  = np.fromfile(fga, dtype='>f4', count=jdm)
  nbte = np.fromfile(fga, dtype='>i4', count=1)[0]

  fga.close()

  return kdm, idm, jdm, ZZ, LON, LAT

def readTS_WODuid(dflnm):
  ZZ = read_ncfld(dflnm, 'z', array=True, finfo=False)
  ZZ = -abs(ZZ)
  S  = read_ncfld(dflnm, 'Salinity')
  T  = read_ncfld(dflnm, 'Temperature')

  return ZZ, T, S

def select_WODdepth(pthdata, UID, Hmin, LON, LAT, HH, \
                    qflag = True, nmin_obs=5, zmin_obs=-100.):
  """
    Subset WOD observations from the list of UID
    based on min bottom depth requirement
    LON, LAT, HH - grid and topography that are used as a reference
    qflag =  discard flagged observations
    nimn_obs - min # of observations in the profile
    also discard too shallow profiles > zmin_obs

    WODf:flag_values = 0s, 1s, 2s, 3s, 4s, 5s, 6s, 7s, 8s, 9s ;
    WODf:flag_meanings = "accepted range_out inversion gradient anomaly gradient+inversion range+inversion range+gradient range+anomaly range+inversion+gradient" ;

  """
  import os
  import mod_utils as mutil

  nrec  = len(UID)
  Jin   = []
  UIDs  = []
  nrej1 = 0
  nrej2 = 0
  nrej3 = 0
  nrej4 = 0
  for ii in range(nrec):
    if ii%200==0:
      print(f' subsetting {float(ii/nrec)*100.:.1f}% done ...')

    flnm = f'wod_0{UID[ii]}O.nc'
    dflnm = os.path.join(pthdata, flnm)

    uid   = read_ncfld(dflnm, 'wod_unique_cast', finfo=False) # WODunique id
    ZZ    = read_ncfld(dflnm, 'z', array=True, finfo=False)  
    ZZ   = -abs(ZZ)
    Hmin = -abs(Hmin)
    zmin = min(ZZ)
    if zmin > zmin_obs or len(ZZ) < nmin_obs:
#      print(f'Min depth obs = {zmin:.3f}')
      nrej1 += 1
      continue 

    sflag = read_ncfld(dflnm, 'Salinity_WODflag', finfo=False) # 0 - ok, 
    tflag = read_ncfld(dflnm, 'Temperature_WODflag', finfo=False)
    T     = read_ncfld(dflnm, 'Temperature', finfo=False)
    S     = read_ncfld(dflnm, 'Salinity', finfo=False)

# Check QC flags:
    IT = np.where(tflag > 3)[0]
    IS = np.where(sflag > 3)[0]

    if len(IT) > 0:
      print(f'{uid}: err flags T  = {len(IT)/len(ZZ)*100:.3f}%, ' +\
            f'S = {len(IT)/len(ZZ)*100:.3f}%')

    len_cutoff = int(0.2*len(ZZ))
    if len(IT) > len_cutoff or len(IS) > len_cutoff:
      nrej2 += 1
      continue

    IT = np.where(tflag <= 3)[0]
    IS = np.where(sflag <= 3)[0]
#    print(f'{uid}: N good flags = {len(IT)/len(ZZ)*100:.3f}% ' +\
#          f'S = {len(IT)/len(ZZ)*100:.3f}%')
    if len(IT) < nmin_obs or len(IS) < nmin_obs:
      nrej3 += 1
      continue

# Check local depth:
    y0    = read_ncfld(dflnm,'lat', finfo=False)
    x0    = read_ncfld(dflnm,'lon', finfo=False)
#    print(f'uid={uid} x0={x0} y0={y0}')
    inx, jnx = mutil.find_indx_lonlat(x0, y0, LON, LAT)
    Hbtm = HH[jnx, inx]
    if Hbtm > Hmin:
      nrej4 += 1
      continue

#    print(f'appending {uid}')
    UIDs.append(uid)
    Jin.append(ii)

  print('===============')
  print(f'Rejected shallow profile: {nrej1}')
  print(f'Rejected QC flag:         {nrej2}')
  print(f'Rejected few good observ: {nrej3}')
  print(f'Rejected shallow bottom:  {nrej4}')

  return Jin, UIDs










 

