"""
# Subroutines for Rossby radius calculation
#
# Dmitry Dukhovskoy, NOAA NESDIS NCEI
# July 2022
#
"""
import importlib
import timeit
import binascii
import struct
import sys
import pdb
import numpy as np
#import netCDF4
from netCDF4 import Dataset as ncFile
from os.path import exists
#
# WOD reading modules:
import mod_utils_wod as uwod
#importlib.reload(uwod)
#from mod_utils_wod import read_header

#
def find_recs_prcs(prbnm,YR1=0,YR2=0,lon1=-999.,lon2=-999.,\
                    lat1=-999.,lat2=-999.):
# Specify time, lon/lat range to limit the selected profiles
# if needed
  print('Finding sequential records for {0} {1}-{2}'.format(prbnm,YR1,YR2))
  ion = 1
  krcrd = 0
  Lst_rcrd = []
  while ion > 0:
    krcrd += 1

    if (krcrd % 100000) == 0:
      print(' processed rec# {0}'.format(krcrd))
    
    PRBINFO, PRB = uwod.read_header(prbnm,krcrd) 
    if not PRB:
      print('{0} EoF reached, tot records={1}'.format(prbnm,krcrd))
      ion = 0
      break
#
# Time limts:
    if YR1 > 1800 and YR2 > 1800:
      yr = PRB.yr
# Assumption is that in WOD sequential records follow the year
# quit the loop if the time of the record is newer than requested years
      if yr > YR2:
        break

      if not (yr >= YR1 and yr <= YR2):
        continue
       
#
# lon/lat range: note has not been tested
    if abs(lon1) <= 360. and abs(lon2) <= 360. and\
       abs(lat1) <= 90. and abs(lat2) <= 90.:
      if lon1 > 180.:
        lon1 = lon1-360.
      if lon2 > 180:
        lon2 = lon2-360.
      if lon2 < lon1:
        dm1 = lon1
        dm2 = lon2
        lon1 = dm2
        lon2 = dm1

      lon = PRB.lon
      lat = PRB.lat
      if lon > 360.:
        lon = lon-180.

      if not ((lon >= lon1 and lon <=lon2) and (lat >= lat1 and lat <= lat2)):
        continue

# Check that both profiles are available:
    add_cast = 1
    for varnm in ["Temperature", "Salinity"]:
      vindx = uwod.find_var_indx(PRB,varnm) 
      pntrv = PRB.pobs[vindx]
      if pntrv < 0:
        add_cast = 0

    if add_cast:
#      print('Record in: {0} yr={1} p1={2} p2={3}'.\
#           format(krcrd,PRB.yr,PRB.pobs[0],PRB.pobs[1]))
      Lst_rcrd.append(krcrd)

  print(' {0} records found for {1}-{2}\n'.format(len(Lst_rcrd),YR1,YR2))
  return Lst_rcrd

def read_field(furl,varnm):
  """
    Read OpenDAP netcdf file
  """
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].squeeze()
  dmm = np.copy(dmm0)
  return dmm

def read_1prof(furl,varnm,j0,i0):
  """
    Read 1 vertical profile for specified location
    Note indices in netcdf run 0,... - same as in Python
  """
  nc=ncFile(furl)
  dmm0 = nc.variables[varnm][0,:,j0,i0].squeeze()
  dmm = np.copy(dmm0)
  return dmm

def get_dim(furl):
  """
    Get dimensions for grid/lat vavriables 
    from OpenDAP netcdf WOA files
  """
  nc=ncFile(furl)
  kdm = nc.dimensions['depth'].size
  jdm = nc.dimensions['lat'].size
  idm = nc.dimensions['lon'].size

  return idm, jdm, kdm


def read_WOAgrid025(woa='woa18', grd=0.25, seas=13):
  """
    Read WOA18 0.25 grid and depths
  """
  urlT='https://www.ncei.noaa.gov/thredds-ocean/'+\
       'dodsC/ncei/woa/temperature/decav/{0}/'.\
       format(grd)
  urlS='https://www.ncei.noaa.gov/thredds-ocean/'+\
       'dodsC/ncei/woa/salinity/decav/{0}/'.\
       format(grd)

  if grd==0.25:
    cgrd=4

  tfnm='{0}_decav_t{1:02d}_{2:02d}.nc'.format(woa,seas,cgrd)
  sfnm='{0}_decav_s{1:02d}_{2:02d}.nc'.format(woa,seas,cgrd)

  furl = urlT+tfnm
  ZZ  = read_field(furl,'depth')
  ZZ  = -abs(ZZ)
  lat = read_field(furl,'lat')
  lon = read_field(furl,'lon')

  return lon, lat, ZZ


def read_topo_WOA025(idm=1440, jdm=720):
  """
    Read ETOPO2 interpolated onto 0.25 WOA grid
    + lon/lat
    see interp_etopo2woa.py
  """
  pthout = '/data/ncei2/work/dmitry.dukhovskoy/topo/woa/'
  fout = pthout+'etopo2woa025_bilinear.dat'
  ijdm = idm*jdm
  fid = open(fout,'rb')
  try:
    fid.seek(0)
    HW = np.fromfile(fid,dtype='<f4', count=ijdm)
  finally:
    fid.close()
  HW = HW.reshape(jdm,idm)
  
  fout = pthout+'lonWOA025.dat'
  fid = open(fout,'rb')
  try:
    fid.seek(0)
    LONW = np.fromfile(fid,dtype='<f4', count=idm)
  finally:
    fid.close()
 
# Make sure lon is -180,180
  LONW=np.where(LONW>180.,LONW-360.,LONW)
 
  fout = pthout+'latWOA025.dat'
  fid = open(fout,'rb')
  try:
    fid.seek(0)
    LATW = np.fromfile(fid,dtype='<f4', count=jdm)
  finally:
    fid.close()

  return HW,LONW,LATW

def find_depth(HW,LONW,LATW,lon0,lat0,dltX0=0.25):
# Longitudes are in the (-180,180) range
  if lon0 > 180.:
    lon0 = lon0-360.

  dx = np.abs(LONW-lon0)
  dy = np.abs(LATW-lat0)
  ix0 = np.where(dx==np.min(dx))[0][0] 
  jx0 = np.where(dy==np.min(dy))[0][0] 
  hb = HW[jx0,ix0]

# Sanity check:
# If lon or lat is out of range - 
# Let it pass, this should be caught at QC-check step
  if (abs(lon0) < 360. and abs(lat0) <= 90.) and \
      (np.min(dx) > dltX0 or np.min(dy) > dltX0):
    print('HEADER ERR: lon={0:.3f} lat={1:.3f}'.format(lon0,lat0))
    raise NameError('ERR: dX={1:.3f} or dY={2:.3f} > dltX0={0:.3f} ETOPO2 '.\
      format(dltX0,np.min(dx),np.min(dy))) 
#    print('HEADER ERR: dX={1:.3f} or dY={2:.3f} > dltX0={0:.3f} ETOPO2 '.\
#      format(dltX0,np.min(dx),np.min(dy))) 
#    hb = 9.e10

  return hb


def add_WOA(Tin,Sin,Zin,PRB,LONW,LATW,ZW):
  """
    Add WOA seasonal climatology to observed profiles
    for missing deep ocean section
    This reduces the error in Rossby R. calculation
    In some cases, mismatch ETOPO2 topography and 
    local bottom depth in observations, if WOA depth is 
    shallower than the deepest observed depth  - WOA is not added
    
    Input: T,S,Z - observed values
           zbtm - local bottom depth
           LONW,LATW,ZW - grid/ depth of WOA18
  """
  grd=0.25
  woa='woa18'
  urlT='https://www.ncei.noaa.gov/thredds-ocean/'+\
       'dodsC/ncei/woa/temperature/decav/{0}/'.\
       format(grd)
  urlS='https://www.ncei.noaa.gov/thredds-ocean/'+\
       'dodsC/ncei/woa/salinity/decav/{0}/'.\
       format(grd)

  if grd==0.25:
    cgrd=4

  mo = PRB.mo

  if mo >=1 and mo <=4:
    seas = 13
  elif mo > 4 and mo <= 8:
    seas = 14
  else:
    seas = 15

  tfnm='{0}_decav_t{1:02d}_{2:02d}.nc'.format(woa,seas,cgrd)
  sfnm='{0}_decav_s{1:02d}_{2:02d}.nc'.format(woa,seas,cgrd)

  lon0 = PRB.lon
  lat0 = PRB.lat
# Find WOA grid index corresponding to obs coord:
  dmm = abs(LATW-lat0)
  jplt = np.argmin(dmm)
  dmm = abs(LONW-lon0)
  iplt = np.argmin(dmm)

  print('>> Adding {0} T/S prof: {1:0.2f}N {2:0.2f}E season={3}'.\
       format(woa,LATW[jplt],LONW[iplt],seas))

# Get T/S/Z profiles from WOA:
  furlT=urlT+tfnm
  furlS=urlS+sfnm

  dmm = read_1prof(furlS,'s_an',jplt,iplt)
  dmm[dmm>1.e10]=np.nan
  Swoa = dmm[np.where(~np.isnan(dmm))]

  dmm = read_1prof(furlT,'t_an',jplt,iplt)
  dmm[dmm>1.e10]=np.nan
  Twoa = dmm[np.where(~np.isnan(dmm))]

  kobs = Twoa.shape[0]
  Zwoa = ZW[:kobs]
  if kobs == 0:
    Zwoa = np.array([100.])  # land point

# Add WOA climatology to missing deep section of 
# of the observed casts
  zlst = Zin[-1]
  Iwoa = np.where(Zwoa < zlst)[0]

  if Iwoa.shape[0] == 0:
    print('WARNING: WOA ETOPO bottom {0:.1f} > obs. depth {1:.1f}, skipping WOA'.\
          format(np.min(Zwoa),zlst))
    return Tin, Sin, Zin
#    raise NameError('Empty WOA depth array') 

  Tupdt = np.concatenate([Tin,Twoa[Iwoa]])
  Supdt = np.concatenate([Sin,Swoa[Iwoa]])
  Zupdt = np.concatenate([Zin,Zwoa[Iwoa]])

  return Tupdt, Supdt, Zupdt

def run_QC(Zobs,Tobs,Sobs,zbtm,PRB, f_addWOA=0):
  """
    QC checks - the function returns QC > 0 if there is a problem
    QC flags are not written in order and can be moved around
    I tried to keep QC flags that caused by errors in the input
    header information at the beginning

    Not all QC flags are checked, the first that raised an error
    stops the subroutine and returns the QC flag

    f_QC = 1:  too few observations
         = 2:  # of obs in T and S profiles is too different
         = 3:  the profiles are too shallow wrt local depth
         = 4:  coastal water < 10m
         = 5:  First obs. depth is too deep missing upper ocean T/S structure
         = 6:  First T/S obs are missing and topmost available are too deep
         = 7:  Bottom values are missing and the next available obs
               to be copied are too far apart 
         = 8:  discard the casts if missing section
               is too big wrt to local depth
         = 9:  T/S profiles have sudden "jumps" due to suspecious
               anomalous values (but within the possible range) 
         = 10: Error in lon/lat from input header, values are out of range
         = 11: T or S profile is empty 
  """
  z_shallow = -20.
  f_QC = 0

# Note: the QC flags are not in order
# (10) Input lon/lat are out of range - errors in header data
  lon0 = PRB.lon
  lat0 = PRB.lat

  if abs(lon0) > 360. or abs(lat0) > 90.:
    f_QC = 10
    return f_QC

# (11) T or S profiles are empty - no obsevations
  nt = Tobs.shape[0]
  ns = Sobs.shape[0]
  if ns == 0 or nt == 0:
    f_QC = 11
    return f_QC

# (1) Check total # of observations, too few obs 
# cannot solve e/value problem
  len_min = 7  # min # of observations in profile
  nz = Zobs.shape[0]
  if nz < len_min:
    f_QC = 1
    return f_QC

# (2) Compliteness:
# both T and S profiles should be ~same length
  nt = np.where(~np.isnan(Tobs))[0].shape[0]
  ns = np.where(~np.isnan(Sobs))[0].shape[0]
  reps = 0.2
  err = abs(nt-ns)/(max(nt,ns))
  if err > reps:
    f_QC = 2 
    return f_QC

# (3) The cast is too shallow
# in the deep ocean, it should be at least deeper than 
# the main pycnocline (-500 m)
  z0_deep = -500. # required depth for profile in the deep ocean to pass QC
  deps = 0.25
  zmin = min(Zobs)
  if abs(zmin) < abs(zbtm):
    derr = (1.-abs(zmin)/abs(zbtm))
    if abs(zbtm) < 500.:     # shelf
      if derr > deps:
        f_QC = 3
        return f_QC
    else:
      z0 = min(abs(z0_deep),0.6*abs(zbtm))
      if abs(zmin) < z0:
        f_QC = 3
        return f_QC

# (4) Coastal region with shallow depth
  if abs(zbtm) <= abs(z_shallow):
    f_QC = 4
    return f_QC

# (5) The profile starts too deep below the surface
#     For the radius calculation, near-surface profiles are critical
#     therefore, profiles that miss the upper ~100 m are not useful
  zsrf = -50. 
  zmax = max(Zobs)
  if abs(zmax) > abs(zsrf):
    f_QC = 5
    return f_QC

# (6) Top values are missing on both profiles
#     if missing, the next existing obs is too deep
#     to be copied over
#     Assumed that the 1st obs point is within the surface layer
#     zsrf, then the next available obs. should also be within this
#     layer for smaller error 
  i1 = min(np.where(~np.isnan(Tobs))[0])
  i2 = min(np.where(~np.isnan(Sobs))[0])
  if i1 > 0 and abs(Zobs[i1]) > abs(zsrf):
    f_QC = 6 
    return f_QC
  if i2 > 0 and abs(Zobs[i2]) > abs(zsrf):
    f_QC = 6
    return f_QC
    
# (7) The last bottom value is missing and the next available obs
#     in the profile to be copied is too far up 
  dZdeep = abs(0.1*zbtm)
  i1 = max(np.where(~np.isnan(Tobs))[0])
  i2 = max(np.where(~np.isnan(Sobs))[0])
  dZ1 = abs(Zobs[-1] - Zobs[i1])
  dZ2 = abs(Zobs[-1] - Zobs[i2])
  if dZ1 > dZdeep or dZ2 > dZdeep and f_addWOA == 0:
    f_QC = 7
    return f_QC

# (8) casts that miss most of the water column
#     have to be discarded 
  Hmin = 0.8*zbtm       # min depth to be observed
  if abs(Zobs[-1]) < abs(Hmin) and f_addWOA == 0:
    f_QC = 8
    return f_QC

# (9) Unrealistically large T/S gradients ("jumps")
#     caused by erroneous T/S measurements due to 
#     instrument misfunction when values are within 
#     the possible range but anomalously high/low compared to 
#     local T/S values
  dTdZ = np.abs(np.diff(Tobs)/np.diff(Zobs))
  dSdZ = np.abs(np.diff(Sobs)/np.diff(Zobs))
  if np.max(dTdZ)>5. or np.max(dSdZ)>15.:
    f_QC = 9
    return f_QC
 
  return f_QC 
        

def interp_gaps_TS(Zobs,Tobs,Sobs,zbtm,add_srf=0):
  """
    Fill nans and (if add_srf=1) expand to the surface
    keep obs on observed depth layers
    When extending profile to the surface:
    Assumed that QC has already checked if available T/S obs
    are not too deep to be extrapolated to the surface (QC flag = 6)
  """ 
  from scipy.interpolate import interp1d

  if zbtm > 0.:
    raise ValueError('Bottom depth > 0 zbtm={0}'.format(zbtm))

  if np.min(Zobs) >= 0.:
    Zobs=-Zobs

  if zbtm > Zobs[-1]:
    raise ValueError('Bottom {0} shallower than last obs depth {1}'.\
       format(zbtm,Zobs[-1]))

# Add surface obs
  zML = -51. 
  i1 = min(np.where(~np.isnan(Tobs))[0])
  i2 = min(np.where(~np.isnan(Sobs))[0])
  if add_srf == 1 and abs(Zobs[0]) > 1.e-12:
    Zobs = np.insert(Zobs,0,0.)
    Tobs = np.insert(Tobs,0,Tobs[i1])
    Sobs = np.insert(Sobs,0,Sobs[i2])

  if add_srf == 0:
    if np.isnan(Tobs[0]):
      Tobs[0] = Tobs[i1]
    if np.isnan(Sobs[0]):
      Sobs[0] = Sobs[i1]

# Fill nans
# at the last point - see QC flag = 7
# It should discard T/S profiles where the bottom gap
# is too big
  i1 = max(np.where(~np.isnan(Tobs))[0])
  i2 = max(np.where(~np.isnan(Sobs))[0])
  if np.isnan(Tobs[-1]):
    Tobs[-1] = Tobs[i1]
  if np.isnan(Sobs[-1]):
    Sobs[-1] = Sobs[i2]

  Zintrp = Zobs

# Depth layers for T and S may start/end at different depths
# Also there could be missing values in the profiles at
# observed depths
  for var in ["Temp","Sal"]:
    if var == "Temp":
      V = Tobs
    else:
      V = Sobs

# Discard points if they are nans
    II = np.where(~np.isnan(V))[0]
    V = V[II]
    Zorg = Zobs[II]

#    fpol = interp1d(abs(Zorg),V,kind='cubic')
    fpol = interp1d(abs(Zorg),V,kind='linear')
    Vintrp = fpol(abs(Zintrp)) 

    if var == "Temp":
      Tintrp = Vintrp
    else:
      Sintrp = Vintrp

  return Tintrp, Sintrp, Zintrp

 
def interp_TS_zlev(Zobs,Tobs,Sobs,zbtm):
  """
    The function fills nans, expands obs to the surface 
    and interpolated onto fixed depth layers
    
    Interpolate T/S observations onto high-res Z-levels
    This makes inversion of matrix A (= 1/N2*A) more stable
    Fixes the problem of missing T or S values at observed depths
    Zobs - observed depths
    Tobs/Sobs - observed T/S 
    zbtm - local bottom depth
   
  """
  from scipy.interpolate import interp1d

  if zbtm > 0.:
    raise ValueError('Bottom depth > 0 zbtm={0}'.format(zbtm))

  if np.min(Zobs) >= 0.:
    Zobs=-Zobs

  if zbtm > Zobs[-1]:
    raise ValueError('Bottom {0} shallower than last obs depth {1}'.\
       format(zbtm,Zobs[-1]))

# Add surface obs
# Assumed that QC has already checked if available T/S obs
# are not too deep to be extrapolated to the surface (QC flag = 6)
  zML = -51. 
  i1 = min(np.where(~np.isnan(Tobs))[0])
  i2 = min(np.where(~np.isnan(Sobs))[0])
  if abs(Zobs[0]) > 1.e-12:
    Zobs = np.insert(Zobs,0,0.)
    Tobs = np.insert(Tobs,0,Tobs[i1])
    Sobs = np.insert(Sobs,0,Sobs[i2])

# Interpolate obs profiles onto z levels 
# Fill nans
  ar1 = np.arange(0.,50,1)
  ar2 = np.arange(50,100,2)
  ar20= np.arange(100,200,5)
  ar3 = np.arange(200,1000,10)
  ar4 = np.arange(1000,2000,20)
  ar5 = np.arange(2000,5000,50)
  ar6 = np.arange(5000,9501,100)

  ZZ = np.concatenate((ar1,ar2,ar20,ar3,ar4,ar5,ar6))
  ZZ = -ZZ

# Construct array of interpolation depths:
  dmm = abs(ZZ - zbtm)
  kE = np.argmin(dmm)
  if ZZ[kE] <= zbtm:
    kE += -1

  Zintrp = ZZ[0:kE+1]
  Zintrp = np.append(Zintrp,zbtm)
 

# Depth layers for T and S may start/end at different depths
# Also there could be missing values in the profiles at
# observed depths
  for var in ["Temp","Sal"]:
    if var == "Temp":
      V = Tobs
    else:
      V = Sobs

# Discard start/end points if they are nans
#    iS = min(np.where(~np.isnan(V))[0])  
#    iE = max(np.where(~np.isnan(V))[0])
    II = np.where(~np.isnan(V))[0]
    V = V[II]
    Zorg = Zobs[II]

# Add ghost points for interpolation
# at the surface and bottom
    zS = Zintrp[0]
    if abs(zS-Zorg[0])>1.e-6:
      Zorg = np.insert(Zorg,0,zS)
      V = np.insert(V,0,V[0]) 

    zE = Zintrp[-1]
    if abs(zE-Zorg[-1])>1.e-6:
      Zorg = np.append(Zorg,zE)
      V = np.append(V,V[-1])

#    pdb.set_trace()
#    fpol = interp1d(abs(Zorg),V,kind='cubic')
    fpol = interp1d(abs(Zorg),V,kind='linear')
    Vintrp = fpol(abs(Zintrp)) 

    if var == "Temp":
      Tintrp = Vintrp
    else:
      Sintrp = Vintrp

  return Tintrp, Sintrp, Zintrp

def update_PRCSD(prbnm,f_QC,krcrd,PRCSD):
# Update dict with Processed sequential records and QC flags
# by probes 
  dmm=PRCSD.get(prbnm)
  dmm.append(krcrd)
  PRCSD[prbnm]=dmm

  qcnm = prbnm+'QC'
  dmm=PRCSD[qcnm]
  dmm.append(f_QC)
  PRCSD[qcnm]=dmm

  return PRCSD

def subset_TS2Z(T0,S0,Z0,Zintrp):
  """
    Interpolate/subset T/S from Z0 onto Zintrp points
  """
  from scipy.interpolate import interp1d

  fpol = interp1d(abs(Z0),T0,kind='cubic')
  Tintrp = fpol(abs(Zintrp))

  fpol = interp1d(abs(Z0),S0,kind='cubic')
  Sintrp = fpol(abs(Zintrp))

  return Tintrp, Sintrp

def write_1output(file_out,value,PRB):
  """
    write 1 output value to a binary file
    write out the resultant value to a file in sys.inf 
    which will contain the values from all probe types. 
    The file in sys.inf contains only the calculated values 
    (a real number). The file is indexed by unique cast number, 
    so position 10029370 in the file will contain the calculated 
    value for unique cast 10029370, the first sequential cast in 
    the CTD database.  
  """
  uniqnmb = PRB.uniqnmb
  strval = struct.pack('f',value)
  nskip = struct.calcsize('f'*(uniqnmb-1)) 
  if exists(file_out):
    fid = open(file_out,'rb+')
  else:
    fid = open(file_out,'wb')

  try:
    fid.seek(nskip)
    fid.write(strval)
  finally:
    fid.close()

  return

def update_rcrds(PRCSD,RCRDS):
  """
    Update list of seq. records for given probe 
    removing those that have been already processed
  """
  print('Updating List of records using list of processed obs:')
  prbR = list(RCRDS.keys())
  keysP = list(PRCSD.keys())
  for prbnm in prbR:
# First check if the probename is in the PRCSD
# skip if not
    try:
      ix = keysP.index(prbnm)
    except:
      ix = -1
      continue

    Lst_in  = RCRDS.get(prbnm)
    Lst_out = PRCSD.get(prbnm)

    print('   Probe {0} processed={1}'.format(prbnm,len(Lst_out)))
    for iq in range(len(Lst_out)):
      nmb = Lst_out[iq]
      if nmb in Lst_in:
        Lst_in.remove(nmb)

    RCRDS[prbnm] = Lst_in

  return RCRDS



