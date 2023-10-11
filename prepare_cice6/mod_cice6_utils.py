"""
  Utility subroutines for cice 6
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pdb
import importlib
#import struct
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
#from mpl_toolkits.basemap import Basemap, shiftgrid

#sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
#sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
#sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
#import mod_read_hycom
#importlib.reload(mod_read_hycom)
#from mod_read_hycom import read_grid_topo, read_hycom, read_topo
#from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text


def param_cice4(nx=4500, ny=3298):
  CICEP = {
    "ncat"   : 5,
    "nilyr"  : 4,
    "nslyr"  : 1,
    "nx"     : nx,
    "ny"     : ny
    }

# # of ice layers in all cat
  ntilyr = CICEP["nilyr"]*CICEP["ncat"]
# # of snow layers in all cat
  ntslyr = CICEP["nslyr"]*CICEP["ncat"]  
  CICEP.update({"ntilyr" : ntilyr})
  CICEP.update({"ntslyr" : ntslyr})

  return CICEP

class cice4():
  def __init__(self, nx=4500, ny=3297):
    self.ncat   = 5
    self.nilyr  = 4
    self.nslyr  = 1
    self.nx     = nx
    self.ny     = ny
    self.ntilyr = self.ncat*self.nilyr
    self.ntslyr = self.ncat*self.nslyr

class cice6():
  def __init__(self, nx=4500, ny=3297):
    self.ncat   = 5
    self.nilyr  = 7
    self.nslyr  = 1
    self.nx     = nx
    self.ny     = ny
    self.ntilyr = self.ncat*self.nilyr
    self.ntslyr = self.ncat*self.nslyr

def plot_2d(A, fgnmb=1):
  """
    Quick plot 2D field
  """
  plt.ion()
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im1 = ax1.pcolormesh(A)

def modify_fld_nc(infile,var_name,AA):
  """
    Modify field in netCDF
    netCDF exists (template) 
  """
  ncdata = ncFile(infile,'r+')
  ncdata[var_name][:] = AA
  ncdata.close()
  
  return 

def modify_glattr_nc(infile, var_name, var_value):
  """
    Modify global attribute in restart netCDF
    such as restart time 
  """
  ncdata = ncFile(infile,'r+')
  ncdata.setncattr(var_name, var_value)  
  ncdata.close()

  return

def addnew_glattr_nc(infile, var_name, var_value):
  """
    Add new global attribute in restart netCDF
    such as information about restart files/ authors, etc/

  Note: the code freezes up trying to add a new attribute
  Need to debug why
  """
  ncdata = ncFile(infile,'r+')
  ncdata.setncattr(var_name, var_value)  
  ncdata.close()

  return


def read_cice4_grid(fl_grid, fld_read, IDM=4500, JDM=3297):
  """
   read bathymetry/grid from cice.regional.r 

   in grid2cice.f of HYCOM-tools
   ice grid is written as direct access array:
   each record is IDMo x JDMo double precision (real*8)
          kmt    land mask array (0,1)
          ulati  latitude  of u-cell centers (radians)
          uloni  longitude of u-cell centers (radians)
          htn    length of northern edge of t-cell (m)
          hte    length of eastern  edge of t-cell (m)
          anglet conversion on t-cell between cice and lat-long grids (radians)
          tlati  latitude  of t-cell centers (radians)
          tloni  longitude of t-cell centers (radians)
  """
  print('Reading ' + fld_read + ' from CICE4 grid file' + fl_grid)
  print(' Domain dimensions: IDM={0} JDM={1}'.format(IDM,JDM))

  FLDS = ['kmt','ulati','uloni','htn','hte','anglet','tlati','tloni']
  try:
    iFld = FLDS.index(fld_read)
  except:
    print("Field " + fld_read + " is not in " + fl_grid)

  IJDM = IDM*JDM
  fga  = open(fl_grid,'rb')
  fga.seek(0)

  
  fga.seek(iFld*(8*IJDM),0)
  AA = np.fromfile(fga, dtype='>f8', count=IJDM)
  AA = np.reshape(AA,(JDM,IDM), order='C')

  fga.close()

  return AA

def read_ncfile(flname, fldname, fsilent=False):
  if not fsilent:
    print('reading ' + 'fldname' + ' from ' + flname)
  ncdata = ncFile(flname, 'r')
  AA    = ncdata[fldname][:].data.squeeze()
  return AA 

def add_fld2D_ncfile(flname, fldname):
  """
    Add new variable to existing netCDF
  """
  print('Adding ' + fldname + ' --> ' + flname)
  infile = ncFile(flname, 'r+')
  nx_nc  = infile.dimensions['ni'].name
  ny_nc  = infile.dimensions['nj'].name
  newvar = infile.createVariable('coszen','f8',(ny_nc,nx_nc))
  infile.close() 

  return
 
def grid_rad2dgr(ulat,ulon, f180 = True):
  """
  Convert CICE coordiantes of the grid from
  radians to degrees
  if f180 true - make (-180 <= lon <= 180)
  """
  rdn2dgr = 180./np.pi
  ulat = ulat*rdn2dgr  
  ulon = ulon*rdn2dgr
  if f180:
    ulon = np.where(ulon>180.,ulon-360.,ulon)
    ulon = np.where(ulon<-180.,ulon+360.,ulon)

  ulat = np.where(ulat > 89.99999, 89.99999, ulat)

  return ulat, ulon


def check_cice_grids(ulati4, uloni4, ulati6, uloni6, frad=True, eps0=0.05):
  """
    Check CICE6 and CIC4 lon/lat 
    lon/lat are in radians by default
  
    There seems to be an error in CICE4 regional grid:
    last column (from ~ in lat is repeated (end-1) column
    in CICE6, these are different columns
    max error is 0.05451 in this column and anywhere else < 1.e-5

  """
  print('Checking CICE4 & 6 grids')
  rdn2dgr = 180./np.pi
  if frad:
    ulati4, uloni4 = grid_rad2dgr(ulati4, uloni4)
    ulati6, uloni6 = grid_rad2dgr(ulati6, uloni6)

  DU = abs(ulati4-ulati6)
  DN = abs(uloni4-uloni6)

  print('Max latitude difference |CICE4-CICE6| = {0} dgr'.\
         format(np.max(DU)))
  print('Max longitude difference |CICE4-CICE6| = {0} dgr'.\
         format(np.max(DN)))

  if np.max(DU) and np.max(DN) > eps0:
    print(' Max lat difference exceeds threshold, check CICE4/CICE6 grids')
    print(' Check if both grids are Arakawa B grid ')
    print(' If not - need to add interpolation algorithm to map ')
    print(' CICE4 B grid ---> CICE6 C grid ')
    raise Exception ('STOPPING: CICE4/CICE6 grid mismatch')

  return

def compute_coszen(ulat, ulon, dnmb, frad=True, time_zone=0):
  """
    Compute cos of solar zenith angle (angle between the sun and the vertical)
    Simplified calculation of declination angle (assumption of
    perfect circular orbit of the sun)
    For more accurate - see icepack code icepack_orbital.F90

    ulat, ulon - geogr. coordinates, CICE grid, radians - default
    dnmb - day number in matlab format
     Use NOAA Global Monitoring Division algorithm
     Low accuracy equations
     General Solar Position Calculations
     https://gml.noaa.gov/grad/solcalc/solareqns.PDF 

    Typo fixed in:
    https://github.com/pvlib/pvlib-python/blob/master/pvlib/solarposition.py
    correction for a constant in eqtime (0.000075 should be
    0.0000075)
  """
  import mod_time as mtime

  dgr2rdn = np.pi/180.
  rdn2dgr = 180./np.pi

  if frad:
    ulat = ulat*rdn2dgr
    ulon = ulon*rdn2dgr
    ulon = np.where(ulon>180.,ulon-360.,ulon)
    ulon = np.where(ulon<-180.,ulon+360.,ulon)

  # Solar declination angle:
  DV   = mtime.datevec(dnmb)
  hr   = DV[3]
  if len(DV) > 4:
    mn = DV[4]
    sc = 0.
  else:
    mn = 0.
    sc = 0.

  jday = mtime.date2jday(DV) 
  jd31 = mtime.date2jday([DV[0],12,31])  # # days in this year
  dnmb_wsol = mtime.datenum([DV[0],12,22])
  dnmb_d31  = mtime.datenum([DV[0],12,31])
  days_offset = dnmb_d31 - dnmb_wsol + 1.
# Simple formula for solar declination (degr) for quik check
#  dlt = -23.45*np.cos(dgr2rdn*(360./jd31*(jday + days_offset)))
  dlt0 = -23.45*np.cos(dgr2rdn*(360./365.*(jday + 10.)))

 
# The fractional year (rad):
  gamma = 2.*np.pi/jd31*(jday - 1. + (hr - 12.)/24.)

# Eq. of time (minutes) - note corrected 7.5e-6 instead of 7.5e-5:  
  eqtime = 229.18*(7.5e-6 + 0.001868*np.cos(gamma) - 0.032077*np.sin(gamma) \
           - 0.014615*np.cos(2.*gamma) - 0.040849*np.sin(2.*gamma))

# Solar declination angle (rad), should be comparable to dlt0:
  dlt = 0.006918 - 0.399912*np.cos(gamma) + 0.070257*np.sin(gamma) \
        - 0.006758*np.cos(2.*gamma) + 0.000907*np.sin(2*gamma) \
        - 0.002697*np.cos(3.*gamma) + 0.00148*np.sin(3.*gamma)

# The true solar time:
# Time offset, minutes, longitudes should be in degrees
# time_zone = hours from UTC (US MST = -7hr)
# Here, default - all time is in UTC, time_zone=0
# The factor of 4 minutes comes from the fact that the Earth rotates 1° 
# every 4 minutes
  time_offset = eqtime + 4.*ulon - 60.*time_zone
# True solar time, minutes
  TST = hr*60. + mn + sc/60. + time_offset
# 
# The solar hour angle (degrees):
  SHA = TST/4. - 180.

# Solar zenith angle:
# cos(phi) < 0 - sun below the horizon (e.g., Polar regions in winter)
  coszen = np.sin(ulat*dgr2rdn)*np.sin(dlt) + \
           np.cos(ulat*dgr2rdn)*np.cos(dlt)*np.cos(SHA*dgr2rdn) 

  return coszen

def cice6_newfile(fl_restartT, fl_restart6, fovrd=False):
  """
    Create a new restart CICE6 file from some template file
    if fovrd - rewrite existing fl_restart6
    otherwise - do not do anything
  """
  import shutil

  if os.path.exists(fl_restart6):
    print(fl_restart6 + ' exists, use this file')
    if fovrd:
      print(fl_restart6 + ' will be overode')
    else:
      return 

  print(fl_restartT + ' ---> ' + fl_restart6)
  shutil.copy(fl_restartT, fl_restart6)

  return
 
def read_rcrd_cice4(fid, nx, ny):
  """
    Read 2D record from open file
    secuential binary file
  """
  fdump = np.fromfile(fid, dtype='>i4', count=1)[0]
  A     = np.fromfile(fid, dtype='>f8', count=nx*ny)
  A     = np.reshape(A,(ny,nx), order='C')
  fdump = np.fromfile(fid, dtype='>i4', count=1)[0]

  return A

def read_cice4_restart(fl_restart4, cice4, fld_read):
  """
    Read CICE4 restart fields
    unformatted binary big endian

    Return fld_read = 'aice','vicen','vsnon','qsnon',...

  """
# ----------------------------

  def print_minmax(sfld,A):
    print('   {2} min/max:  {0}/{1}'.format(np.min(A),np.max(A),sfld))

    return
# ----------------------------
  spval = 1.e30

  try:
    fid = open(fl_restart4, 'rb')
  except:
    print('Could not open '+fl_restart4)
    raise Exception('ERR: restart not found')

  print('Reading restart: ' + fl_restart4)
  print('Field to read: ' + fld_read)

  fid.seek(0)
# Read Fortran binary
  fdump   = np.fromfile(fid, dtype='>i4', count=1)[0]
  istep   = np.fromfile(fid, dtype='>i4', count=1)[0]
  runtime = np.fromfile(fid, dtype='>f8', count=1)[0]  # total elapsed time, sec
  frtime  = np.fromfile(fid, dtype='>f8', count=1)[0]  # forcing time, sec
  fdump   = np.fromfile(fid, dtype='>i4', count=1)[0]

  print('Restart: step={0}, total time(yrs)={1}, forcing last update(hrs)={2}'\
         .format(istep,runtime/(3600*24*365.25),frtime/3600.))

# Read state variables:
# Tsfc is the only tracer read in this file
  nx     = cice4.nx
  ny     = cice4.ny
  ncat   = cice4.ncat
  ntilyr = cice4.ntilyr  # total # of icelrs * cat 
  ntslyr = cice4.ntslyr

  aicen = np.zeros((ncat,ny,nx), dtype='float64')
  vicen = np.zeros((ncat,ny,nx), dtype='float64')
  vsnon = np.zeros((ncat,ny,nx), dtype='float64')
  trcrn = np.zeros((ncat,ny,nx), dtype='float64')

  print('Reading state variables: aicen, vicen, vsnon, trcrn')
  for n in range(ncat):
#    print(' Category {0}'.format(n+1))
# Read ice area for category n
    A = read_rcrd_cice4(fid, nx, ny)
    aicen[n,:,:] = A 
#    print_minmax('ice area',A)
  
# Read ice volume/m2 for category n
    A = read_rcrd_cice4(fid, nx, ny)
    vicen[n,:,:] = A
#    print_minmax('ice vol',A)

# Read snow volume/m2 for category n
    A = read_rcrd_cice4(fid, nx, ny)
    vsnon[n,:,:] = A
#    print_minmax('snow vol',A)

# Read tracer 1 (surf T) - only 1 in CICE4
    A = read_rcrd_cice4(fid, nx, ny)
    trcrn[n,:,:] = A
#    print_minmax('surf T',A)

  if fld_read == 'aicen':
    fid.close()  
    return(aicen)

  if fld_read == 'vicen':
    fid.close()  
    return(vicen)

  if fld_read == 'vsnon':
    fid.close()  
    return(vsnon)

  if fld_read == 'surft':
    fid.close()  
    return(trcrn)

# Internal ice layer energy:
  eicen = np.zeros((ntilyr,ny,nx), dtype='float64')
  print('\n Ice energy eicen')

  for k in range(ntilyr):
    A = read_rcrd_cice4(fid, nx, ny)
    eicen[k,:,:] = A
#    print_minmax('{0} eicen'.format(k+1),A) 

  if fld_read == 'eicen':
    fid.close()  
    return(eicen)

# Snow energy:
  esnon = np.zeros((ntslyr,ny,nx), dtype='float64')
  print('\n Snow energy esnon')

  for k in range(ntslyr):
    A = read_rcrd_cice4(fid, nx, ny)
    esnon[k,:,:] = A
#    print_minmax('{0} esnon'.format(k+1),A)

  if fld_read == 'esnon':
    fid.close()  
    return(esnon)

# Velocities:
  print('\n Velocity uvel')
  A = read_rcrd_cice4(fid, nx, ny)
  uvel = A.copy()
#  print_minmax('U vel',A)
  if fld_read == 'uvel':
    fid.close()  
    return(uvel)

  print('\n Velocity vvel')
  A = read_rcrd_cice4(fid, nx, ny)
  vvel = A.copy()
#  print_minmax('V vel',A)
  if fld_read == 'vvel':
    fid.close()  
    return(vvel)

# Radiation fields
# 4 radiative categories
# for calculating albedo for visible and IR wavelengths
# and penetrating sh/wave

# Scale factor to change MKS units
# for shortwave components
# default = 1
  print('\n Radiation fields, W/m2: ')
  A = read_rcrd_cice4(fid, nx, ny)
  scale_factor = A.copy()
#  print_minmax('Scale Factor',A)
  if fld_read == 'scale_factor':
    fid.close()  
    return(scale_factor)

  A = read_rcrd_cice4(fid, nx, ny)
  swvdr = A.copy()
#  print_minmax('Sh/wave down vis. direct',A)
  if fld_read == 'swvdr':
    fid.close()  
    return(swvdr)

  A = read_rcrd_cice4(fid, nx, ny)
  swvdf = A.copy()
#  print_minmax('Sh/wave down vis. diff',A)
  if fld_read == 'swvdf':
    fid.close()  
    return(swvdf)

  A = read_rcrd_cice4(fid, nx, ny)
  swidr = A.copy()
#  print_minmax('Sh/wave down near IR dir',A)
  if fld_read == 'swidr':
    fid.close()  
    return(swidr)

  A = read_rcrd_cice4(fid, nx, ny)
  swidf = A.copy()
#  print_minmax('Sh/wave down near IR diff',A)
  if fld_read == 'swidf':
    fid.close()  
    return(swidf)

# Ocean stress, N/m2
  print('\n Ocean stress components, N/m2:')

  A = read_rcrd_cice4(fid, nx, ny)
  strocnxT = A.copy()
#  print_minmax('ocean stress x-comp',A)
  if fld_read == 'strocnxT':
    fid.close()  
    return(strocnxT)

  A = read_rcrd_cice4(fid, nx, ny)
  strocnyT = A.copy()
#  print_minmax('ocean stress y-comp',A)
  if fld_read == 'strocnyT':
    fid.close()  
    return(strocnyT)

# Internal stress, stress tensor kg/s2
# (1) northeast, (2) northwest, (3) southwest, (4) southeast
  print('\n Internal stress, kg/s2: ')
 
  A = read_rcrd_cice4(fid, nx, ny)
  stressp_1 = A.copy()
#  print_minmax('stressp_1',A)
  if fld_read == 'stressp_1':
    fid.close()  
    return(stressp_1)
 
  A = read_rcrd_cice4(fid, nx, ny)
  stressp_3 = A.copy()
#  print_minmax('stressp_3',A)
  if fld_read == 'stressp_3':
    fid.close()  
    return(stressp_3)

  A = read_rcrd_cice4(fid, nx, ny)
  stressp_2 = A.copy()
#  print_minmax('stressp_2',A)
  if fld_read == 'stressp_2':
    fid.close()  
    return(stressp_2)

  A = read_rcrd_cice4(fid, nx, ny)
  stressp_4 = A.copy()
#  print_minmax('stressp_4',A)
  if fld_read == 'stressp_4':
    fid.close()  
    return(stressp_4)

  A = read_rcrd_cice4(fid, nx, ny)
  stressm_1 = A.copy()
#  print_minmax('stressm_1',A)
  if fld_read == 'stressm_1':
    fid.close()  
    return(stressm_1)
 
  A = read_rcrd_cice4(fid, nx, ny)
  stressm_3 = A.copy()
#  print_minmax('stressm_3',A)
  if fld_read == 'stressm_3':
    fid.close()  
    return(stressm_3)

  A = read_rcrd_cice4(fid, nx, ny)
  stressm_2 = A.copy()
#  print_minmax('stressm_2',A)
  if fld_read == 'stressm_2':
    fid.close()  
    return(stressm_2)

  A = read_rcrd_cice4(fid, nx, ny)
  stressm_4 = A.copy()
#  print_minmax('stressm_4',A)
  if fld_read == 'stressm_4':
    fid.close()  
    return(stressm_4)

  A = read_rcrd_cice4(fid, nx, ny)
  stress12_1 = A.copy()
#  print_minmax('stress12_1',A)
  if fld_read == 'stress12_1':
    fid.close()  
    return(stress12_1)
 
  A = read_rcrd_cice4(fid, nx, ny)
  stress12_3 = A.copy()
#  print_minmax('stress12_3',A)
  if fld_read == 'stress12_3':
    fid.close()  
    return(stress12_3)

  A = read_rcrd_cice4(fid, nx, ny)
  stress12_2 = A.copy()
#  print_minmax('stress12_2',A)
  if fld_read == 'stress12_2':
    fid.close()  
    return(stress12_2)

  A = read_rcrd_cice4(fid, nx, ny)
  stress12_4 = A.copy()
#  print_minmax('stress12_4',A)
  if fld_read == 'stress12_4':
    fid.close()  
    return(stress12_4)

# Ice mask for dynamics
  print('\n Ice Mask for Dynamics: ')

  A = read_rcrd_cice4(fid, nx, ny)
  iceumask = A.copy()
#  print_minmax('iceumask',A)
  if fld_read == 'iceumask':
    fid.close()  
    return(iceumask)

# For trully coupled HYCOM-CICE these fields
# are not needed
# if defined ocean mixed layer in CICE
# This is for ocean mixed layer defined (for GOFS )
  print('\n Ocean mixed layer: \n');

  A = read_rcrd_cice4(fid, nx, ny)
  sst = A.copy()
#  print_minmax('sst',A)
  if fld_read == 'sst':
    fid.close()  
    return(sst)
  
  A = read_rcrd_cice4(fid, nx, ny)
  frzmlt = A.copy()
#  print_minmax('frzmlt',A)
  if fld_read == 'frzmlt':
    fid.close()  
    return(frzmlt)

  fid.close()  

  print(fld_read + ' not found in restart ')
  return

def sub_region2D(AA, region='Arctic'):
  """
    Subsample polar regions from the global grid
    input: AA - 2D scalar field

  """
  imd = AA.shape[1]
  jmd = AA.shape[0]

  if region == 'Arctic':
    jS = 2228
    jE = jmd
  elif region == 'Antarctic':
    jS = 0
    jE = 750
  else:
# Assume Antarctic
    jS = 0
    jE = 750

  AS = AA[jS:jE,:]

  return AS

def colormap_conc():
  """
   Prepare colormap for sea ice conc maps
  """
  import mod_colormaps as mclrs
  CLR = [[238, 226, 215],
         [223, 206, 189],
         [216, 162, 107],
         [208, 131,  54],
         [232, 177,  59],
         [232, 208,  50],
         [250, 241, 110],
         [219, 240,  94],
         [210, 250, 162],
         [157, 246,  97],
         [97,  246, 127],
         [35,  202, 157],
         [122, 238, 246],
         [4,   173, 185],
         [25,  154, 253],
         [8,    80, 174],
         [255, 255, 255]]

  CLR = np.array(CLR)/255.
  CLR = np.flip(CLR, axis=0)
  CMP = mclrs.create_colormap(CLR, 200)

  return CMP

def colormap_enthlp():
  """
   Prepare colormap for ice/snow enthalpy/energy
   Note: enthalpy < 0

   Creates ListedColormap object
   to get colors:
   cmpice.colors
   cmpice.N - # of colors
  """
  import mod_colormaps as mclrs
  CLR = [[255, 255, 255],
         [233, 225, 216],
         [247, 210, 173],
         [208, 131,  54],
         [232, 177,  59],
         [232, 208,  50],
         [250, 241, 110],
         [219, 240,  94],
         [210, 250, 162],
         [157, 246,  97],
         [97,  246, 127],
         [35,  202, 157],
         [122, 238, 246],
         [4,   173, 185],
         [25,  154, 253],
         [11,   99, 250],
         [53,   72, 144]]

  CLR = np.array(CLR)/255.
  CLR = np.flip(CLR, axis=0)
  CMP = mclrs.create_colormap(CLR, 200)

  return CMP


def plot_polar_2D(LON, LAT, A2D, region='Arctic', nfg=1, \
           rmin=-1.e30, rmax=-1.e30, cmpice=[], stl='CICE6 output', \
           cntr1=[], clr1=[], cntr2=[], clr2=[]):
  """
    Plot CICE6 field in polar projection
    alows plotting contours for 2 types of data, e.g. pos/neg SSH
    specify cntr1, cntr2 and colors
  """
  from mpl_toolkits.basemap import Basemap, cm
  import matplotlib.colors as colors
  import matplotlib.mlab as mlab
  from matplotlib.colors import ListedColormap

  LONA = sub_region2D(LON, region=region)
  LATA = sub_region2D(LAT, region=region)
  A2D  = sub_region2D(A2D, region=region)

# Check if colormap object is provided:
  try:
    Nclr = cmpice.colors
  except:
    cmpice = colormap_conc()  
    cmpice.set_under(color=[1, 1, 1])

  cmpice.set_bad(color=[0.3, 0.3, 0.3])

# Projection for N. Hemisphere:
  if region == 'Arctic':
    m = Basemap(projection='npstere',boundinglat=50,lon_0=0,resolution='l')
  else:
    m = Basemap(projection='spstere',boundinglat=-50,lon_0=0,resolution='l')


  plt.ion()
  fig1 = plt.figure(nfg,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.1, 0.15, 0.75, 0.75])

  # draw parallels.
  if region == 'Arctic':
    parallels = np.arange(0.,90,10.)
  else:
    parallels = np.arange(-90,-10,10.)

  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
  # draw meridians
  meridians = np.arange(-360,360.,45.)
  m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)


  ny = A2D.shape[0]; nx = A2D.shape[1]

  lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly spaced grid.
  x, y = m(lons, lats) # compute map proj coordinates.
  xh, yh = m(LONA,LATA)  # hycom coordinates on the projections
  im1 = ax1.pcolormesh(xh, yh, A2D, shading='flat',
                       cmap=cmpice)

  if rmin > -1.e30 and rmax > -1.e30:
    im1.set_clim(rmin,rmax)

  if len(cntr1) > 0:
    if len(clr1) == 0:
      clr1=[(0,0,0)]

    ax1.contour(xh, yh, A2D, cntr1, linestyles='solid',
                colors=clr1, linewidths=1)

  if len(cntr2) > 0:
    if len(clr2) == 0:
      clr2=[(0.5,0.5,0.5)]

    ax1.contour(xh, yh, A2D, cntr2, linestyles='solid',
                colors=clr2, linewidths=1)

  ax1.set_title(stl) 
 
  ax2 = fig1.add_axes([ax1.get_position().x0, ax1.get_position().y0-0.085,
                       ax1.get_position().width,0.02])
  clb = plt.colorbar(im1, cax=ax2, orientation='horizontal')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  return

def sice_lr_cice4(klr, Ni, aice_lr, Smax=3.2, a=0.407, b=0.573):
  """
    Compute the vertical S profile for CICE4 
    following CICE4 thermodynamics 
    Each column is divided into N ice lyaers dlth = thickness hice/Ni
    S profile is prescribed and is unchanging in time, the snow is assumed
    to be fresh and the midpoint Sik is computed from a formula
    S goes from S=0 at the surf to Smax

    The algorithm yields 2D fields of S at mid-point of layer = klr for the total
    number of ice layers = Ni
  """
  if klr > Ni:
    raise Exception ('klr {0} cannot be > Ni {1}'.format(klr,Ni))

  sice_lr = aice_lr.copy()*0.0
  z       = (klr-0.5)/float(Ni)
  sice    = 0.5*Smax*(1.-np.cos(np.pi*z**(a/(z+b))))
  sice_lr = np.where(aice_lr < 1.e-10, 0.0, sice)

  return sice_lr
 
 
def remap_enthalpy_bins(qicen4, nilyrs4, nilyrs6, eps_dq=1.e-5):
  """
    For remapping snow or ice enthalpies from N layers
    to K (K>N) layers
    Note that enthalpies have to be defined uniqnuely for this step
    i.e. J/m3 (CICE6) enthalpy per vol
    In CICE4 --> CICE6,
    qicen4 already converted CICE4 esnon (J/m2) to qsnon (J/m3)
    then remap onto K snow/ice layers

    eps_dq - tolerance error for interpolated enthalpie
             small error (<eps_dq) is distributed
             across the cice6 layers

    Remapping preserves total enthalpy integrated over all layers
    Remapping done by binning
   for mapping CICE4 ---> CICE6

     cice4                    cice 6

    k=5 ------  z=1      k=8 --------  z=1
                              lr=7
                         k=7 -------  z=0.857
    k=4 -----  z=0.75         lr=6
                         k=6 -------  z=0.714
                              lr=5
                         k=5 -------  z=0.714
    k=3 -----  z=0.5


    qicen4 has shape (nilyrs4, ncat, ny, nx)
  """
  zz4    = np.zeros((nilyrs4+1)) # interface coordinates
  zz6    = np.zeros((nilyrs6+1))
  dh4    = 1./float(nilyrs4)
  dh6    = 1./float(nilyrs6)
  ny     = qicen4.shape[2]
  nx     = qicen4.shape[3]
  ncat   = qicen4.shape[1]
  qicen6 = np.zeros((nilyrs6, ncat, ny, nx))

  for ik in range(1,nilyrs4+1):
    zz4[ik] = zz4[ik-1] + dh4

  for ik in range(1,nilyrs6+1):
    zz6[ik] = zz6[ik-1] + dh6

  for n in range(ncat):
    eta_tot = 0.       
    qtot = qicen4[0,0,:,:]*0.0  # check total enthalpy must be conserved
    for k in range(nilyrs6):
      zbtm = zz6[k]
      ztop = zz6[k+1]
  
# Find cice4 ice layer interfaces for cice6 layer k
# cice6 bottom interface
      i4btm_btm = np.where(zz4 <= zbtm)[0][-1] # btm cice4 for cice6 btm
      i4btm_top = i4btm_btm + 1       # top cice4 for cice6 btm
      z4btm_btm = zz4[i4btm_btm]      # depth of cice4 btm intrf
      z4btm_top = zz4[i4btm_top]      # depth of cice4 srf intrf
    
# cice6 top interface:
      i4top_btm = np.where(zz4 <= ztop)[0][-1] # btm cice4 
      i4top_top = i4top_btm + 1       # top cice4 for cice6 top
      z4top_btm = zz4[i4top_btm]
      z4top_top = zz4[i4top_top]
    
# eta1 = thickn from zbtm cice6 to cice4 top interface
# if cice6 lr is inside 1 cice4 lr eta1 = dh6 (cice6 lr thickness)
      eta1   = min([z4btm_top,ztop]) - zbtm
# eta2 = thickn from ztop cice6 to zbtm cice4
# Check if cice6 ice layer is within cice 4 layer
      if i4btm_btm == i4top_btm and i4btm_top == i4top_top:
        eta2   = 0.
      else:
        eta2 = ztop - max([z4top_btm,zbtm])

# Compute enthalpie in cice6 layer:
      k4_btm = i4btm_btm       # ice layer cice4 where cice6 btm
      k4_top = i4top_btm       # ice layer cice4 where cice6 top
      qn6 = (qicen4[k4_btm, n, :, :]*eta1 + \
            qicen4[k4_top, n, :, :]*eta2) / dh6

      qicen6[k,n,:,:] = qn6
      qtot    = qtot + qn6*dh6
      eta_tot = eta_tot + eta1 + eta2
      print('cat={0}, lr={1}, total_fraction={2}'.format(n+1,k+1,eta_tot))

    q4tot = np.sum(np.squeeze(qicen4[:,n,:,:]), axis=0) * dh4
    dq    = np.abs(q4tot-qtot)
    print('Cat={0} Max abs dlt enthalipes in cice4 and cice6: {1}'.\
           format(n,np.max(dq)))
    if np.max(dq) > eps_dq:
      raise Exception ('Large Error in interpolated enthalpy to cice6 ')

  i_chck = False
  if i_chck:
    ii = 1589
    jj = 32

    zi4  = np.linspace(dh4/2., 1.-dh4/2., num=4)
    zi6  = np.linspace(dh6/2., 1.-dh6/2., num=7)
    zzi4 = np.linspace(0., 1., num=5)
    zzi6 = np.linspace(0., 1., num=8)

    q4  = qicen4[:,n,jj,ii]
    q6  = qicen6[:,n,jj,ii]
    dqij = np.sum(q4*dh4) - np.sum(q6*dh6)
    q6c = q6 + dqij

    ll1 = np.min(q4)
    ll2 = np.max(q4)
    
    fig1 = plt.figure(1,figsize=(9,8))
    fig1.clf()
# Enthalpie of cat=n by layers J/m3 cice4
    for ik in range(nilyrs4):
      qq = q4[ik]
      s1 = zzi4[ik]
      s2 = zzi4[ik+1]
      plt.plot([qq, qq],[s1, s2],'-', color=[0.8,0.2,0])

# Enthalpie of cat=n by layers J/m3 cice6
    for ik in range(nilyrs6):
      qq = q6[ik]
      s1 = zzi6[ik]
      s2 = zzi6[ik+1]
      plt.plot([qq, qq],[s1, s2],'-', color=[0.,0.2,0.9])

# Plot sea ice layers for cice4
    for ik in range(nilyrs4+1):
      plt.plot([ll1,ll2],[zzi4[ik], zzi4[ik]],'--', color=[0.5,0.5,0.5])
  
    plt.title('Ice enthalpy, J/m3, CICE4(r) and CICE6(b), cat={0}'.format(n))
    btx = 'mod_cice6_utils.py'
    bottom_text(btx) 


  return qicen6

