"""
  Print min/max values from the restart
  can be compared with the info in ice_diag.d
  from CICE6 restart
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
from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime


ncat = 5
nilyrs = 7

# CICE4 restart date:
YR  = 2000
MM  = 1
MD  = 27
HR  = 0


pthrst      = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6/RESTART/'
#cicerst6    = 'iced.{0}-{1:02d}-{2:02d}-{3:05d}.nc'.\
#               format(YR,MM,MD,HR)
cicerst6    = 'iced.{0}-{1:02d}-{2:02d}-gofs3.2.nc'.\
               format(YR,MM,MD)
fl_restart6 = pthrst + cicerst6

puny      = 1.e-11
c0        = 0.0
c1        = 1.0
c2        = 2.0
p5        = 0.5
Lsub      = 2.835e6    # latent heat sublimation fw (J/kg)
Lvap      = 2.501e6    # latent heat vaporization fw (J/kg)
Lfresh    = Lsub - Lvap # latent heat of melting of fresh ice (J/kg)
cp_ice    = 2106.       # specific heat of fresh ice (J/ kg/K)
rhos      = 330.        # density of snow (kg/m3)
hs_min    = 1.e-4       # min snow thickness for computing Tsno (m)
nsal      = 0.407
msal      = 0.573
min_salin = 0.1      # threshold for brine pocket treatment
saltmax   = 3.2        # max S at ice base


ncdata = ncFile(fl_restart6,'r')

def print_minmax(sfld,A):
  print('   {2} min/max:  {0}/{1}'.format(np.nanmin(A),np.nanmax(A),sfld))

  return

def print_minmaxsum(sfld,A):
  print(' {3} min, max, sum =  {0}, {1}, {2}'.format(np.nanmin(A), np.nanmax(A), np.nansum(A), sfld))

  return

def print_2D(fld, ltrunc_less0 = False, ltrunc_grt0 = False):
  AA    = ncdata[fld][:,:]
  print_minmaxsum(fld, AA)

  return

def print_3D(fld, ncat=5, ltrunc_less0 = False, ltrunc_grt0 = False):
  """
    Read 3D fields from restart
    and print min/max
  """
  for kcat in range(ncat):
    lsnow = 0
    AA    = ncdata[fld][kcat,:,:]

# Keep only <0 values (e.g., enthalpy):
    if ltrunc_less0:
      AA = np.where(AA >= 0., np.nan, AA)
# Keep only >0 values (e.g., aice):
    if ltrunc_grt0:
      AA = np.where(AA <= 0., np.nan, AA)

    print_minmaxsum(fld+'_{0}'.format(kcat+1), AA)

  return
    
print_3D('aicen')
print_3D('vicen')
print_3D('vsnon')
print_3D('Tsfcn', ltrunc_less0 = True)
for ilyr in range(nilyrs):
  fld='sice{0:03d}'.format(ilyr+1)
  print_3D(fld) 
for ilyr in range(nilyrs):
  fld='qice{0:03d}'.format(ilyr+1)
  print_3D(fld, ltrunc_less0 = True)
fld='qsno{0:03d}'.format(1)
print_3D(fld, ltrunc_less0 = True)

print_2D('uvel')
print_2D('vvel')
print_2D('coszen')
print_2D('scale_factor')
print_2D('swvdr')
print_2D('swvdf')
print_2D('swidr')
print_2D('swidf')
print_2D('strocnxT')
print_2D('strocnyT')






