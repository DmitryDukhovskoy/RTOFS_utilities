"""
  Follow CICE icepack icepack_itd.F90
  zap_snow_temperature subroutine
  to compute Tmax and zTsn (snow T) from snow enthalpy
  Note that saved enthalpy is layer internal energy in J/m3
  the subroutine uses also  J/m3
 
  See comments in icepack_therm_vertical.F90
      !-----------------------------------------------------------------
      ! Tmax based on the idea that dT ~ dq / (rhos*cp_ice)
      !                             dq ~ q dv / v
      !                             dv ~ puny = eps11
      ! where 'd' denotes an error due to roundoff.
      !-----------------------------------------------------------------
  
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


pthrst      = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6_003/RESTART/'
#cicerst6    = 'iced.{0}-{1:02d}-{2:02d}-{3:05d}.nc'.\
#               format(YR,MM,MD,HR)
#cicerst6    = 'iced.{0}-{1:02d}-{2:02d}-gofs3.2.nc'.\
#               format(YR,MM,MD)
cicerst6    = 'iced.2020-01-01-00000.nc'
#cicerst6    = 'iced.2022-01-24-00000.nc'
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

# Read vsnon:
for kcat in range(ncat):
  vsnon  = ncdata['vsnon'][kcat,:,:].data
  fld    = 'qsno{0:03d}'.format(1)
  qsnon  = ncdata[fld][kcat,:,:].data
  aicen  = ncdata['aicen'][kcat,:,:].data
#  snow thickness
  aicen  = np.where(aicen < puny, np.nan, aicen)
  vsnon  = np.where(np.isnan(aicen), np.nan, vsnon)
  hsn    = vsnon / aicen
#  hsn    = np.where(np.isnan(aicen), 1.e-11, hsn)
  
#
# Convert enthalpy to layer internal energy
# to make it similar to icepack_itd code
  nslyr = 1.
#  zqsn  = qsnon*vsnon/nslyr
  zqsn  = qsnon
  zqsn  = np.where(zqsn >= -0.0, np.nan, zqsn)
  zqsn  = np.where(np.isnan(aicen), np.nan, zqsn)
#  zqsn  = np.where(vsnon < 1.e-32, np.nan, zqsn)
  zqsn  = np.where(hsn < hs_min, -rhos*Lfresh, zqsn)

  zTsn = (Lfresh + zqsn/rhos)/cp_ice 
  Tmax = -zqsn*puny*nslyr / (rhos*cp_ice*vsnon)

  J0, I0 = np.where(zTsn > 1.e-11)
  j0,i0 = np.where(zqsn == np.nanmax(zqsn))
  j0=j0[0]
  i0=i0[0]
  zq0   = zqsn[j0,i0]
  vs0   = vsnon[j0,i0]
  hs0   = hsn[j0,i0]
  ai0   = aicen[j0,i0]
  zT0   = (Lfresh + zq0/rhos)/cp_ice

  print('cat={0}, max zqsn pnt: zq0={1}, \n vs0={2}, hs0={3}'.\
         format(kcat+1, zq0, vs0, hs0))
  print('  ai0={0}, zT0={1}\n'.format(ai0, zT0))




