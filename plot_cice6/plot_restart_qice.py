"""
  Plot CICE6 3D ice enthalpy qice[ncat, ny, nx]
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

# CICE4 restart date:
YR  = 2000
MM  = 1
MD  = 27
HR  = 0


pthrst      = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6/RESTART/'
#cicerst6    = 'iced.{0}-{1:02d}-{2:02d}-{3:05d}.nc'.\
#               format(YR,MM,MD,HR)
# Generated restart from GOFS3.2 CICE4
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


import mod_cice6_utils as mc6util

# Grid, use CICE6 grid - should be similar to CICE4
pthgrd  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6/'
grdfl   = 'grid_cice_NEMS_mx008.nc'
fgrdin  = pthgrd + grdfl

pthdpth = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6/INPUT/'
dpthfl  = 'depth_GLBb0.08_09m11ob2_mom6.nc'
fdpthin = pthdpth + dpthfl

Lmsk  = mc6util.read_ncfile(fdpthin,'wet')
cice6 = mc6util.cice6()

ulati6 = mc6util.read_ncfile(fgrdin, 'ulat')
uloni6 = mc6util.read_ncfile(fgrdin, 'ulon')
LAT, LON = mc6util.grid_rad2dgr(ulati6, uloni6)



ncdata = ncFile(fl_restart6,'r')

#importlib.reload(mc6util)
btx = 'plot_restart_qice.py'
rmin   = -350.  # MJ
rmax   = -200.  # max enthalpy to give T snow <= 0C
cmpice = mc6util.colormap_enthlp()

kcat  = 0
ilyr  = 6
cff   = 1.e-6
fld   = 'qice{0:03d}'.format(ilyr+1)
A2D    = ncdata[fld][kcat,:,:].data*cff
A2D    = np.where(Lmsk==0, np.nan, A2D)
rgn    = 'Arctic'

stl = 'CICE6, {0}/{1}/{2}, {3} J/m^3*{4}, icat={5}, lyr={6}\n {7}'.\
        format(YR, MM, MD, fld, cff, kcat+1, ilyr+1, fl_restart6)
mc6util.plot_polar_2D(LON, LAT, A2D, region=rgn,  \
                rmin=rmin, rmax=rmax, cmpice=cmpice, stl=stl)
bottom_text(btx, pos=[0.05, 0.02])


for kcat in range(5):
  lice  = 0
  fld   = 'qice{0:03d}'.format(lice+1)
  AA    = ncdata[fld][kcat,:,:].data


  # Check max snow T from max snow enthalpy:
  AA = np.where(AA >= 0., np.nan, AA)

  J,I = np.where(AA < 0.)
  qsn = np.max(AA[J,I])
  zTsn = (Lfresh + qsn/rhos)/cp_ice
  j0, i0 = np.where(AA == qsn)
  j0 = j0[0]
  i0 = i0[0]

  print('{0}, i={5}, j={6}, Cat: {1}, Lr: {2}, max enth={3}, Temp={4}'.\
        format(fld, kcat+1, lice+1, qsn, zTsn,i0,j0))





