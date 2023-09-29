"""
  Plot CICE6 3D ice state variables aicen, vicen, vsnon [ncat, ny, nx]
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

nx       = cice6.nx
ny       = cice6.ny
ncat     = cice6.ncat
ulati6   = mc6util.read_ncfile(fgrdin, 'ulat')
uloni6   = mc6util.read_ncfile(fgrdin, 'ulon')
LAT, LON = mc6util.grid_rad2dgr(ulati6, uloni6)


ncdata = ncFile(fl_restart6,'r')

#importlib.reload(mc6util)
btx = 'plot_restart_icestate.py'
rmin   = 0.  # min conc
rmax   = 1.  # max conc
cmpice = mc6util.colormap_conc()

fld   = 'aicen'
# To plot values aggregated over all categories
# make kcat_plot > ncat:
kcat_plot = 100
if kcat_plot < ncat:
  kcat     = kcat_plot - 1
  str_kcat = '{0}'.format(kcat_plot)
  A2D      = ncdata[fld][kcat,:,:].data
  A2D      = np.where(Lmsk==0, np.nan, A2D)
else:
# Aggregated:
  print('aice aggregated over all categories')
  A2D  = np.zeros((ny,nx))
  str_kcat = 'ALL cat'
  for kcat in range(ncat):
    dmm = ncdata[fld][kcat,:,:].data
    dmm = np.where(Lmsk==0, np.nan, dmm)
    A2D = A2D + dmm

rgn    = 'Arctic'

stl = 'CICE6, {0}/{1}/{2}, {3}, icat={4} \n {5}'.\
        format(YR, MM, MD, fld, str_kcat, fl_restart6)
mc6util.plot_polar_2D(LON, LAT, A2D, region=rgn,  \
                rmin=rmin, rmax=rmax, cmpice=cmpice, stl=stl)
bottom_text(btx, pos=[0.05, 0.02])





