# Plot ice layer enthalpy converted to layer internal energy (J/m3)
# to compare with CICE6
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
YR  = 2020
MM  = 1
MD  = 1
HR  = 9


pthrst   = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6/RESTART/'
pthglb   = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/GLBb0.08_expt93.0/'
cicerst4 = 'cice.restart.{0}{1:02d}{2:02d}{3:02d}'.format(YR,MM,MD,HR)
fl_restart4 = pthglb + cicerst4

# Grid CICE4 - unformatted binary file
pthgrd4 = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_paraCd/fix/'
grdfl4  = 'regional.cice.r'
fgrdin4 = pthgrd4 + grdfl4


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
spval     = 1.e30


import mod_cice6_utils as mc6util
importlib.reload(mc6util)

# CICE4 params and restart fields
# Note in GOFS3.1 CICE has 1 row less than HYCOM
# ny = 3297 (and it is 3298 in HYCOM)
#cice4 = mc6util.param_cice4()
cice4 = mc6util.cice4(nx=4500, ny=3297)


aicen = mc6util.read_cice4_restart(fl_restart4, cice4, 'aicen')
aicen = np.where(aicen > 0.5*spval, 0., aicen)

vicen = mc6util.read_cice4_restart(fl_restart4, cice4, 'vicen')
vicen = np.where(vicen > 0.5*spval, 0., vicen)

fld   = 'eicen'
eicen = mc6util.read_cice4_restart(fl_restart4, cice4, fld)
eicen = np.where(eicen > 0.5*spval, 0., eicen)


# Convert layer enthalpy (J/m2) --> layer internal energy (J/m3):
nilyr = cice4.nilyr
nslyr = cice4.nslyr
ncat  = cice4.ncat
nx    = cice4.nx
ny    = cice4.ny
ilyr1 = np.zeros((ncat), dtype='int')
slyr1 = np.zeros((ncat), dtype='int')
ilyr1[0] = 0   # 1st index, python
slyr1[0] = 0
for n in range(1,ncat):
  ilyr1[n] = ilyr1[n-1] + nilyr
  slyr1[n] = slyr1[n-1] + nslyr


# Convert enthalpy CICE4 J/m2 ---> CICE6 J/m3
# on the same vertical ice layers 
qicen = np.zeros((nilyr,ncat,ny,nx))
for n in range(ncat):                 # ice categories
  for k in range(nilyr):
    ai_cat = aicen[n,:,:]
    iicen  = ilyr1[n]+k
    print('eicen index: iicen = {0}'.format(iicen))
    vin    = vicen[n,:,:]               # ice volume per m2 in cat=n
    vin    = np.where(vin < 1.e-30, 1.e-30, vin)
    qin    = eicen[iicen,:,:]*float(nilyr)/vin  # J/m2 --> J/m3
    qin    = np.where(ai_cat <= puny, 0.0, qin)
    qicen[k,n,:,:] = qin


# Lon/lat in CICE4:
ulati4   = mc6util.read_cice4_grid(fgrdin4, 'ulati', IDM=nx, JDM=ny)
uloni4   = mc6util.read_cice4_grid(fgrdin4, 'uloni', IDM=nx, JDM=ny)
LAT, LON = mc6util.grid_rad2dgr(ulati4, uloni4)

pthdpth = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6/INPUT/'
dpthfl  = 'depth_GLBb0.08_09m11ob2_mom6.nc'
fdpthin = pthdpth + dpthfl

Lmsk  = mc6util.read_ncfile(fdpthin,'wet')


btx = 'plotcice4_restart_qice.py'
rmin   = -350.  # MJ
rmax   = -200.  # max enthalpy to give T snow <= 0C
cmpice = mc6util.colormap_enthlp()

kcat  = 0
ilyr  = 3
cff   = 1.e-6
A2D   = np.squeeze(qicen[ilyr, kcat, :, :])*cff
A2D   = np.where(Lmsk==0, np.nan, A2D)
rgn   = 'Arctic'

stl = 'CICE4 GOFS3.2-93.2, {0}/{1}/{2}, {3} J/m^3*{4}, icat={5}, lyr={6}\n {7}'.\
        format(YR,MM,MD,'qicen',cff,kcat+1,ilyr+1,fl_restart4)
mc6util.plot_polar_2D(LON, LAT, A2D, region=rgn,  \
                rmin=rmin, rmax=rmax, cmpice=cmpice, stl=stl)
bottom_text(btx, pos=[0.05, 0.02])





