# Plot scale factor 
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


scfct = mc6util.read_cice4_restart(fl_restart4, cice4, 'scale_factor')
scfct = np.where(scfct > 0.5*spval, 0., scfct)

# Lon/lat in CICE4:
nx       = cice4.nx
ny       = cice4.ny
ncat     = cice4.ncat
ulati4   = mc6util.read_cice4_grid(fgrdin4, 'ulati', IDM=nx, JDM=ny)
uloni4   = mc6util.read_cice4_grid(fgrdin4, 'uloni', IDM=nx, JDM=ny)
LAT, LON = mc6util.grid_rad2dgr(ulati4, uloni4)

pthdpth = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6/INPUT/'
dpthfl  = 'depth_GLBb0.08_09m11ob2_mom6.nc'
fdpthin = pthdpth + dpthfl

Lmsk  = mc6util.read_ncfile(fdpthin,'wet')


btx = 'plotcice4_restart_scalefct.py'
rmin   = 0.  #
rmax   = 10.  # 
cmpice = mc6util.colormap_conc()

A2D = scfct

rgn   = 'Arctic'

stl = 'CICE4 GOFS3.2-93.2, {0}/{1}/{2}, {3}\n {4}'.\
        format(YR, MM, MD, 'scale_factor', fl_restart4)
mc6util.plot_polar_2D(LON, LAT, A2D, region=rgn,  \
                rmin=rmin, rmax=rmax, cmpice=cmpice, stl=stl)
bottom_text(btx, pos=[0.05, 0.02])


# Plot global:
cmpice = mc6util.colormap_conc()
rmin  = 0.
rmax  = 200.
cmpice.set_bad(color=[0.3, 0.3, 0.3])

plt.ion()
fig1 = plt.figure(2,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.15, 0.75, 0.75])
im1 = ax1.pcolormesh(A2D, cmap=cmpice)
im1.set_clim(rmin,rmax)

ax1.contour(Lmsk,[0.01], colors=[(0,0,0)], linewidths=1)
ax1.set_title(stl)

ax2 = fig1.add_axes([ax1.get_position().x0, ax1.get_position().y0-0.085,
                     ax1.get_position().width,0.02])
clb = plt.colorbar(im1, cax=ax2, orientation='horizontal')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.tick_params(direction='in', length=12)

bottom_text(btx, pos=[0.05, 0.02])



