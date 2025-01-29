# Plot grid resolution for NEP domain
#
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
#import struct
import datetime
#import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
import yaml
from netCDF4 import Dataset as ncFile

#PPTHN = '/home/Dmitry.Dukhovskoy/python'
PPTHN = []
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mmom6
import mod_misc1 as mmsc1
#import mod_valid_utils as mvutil
importlib.reload(mcmp)


nrun = "MOM6_NEP"
expt = "test"

with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)
#
# Read MOM6 NEP domain nx=342, ny=816
pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom) 
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom) 
LONM, LATM  = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
HHM         = mmom6.read_mom6depth(dftopo_mom) 
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]

# Read grid resolution:
# dX, dY are on MOM "supergrid" - half grid points
DX, DY = mhycom.dx_dy(LONM, LATM)

# Effective resolution as norm-2
RS2 = np.sqrt(DX**2 + DY**2)*1.e-3  # km
mm, nn = RS2.shape

# Resolution as a max of DX, DY
D = DX-DY
RS1 = np.where(D>=0, DX, DY)*1.e-3 # km

RS1 = np.where(HHM >= 0., np.nan, RS1)
RS2 = np.where(HHM >= 0., np.nan, RS2)


# Set up orthographic projection
from mpl_toolkits.basemap import Basemap, cm
lon0 = 220.
lat0 = 50.
res  = 'l'
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
xR, yR = m(LONM,LATM)

PMsk = ( (xR > 1e20) | (yR > 1e20) )
AA = RS1.copy()
AA = np.where(HHM >= 0, np.nan, AA)
AA[PMsk] = np.nan
xR[PMsk]   = 1.e30
yR[PMsk]   = 1.e30

rmin = 8.
rmax = 11.

plt.ion()

clrmp = mcmp.colormap_warm(nclrs=200)
clrmp.set_bad(color=[0.2,0.2,0.2])

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])
m.drawcoastlines()
im1 = m.pcolormesh(xR, yR, AA, cmap=clrmp, vmin=rmin, vmax=rmax)

m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

stl = 'MOM6-SIS2 grid resolution max(dx,dy), km,  NEP region'
ax1.set_title(stl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


btx = 'plot_resolution_ortho.py'
bottom_text(btx)





