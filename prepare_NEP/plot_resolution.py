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
import mod_mom6 as mom6util
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
LONM, LATM  = mom6util.read_mom6grid(dfgrid_mom)
HHM         = mom6util.read_mom6depth(dftopo_mom) 
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]

# Read grid resolution:
# dX, dY are on MOM "supergrid" - half grid points
nc  = ncFile(dfgrid_mom,'r')
dX  = nc.variables['dx'][:].data
dY  = nc.variables['dy'][:].data
nyp = dX.shape[0]
nx  = dX.shape[1]
ny  = dY.shape[0]
nxp = dY.shape[1]

DX = dX[0:nyp-1:2,0:nx:2] + dX[1:nyp:2,1:nx:2]
DY = dY[0:ny:2, 0:nxp-1:2] + dY[1:ny:2, 0:nxp-1:2] 
#DX = dX[0:nyp-1:2,0:nx:2]
#DY = dY[0:ny:2, 0:nxp-1:2]
RS = np.sqrt(DX**2 + DY**2)*1.e-3  # km
mm = RS.shape[0]
nn = RS.shape[1]

# Check distance:
DSouth = mmsc1.dist_sphcrd(LATM[0,0], LONM[0,0], LATM[0,-1], LONM[0,-1])*1e-3
DLS    = np.nansum(DX[0,0:])*1e-3
DNorth = mmsc1.dist_sphcrd(LATM[-1,0], LONM[-1,0], LATM[-1,-1], LONM[-1,-1])*1e-3
DLN    = np.nansum(DX[-1,0:])*1e-3


RS = np.where(HHM >= 0., np.nan, RS)


plt.ion()

clrmp = mcmp.colormap_temp(nclrs=200)
clrmp.set_bad(color=[0.2,0.2,0.2])

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])
rmin = 5.
rmax = 15.
im1 = ax1.pcolormesh(RS, \
                 cmap=clrmp,\
                 vmin=rmin, \
                 vmax=rmax)

ax1.axis('scaled')
ax1.set_xlim([0, nn-1])
ax1.set_ylim([0, mm-1])

LON1 = LONM.copy()
LON1 = np.where(LON1<0., LONM+360., LONM)
clrg = [(0.95,0.95,0.95)]
ax1.contour(LON1,list(range(100,320,10)),
          colors=clrg,
          linestyles='solid', 
          linewidths=1.0)
ax1.contour(LATM,list(range(0,89,10)),
          colors=clrg,
          linestyles='solid',
          linewidths=1.0)


stl = 'MOM6-SIS2 grid resolution ||dl||$_2$, km,  NEP region'
ax1.set_title(stl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.0f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


btx = 'plot_resolution.py'
bottom_text(btx)





