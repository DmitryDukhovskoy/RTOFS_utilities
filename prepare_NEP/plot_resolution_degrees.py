"""
  Plot grid resolution for NEP domain
  in lon/lat degrees
"""
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
dlt_coord = 'lonlat'  # long or latit or lonlat = Euclidian 

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

dlt_lon = np.diff(LONM, axis=1)
dlt_lon = np.pad(dlt_lon, ((0, 0), (0, 1)), mode='edge')  # Pad last column
dlt_lat = np.diff(LATM, axis=0)
dlt_lat = np.pad(dlt_lat, ((0, 1), (0, 0)), mode='edge')  # Pad last row

dlt_dgr = np.sqrt(dlt_lon**2 + dlt_lat**2)  # 
mm,nn = LONM.shape

if dlt_coord == 'longit':
  RS = dlt_lon
elif dlt_coord == 'latit':
  RS = dlt_lat
elif dlt_coord == 'lonlat':
  RS = dlt_dgr

RS = np.where(HHM >= 0., np.nan, RS)


plt.ion()

clrmp = mcmp.colormap_temp(nclrs=200)
clrmp.set_bad(color=[0.2,0.2,0.2])

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
rmin = 0.
rmax = 0.22

ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
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


stl = f'MOM6-SIS2 grid dlt_{dlt_coord}  NEP region'
ax1.set_title(stl)



ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.3f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


btx = 'plot_resolution_degrees.py'
bottom_text(btx)





