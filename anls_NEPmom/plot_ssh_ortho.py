# Plot ssh fields on ortho
# Plot regional Rossby / max(dx, dy) of the grid cells
# gives # of grid cells per R
#
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#import torch
import sys
import pdb
import netCDF4
import importlib
from netCDF4 import Dataset as ncFile
import timeit
import pickle
import yaml


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

YR      = 1993
pfld    = 'ssh'
expt    = 'NEP_BGCphys'
pthout  = f'/work/Dmitry.Dukhovskoy/run_output/NEP_BGCphys/{YR}'

btx = 'plot_ssh_ortho.py'

with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom)
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom)
LONM, LATM  = mom6util.read_mom6grid(dfgrid_mom)
HHM         = mom6util.read_mom6depth(dftopo_mom)
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]

jday   =  10
dnmb   = mtime.jday2dnmb(YR, jday)
dv     = mtime.datevec(dnmb)
ds     = mtime.datestr(dnmb)
MM     = dv[1]
DM     = dv[2]

flmom  = f'ocean_{YR}_{jday:03d}_12.nc'
dflmom = os.path.join(pthout,flmom)
nc     = ncFile(dflmom,'r')
A2d    = nc.variables[pfld][:].data.squeeze()
A2d    = np.where(A2d > 1.e10, np.nan, A2d)
A2d    = A2d - np.nanmean(A2d)
nx     = A2d.shape[1]
ny     = A2d.shape[0]


ctitle = f'{expt} MOM6 {pfld} {YR}/{MM:02d}/{DM:02d}'
cmpr = mutil.colormap_ssh(nclrs=200)
rmin = -0.5
rmax = 0.5
cmpr.set_bad(color=[0.1, 0.1, 0.1])


# Set up orthographic projection
from mpl_toolkits.basemap import Basemap, cm

# Add extra row/col for plotting with pcolormesh
lonw = LONM.copy()
latw = LATM.copy()
lonw = np.insert(lonw, -1, lonw[:,-1]+0.01, axis=1)
lonw = np.insert(lonw, -1, lonw[-1,:]+0.01, axis=0)
latw = np.insert(latw, -1, latw[-1,:]+0.01, axis=0)
latw = np.insert(latw, -1, latw[:,-1]+0.01, axis=1)

lon0 = 220.
lat0 = 50.
res  = 'l'
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
xR, yR = m(lonw, latw)
PMsk = ( (xR > 1e20) | (yR > 1e20) )
data = A2d.copy()
data = np.insert(data, -1, data[:,-1], axis=1)
data = np.insert(data, -1, data[-1,:], axis=0)
data[PMsk] = np.nan
xR[PMsk]   = 1.e30
yR[PMsk]   = 1.e30

data = data[0:ny, 0:nx]


#================
#    Plotting

plt.ion()

fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
plt.clf()

ax1 = plt.axes([0.04, 0.04, 0.8, 0.8])
m.drawcoastlines()
im1 = m.pcolormesh(xR, yR, data, shading='flat', cmap=cmpr,\
                   vmin=rmin, vmax=rmax)

m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

ax1.set_title(ctitle)

ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=10)

#fig1.colorbar(im,ax=ax1,orientation='horizontal')

bottom_text(btx, pos=[0.02, 0.02])



