"""
  Plot MOM-NEP output T, S, or ssh 
  orthonormal projection
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
#import struct
import datetime
import pickle
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
#import mod_valid_utils as mvutil
importlib.reload(mcmp)

nrun   = "MOM6_NEP"
#expt   = "NEP_BGCphys_GOFS"
expt   = "NEP_BGCphys"
dnmb0  = mtime.datenum([1993,5,1])
varplt = 'salt' # potT salt ssh 
kplt   = 0


with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

# Read MOM6 NEP domain nx=342, ny=816
pthgrid_mom = dct["MOM6_NEP"][expt]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"][expt]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"][expt]["fgrid"]
pthoutp     = dct["MOM6_NEP"][expt]["pthoutp"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom)
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom)
LONM, LATM  = mom6util.read_mom6grid(dfgrid_mom)
HHM         = mom6util.read_mom6depth(dftopo_mom)
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]

DV      = mtime.datevec(dnmb0)
YR      = DV[0]
MM      = DV[1]
DD      = DV[2]
_, jday = mtime.dnmb2jday(dnmb0)
jday    = int(jday)
flmom   = f'ocean_{YR}_{jday:03d}_12.nc'
dflmom  = os.path.join(pthoutp,f"{YR}",f"{MM:02d}",flmom)
nc      = ncFile(dflmom,'r')
f_ssh   = False
f_temp  = False
f_salt  = False
if varplt == "ssh":
  f_ssh = True
elif varplt == "potT":
  f_temp = True
elif varplt == "salt":
  f_salt = True

if f_ssh:
  A2d   = nc.variables[varplt][:].data.squeeze()
else:
  A2d   = nc.variables[varplt][:,kplt,:,:].data.squeeze()

A2d     = np.where(A2d > 1.e10, np.nan, A2d)
if f_ssh:
  A2d   = A2d - np.nanmean(A2d)
  cmpr  = mutil.colormap_ssh(nclrs=200)
  rmin = -0.5
  rmax = 0.5
  cntr1 = [x/10 for x in range(-10, 10, 1)]
elif f_salt:
  cmpr  = mcmp.colormap_haline()
  rmin  = 31.5
  rmax  = 35.5
  if kplt>45:
    rmin = 34.
    rmax = 34.9
  cntr1 = [x/10 for x in range(348, 358, 1)]
elif f_temp:
  cmpr = mcmp.colormap_temp_coldhot()
  rmin  = -2.
  rmax  = 25.
  if kplt>15:
    rmin = -1.
    rmax = 7.
  cntr1 = [x/10 for x in range(-10, 250, 10)]

nx     = A2d.shape[1]
ny     = A2d.shape[0]

ctitle = f'MOM6 {expt} {varplt} Lr={kplt+1} {YR}/{MM:02d}/{DD:02d}'
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
AA = A2d.copy()
AA = np.insert(AA, -1, AA[:,-1], axis=1)
AA = np.insert(AA, -1, AA[-1,:], axis=0)
AA[PMsk] = np.nan
xR[PMsk]   = 1.e30
yR[PMsk]   = 1.e30

AA = AA[0:ny, 0:nx]


#================
#    Plotting

plt.ion()

fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
plt.clf()

ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
im1 = m.pcolormesh(xR, yR, AA, shading='flat', cmap=cmpr,\
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
btx = 'plot_momTS_ortho.py'
bottom_text(btx, pos=[0.02, 0.02])





