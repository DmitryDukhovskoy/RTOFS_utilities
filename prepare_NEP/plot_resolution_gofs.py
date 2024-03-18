# Plot NEP domain on GOFS grid
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
#import mod_valid_utils as mvutil
importlib.reload(mcmp)


nrun = "GOFS3.1"
expt = "93.0"

with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthgrid = dct[nrun][expt]["pthgrid"]
ftopo   = dct[nrun][expt]["ftopo"]
fgrid   = dct[nrun][expt]["fgrid"]

LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

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
DX, DY      = mhycom.dx_dy(LON,LAT)

RS = np.sqrt(DX**2 + DY**2)*1.e-3  # km
mm = RS.shape[0]
nn = RS.shape[1]

RS = np.where(HH >= 0., np.nan, RS)

# Find boundaries of the NEP domain:
IBND = np.array([])
JBND = np.array([])
dstp = 10
print("Searching indices along western boundary ...")
for jj in range(0,jdm,dstp):
  if jj%100 == 0:
    print(f"  {jj/jdm*100:3.1f}% ...")
  x0 = LONM[jj,0]
  y0 = LATM[jj,0]
  ii0, jj0 = mutil.find_indx_lonlat(x0, y0, LON, LAT)

  IBND = np.append(IBND, ii0)
  JBND = np.append(JBND, jj0)

print("Searching indices along northern boundary ...")
for ii in range(0,idm,dstp):
  if ii%100 == 0:
    print(f"  {ii/idm*100:3.1f}% ...")
  x0 = LONM[jdm-1,ii]
  y0 = LATM[jdm-1,ii]
  ii0, jj0 = mutil.find_indx_lonlat(x0, y0, LON, LAT)

  IBND = np.append(IBND, ii0)
  JBND = np.append(JBND, jj0)

print("Searching indices along eastern boundary ...")
for jj in range(0,jdm,dstp):
  if jj%100 == 0:
    print(f"  {jj/jdm*100:3.1f}% ...")
  x0 = LONM[jj,idm-1]
  y0 = LATM[jj,idm-1]
  ii0, jj0 = mutil.find_indx_lonlat(x0, y0, LON, LAT)

  IBND = np.append(IBND, ii0)
  JBND = np.append(JBND, jj0)

print("Searching indices along southern boundary ...")
for ii in range(0,idm,dstp):
  if ii%100 == 0:
    print(f"  {ii/idm*100:3.1f}% ...")
  x0 = LONM[0,ii]
  y0 = LATM[0,ii]
  ii0, jj0 = mutil.find_indx_lonlat(x0, y0, LON, LAT)

  IBND = np.append(IBND, ii0)
  JBND = np.append(JBND, jj0)


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
ax1.set_xlim([1010, 2280])
ax1.set_ylim([1500, 3050])

ax1.plot(IBND, JBND, 'r.')

LON1 = LON.copy()
LON1 = np.where(LON1<0., LON+360., LON)
clrg = [(0.95,0.95,0.95)]
ax1.contour(LON1,list(range(100,181,10)),
          colors=clrg,
          linestyles='solid', 
          linewidths=1.0)
ax1.contour(LON,list(range(-170,-60,10)),
          colors=clrg,
          linestyles='solid',
          linewidths=1.0)
ax1.contour(LAT,list(range(0,89,10)),
          colors=clrg,
          linestyles='solid',
          linewidths=1.0)


stl = 'GOFS3.1 grid resolution ||dl||$_2$, km, NEP region'
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


btx = 'plot_resolution_gofs.py'
bottom_text(btx)





