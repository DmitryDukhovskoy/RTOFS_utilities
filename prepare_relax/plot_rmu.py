"""
  Plot climatology relaxation HYCOM fields 
  rmu.a
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
import importlib
import struct
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
#from mpl_toolkits.basemap import Basemap, shiftgrid

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
import mod_read_hycom as mrhycom

pthrmu = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_paraCd/relax/rmu/'
#flrmu  = 'rmu_tropAtl_30days'
flrmu  = 'rmu_02days_Carib'

fina = pthrmu + flrmu + '.a'
finb = pthrmu + flrmu + '.b'

try:
  fgb = open(finb,'r')
except:
  print('Could not open '+finb)

fgb.seek(0)
nl0 = 0
while nl0 < 100:
  nl0 += 1
  data = fgb.readline().split()
  if len(data) < 2:
    continue
  adim = data[0]
  ii = adim.find('jdm')
  if ii>0:
    break

fgb.close()
if ii<0:
  raise Exception('No idm found: Reading ',finb)

IDM  = int(data[2])
JDM  = int(data[3])
IJDM = IDM*JDM
npad = 4096-IJDM%4096

huge = 0.01*2.**99
fga = open(fina,'rb')
dmm = np.fromfile(fga, dtype='>f',count=IJDM) 
fga.close()

amin = np.min(dmm[np.where(dmm<huge)])
amax = np.max(dmm[np.where(dmm<huge)])
print('Reading {0} min={1:9.4e} max={2:9.4e}'.\
        format(flrmu,amin,amax))
dmm = dmm.reshape((JDM,IDM))


# Read topo
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo   = 'regional.depth'
fgrid   = 'regional.grid'

LON, LAT, HH = mrhycom.read_grid_topo(pthgrid,ftopo,fgrid)

# Convert relax e-folding time to days 
dmm = np.where(dmm==0.0, np.nan, dmm)   # sec ^-1
Rday = (1./dmm)*(1./86400.)             # days

clrmp = copy(plt.cm.gist_ncar_r)
clrmp.set_bad(color=[0.7,0.7,0.7])

Rday = np.where(np.isnan(Rday),0.0,Rday)
Rday[np.where(HH >= 0.)] = np.nan
Rmax = np.ceil(np.nanmax(Rday))

plt.ion()
fgnmb = 1
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])

im1 = ax1.pcolormesh(Rday, vmin=0, vmax=Rmax, cmap=clrmp)
ax1.axis('equal')
ax1.set_xlim([2300,3200])
ax1.set_ylim([1500,2100])

# Colorbar:
ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='max')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.tick_params(direction='in', length=12)





