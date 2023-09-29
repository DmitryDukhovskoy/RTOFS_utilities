"""
  Check min/max values in the water column
  Identify extreme T,S, U/V values
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
import importlib
import struct
import pickle
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
import mod_read_ncoda as rncoda

rdate   = '20220616'
fldplt = 'temp'  # temp salin u-vel. v-vel.

pthbs   = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/'
pth1    = pthbs+'hycom/ncoda_archv_inc/'
pthbin  = pth1
pthbgr  = pthbs+'hycom/ncoda_archv_inc/'  # background field dir
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
pthrsrt = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/rtofs_expts/'+\
          'rtofs_para8.1/navy_restart/'

jday = rncoda.rdate2julian(rdate)
yr, mo, mday, hr = rncoda.parse_rdate(rdate)

# Select - updated or background fields:
f_updated = False
# flupdt - updated fields=background+increments output from ncoda_archv_inc
# flbgr  - background field from HYCOM f/cast
flupdt = 'archv_1_inc.{0}_{1:03d}_00'.format(yr,jday)
flbgr   = 'archv.{0}_{1:03d}_00'.format(yr,jday)

IDM = 4500
JDM = 3298
KDM = 41

get_topo = True
ftopo = 'regional.depth'
fgrid = 'regional.grid'

import mod_read_hycom
#importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo, read_hycom, \
                           read_topo
import mod_utils as utls
#importlib.reload(utls)
from mod_read_hycom import read_2Dbinary
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text

btx = 'plot_minmax.py'

#  HH = read_topo(pthgrid,ftopo,nn,mm)
LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)

huge = 1.e20
rg   = 9806.

# Read layer depths
if f_updated:
# Updated fields with increments:
  fina = pthbin+flupdt+'.a'
  finb = pthbin+flupdt+'.b'
else:
# Background fields:
  fina = pthbin+flbgr+'.a'
  finb = pthbin+flbgr+'.b'

print('Reading layer thickness: '+fina)
fld  = 'thknss'
F,nn,mm,ll = read_hycom(fina,finb,fld,rLayer=1)
F[np.where(F>huge)] = np.nan
F = F/rg
F[np.where(F<0.001)] = 0.

dH = np.zeros((ll,mm,nn))
dH[0,:,:] = F
for kk in range(2,ll+1):
  F,nn,mm,lmm = read_hycom(fina,finb,fld,rLayer=kk)
  F = F/rg
  F[np.where(F>huge)] = np.nan
  F[np.where(F<0.001)] = 0.
  dH[kk-1,:,:] = F

ZZ, ZM = zz_zm_fromDP(dH, f_btm=False)
kdm = ZM.shape[0]
jdm = ZM.shape[1]
idm = ZM.shape[2]

print('Processing '+fina)

fld = fldplt
A3d = np.array([])
for kk in range (1,KDM+1):
  F,n1,m1,l1 = read_hycom(fina,finb,fld,rLayer=kk)
  F[np.where(F>huge)] = np.nan
  if A3d.size == 0:
    A3d = F.copy()
    A3d = np.expand_dims(A3d, axis=0)
  else:
    F = np.expand_dims(F, axis=0)
    A3d = np.append(A3d, F, axis=0)

# plt.clim(a,b) = limits of colorbar
# Check lr thkn gradient:
import mod_utils as mutil
importlib.reload(mutil)

maxA = mutil.find_max3d(A3d,dH)


if fldplt == 'u-vel.' or fldplt == 'v-vel.':
  rmin, rmax = psct.minmax_clrmap(maxA,pmin=5.,pmax=80.,cpnt=0.01,fsym=True)
else:
  rmin, rmax = psct.minmax_clrmap(maxA,pmin=10.,pmax=80.,cpnt=0.1)
rmin = 0.

eps_max = 50.
J,I = np.where(maxA>eps_max)

# Plot max 
fgnmb1 = 1
print('Plotting max map ...')
clrmp = copy(plt.cm.afmhot_r)
clrmp.set_bad(color=[0.7,0.7,0.7])
plt.ion()
fig1 = plt.figure(fgnmb1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
im1 = plt.pcolormesh(maxA,shading='flat',\
                     vmin=0., vmax=50, cmap=clrmp)
plt.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
plt.plot(I,J,'r*')

ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='max')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)

plt.sca(ax1)
ctl = 'max({0}) over Z, {1} rdate={2}'.format(fld,fina[-42:],rdate)
ax1.set_title(ctl)

nmax = min(20,I.shape[0])
ss1  = 'Extreme Values >{0:5.1f} \n I       J\n'.format(eps_max)
for kk in range(nmax):
  ss1 = ss1+'{0}  {1}\n'.format(I[kk],J[kk])


ax1.text(100,10,ss1)

bottom_text(btx)





