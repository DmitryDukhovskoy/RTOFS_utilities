"""
  Read rtofs output fields
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
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
#from mpl_toolkits.basemap import Basemap, shiftgrid

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')

rdate   = '20221103'
pthbin  = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/FOR_DMITRY/'+'rtofs.'+rdate+'/'
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

flnm = 'rtofs_glo.t00z.n-24.archv'
fina = pthbin+flnm+'.a'
finb = pthbin+flnm+'.b'

get_topo = True
ftopo = 'regional.depth'
fgrid = 'regional.grid'


import mod_read_hycom
#importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo, read_hycom, read_topo
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text

print('Processing '+fina)


def print_1D(A,wd=8,prc=2):
  ndim1 = A.shape[0]
  for k in range (ndim1):
    print('{0}: {1:{width}.{precis}f}'.format(k+1,A[k],width=wd,precis=prc))


huge = 1.e20
rg   = 9806.
fld  = 'thknss'
F,nn,mm,ll = read_hycom(fina,finb,fld,rLayer=1)
F[F>huge] = np.nan
F = F/rg
F[np.where(F<0.001)] = 0.

# For example 30S, 30W is  i=3199 j=1112
dH = np.zeros((ll,mm,nn))
dH[0,:,:] = F
for kk in range(2,ll+1):
  F, nn, mm, ll = read_hycom(fina,finb,fld,rLayer=kk)
  F = F/rg
  F[np.where(F>1.e20)] = np.nan
  F[np.where(F<0.001)] = 0.
  dH[kk,:,:] = F

ZZ, ZM = zz_zm_fromDP(dH)

if get_topo:
#  HH = read_topo(pthgrid,ftopo,nn,mm)
  LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)
  get_topo = False  


#plt.ion()
#plt.clf()
#plt.pcolormesh(T)

# Plot J-section
i1 = 3050
i2 = 3250
j1 = 1111

ZZs = ZZ[:,j1,i1:i2]
Hb  = HH[j1,i1:i2]
XX  = LON[j1,i1:i2]
y0  = LAT[j1,i1]
Bmin = np.floor(np.min(Hb))-100.

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
#ax1.plot(XX,Hb)
# Patch bottom:
verts = [(np.min(XX),-8000),*zip(XX,Hb),(np.max(XX),-8000)]
poly = Polygon(verts, facecolor='0.6', edgecolor='0.6')
ax1.add_patch(poly)

# Plot interfaces
for kk in range(ll):
  dmm = ZZs[kk,:]
  ax1.plot(XX,dmm,color=[0,0,0])

#ax1.autoscale(enable=True, axis='both', tight='True')
ax1.set_xlim([np.min(XX), np.max(XX)])
ax1.set_ylim([Bmin, 0])


stl = ('Intrf depths, lat={0:5.2f}, j={1}, i={2}:{3}'.format(y0,j1+1,i1+1,i2+1))
ax1.set_title(stl)

btx = 'read_rtofs.py'
bottom_text(btx)


dZZ = np.diff(ZZs, axis=0)
print_1D(dZZ[:,j0])

