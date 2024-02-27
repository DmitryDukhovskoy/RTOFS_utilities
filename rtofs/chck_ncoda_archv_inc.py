"""
  Quick check of the saved fields
  Plot diagnostic output fields from
  ncoda_archv_inc.f code
  to analyze increments/pressure/density changes
  caused by ncoda incr etc
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

rdate   = '20221117';
pth1    = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/hycom/ncoda_archv_inc'

pthbin  = pth1+'/'
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

frhobg = 'rho_bckgr'    # Rho of background field on NCODA z-levels 

IDM = 4500
JDM = 3298
KDM = 41

get_topo = True
ftopo = 'regional.depth'
fgrid = 'regional.grid'


import mod_read_hycom
#importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo, read_hycom, \
                           read_topo, read_hycom_restart
from mod_read_hycom import read_2Dbinary
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text



def print_1D(A,wd=8,prc=2):
  ndim1 = A.shape[0]
  for k in range (ndim1):
    print('{0}: {1:{width}.{precis}f}'.format(k+1,A[k],width=wd,precis=prc))


huge = 1.e20
rg   = 9806.

# Background rho
fina = pthbin+frhobg+'.a'
finb = ""      # only binary
print('Processing '+fina)

Rho_bg = np.array([])
for kk in range (1,KDM+1):
  F = read_2Dbinary(fina,IDM,JDM,kk)
  F[np.where(F>huge)] = np.nan


get_topo = True
if get_topo:
#  HH = read_topo(pthgrid,ftopo,nn,mm)
  LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)
  get_topo = False  

# ================================
def find_indx_lonlat(x0,y0):
  """
  Find closest grid point to lon/lat coordinate
  """
  if x0 > 180.:
    x0 = x0-360.
  
  dmm = np.sqrt((LON-x0)**2+(LAT-y0)**2)
  jj0, ii0 = np.where(dmm == np.min(dmm))

  return ii0[0], jj0[0]

# ================================
kdm = F.shape[0]
jdm = F.shape[1]
idm = F.shape[2]


f_pltrdp = True
if f_pltrdp:
  fgnmb1 = 1
  print('Plotting rho ...')
  clrmp = copy(plt.cm.afmhot_r)
  clrmp.set_bad(color=[0.7,0.7,0.7])
  plt.ion()
  fig1 = plt.figure(fgnmb1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im1 = plt.pcolormesh(F, vmin=0, vmax=40, cmap=clrmp)
  plt.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
 
  ax2 = fig1.add_axes([ax1.get_position().x1+0.02, 
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='max')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  plt.sca(ax1)

  ctl = 'Inverse min(lr_thkns/dH) wrt {0:3.2f}, {1}'.\
        format(rdp0,rdate+'/'+flnm)
  ax1.set_title(ctl)
 
  btx = 'plot_ncoda_inc.py'
  bottom_text(btx)


# ==========================
# Plot selected W-E profiles  
import plot_sect as psct
importlib.reload(psct)

f_plt = False
if f_plt:
  jp0 = 48
  xsct = "GoM1"
#  xsct = "NAtl1"
  ZZs, Hb, XX, dZZs = psct.plot_Jsection(xsct,ZZ,HH,LON,LAT,rdate,sct_show=ax1,\
                      fgnmb=5,dcntr=1, zstart=-1000., jlrs=-1)
  print_1D(dZZs[:,jp0])
 






