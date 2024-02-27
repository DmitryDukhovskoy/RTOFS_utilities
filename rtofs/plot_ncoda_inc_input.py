"""
  Plot NCODA increment input fields
  used in ncoda_archv_inc.f
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

#rdate   = '2022061400'
rdate  = '2022111700'
pthbin = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/hycom/ncoda_archv_inc/'
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

btx = 'plot_ncoda_inc_input.py'


IDM = 4500
JDM = 3298
# Use either original file names dumped from NCODA or
# soft-links names created in ncoda_archv shell
# but check directories
#flnm = 'lyrprs_pre_1o{0}x{1}_{2}_0000_analinc'.format(IDM,JDM,rdate)
flnm = 'lypr_'+rdate+'_analinc'
fina = pthbin+flnm

get_topo = True
ftopo = 'regional.depth'
fgrid = 'regional.grid'


import mod_read_ncoda as rncoda
#importlib.reload(rncoda)
from mod_read_hycom import read_grid_topo, read_hycom, read_topo
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text

print('Processing '+fina)


huge = 1.e20
rg   = 9806.
onem = 1.0198*rg  # NCODA conversion to Pa ???

# Read NCODA depths:
ZI, kzi, ZM, kzm = rncoda.ncoda_depths()
KDM = kzm
Lrprs = np.array([])
for rlr in range(1,KDM+1):
  F = rncoda.read_ncoda_inc(fina,IDM,JDM,KDM,rLayer=rlr)
#F[F<1.e-10] = np.nan     # no increments
#F = F/rg
  F[np.where(F<1.e-10)] = 0.
  if Lrprs.size == 0:
    Lrprs = F.copy()
    Lrprs = np.expand_dims(Lrprs,axis=0)
  else:
    F = np.expand_dims(F, axis=0)
    Lrprs = np.append(Lrprs, F, axis=0)

if get_topo:
#  HH = read_topo(pthgrid,ftopo,nn,mm)
  LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)
  get_topo = False  


# Plot 2D map of NCODA pressure increments:
f_pltrdp = True
if f_pltrdp:
  rlr = 1
  zlr = ZM[rlr-1]
 
  F = np.squeeze(Lrprs[rlr-1,:,:])
  fgnmb1 = 1
  print('Plotting NCODA incr '+flnm)
#  clrmp = copy(plt.cm.afmhot_r)
  clrmp = copy(plt.cm.gist_ncar_r)
  clrmp.set_bad(color=[0.7,0.7,0.7])
  clrmp.set_under(color=[1,1,1])
  plt.ion()
  fig1 = plt.figure(fgnmb1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im1 = plt.pcolormesh(F, vmin=1, vmax=201, cmap=clrmp)
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

  ctl = 'Layer pressure incr {0:5.1f}m, {1}, {2}'.\
        format(zlr,rdate,flnm)
  ax1.set_title(ctl)

  bottom_text(btx)



# ==========================
# Plot selected W-E profiles  
import plot_sect as psct
importlib.reload(psct)

SCTP = ["GoM1","NAtl1"]
nsct = len(SCTP)
f_pltxsct = True
if f_pltxsct:
  nfg = 0
  for ii in range(nsct):
    nfg += 1
    xsct = SCTP[ii]
    SCT = psct.Jsections()

    stl = ('ncoda pr incr input {0} {1}'.\
           format(xsct,flnm))
    Hb, XX = psct.plot_rho_Jsct(xsct, Lrprs, ZM, HH, LON, LAT,\
                                fgnmb=nfg, stl=stl, sct_show=True, \
                                rmin=1., rmax=100., btx=btx)



