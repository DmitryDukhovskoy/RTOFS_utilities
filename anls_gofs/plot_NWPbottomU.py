"""
  Derive coastal indices along the western coast to plot
  Average bottom velocities along the western US coast
  GOFS reanalysis
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import time
import timeit
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap

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
import mod_gofs31 as mgofs

# Region domain
lat1 = 22.
lat2 = 48.
lon1 = 235.42
lon2 = 251.

pfld = 'ubtm'  # u or v component to plot
#YR = 1994
#MM = 1
# Average over years YR1-YR2 and months = MM1-MM2
YR1 = 1994
YR2 = 2007
MM1 = 1
MM2 = 3

# Western coast for GLBv grid:
is1 = 650
is2 = 886
js1 = 1775
js2 = 2200


pthoutp    = '/work/Dmitry.Dukhovskoy/data/gofs31_btmu_westcoast/'
ftopo_regn = "gofs31_GLBv008_topo11_westcoast.pkl"
#ftopo_regn = "gofs31_GLBv008_11_westcoast.pkl"
dftopo_regn = os.path.join(pthoutp,ftopo_regn)
print(f"Loading region topo, grid --> {dftopo_regn}")
with open(dftopo_regn,'rb') as fid:
  lon, lat, HH = pickle.load(fid)

LON, LAT = np.meshgrid(lon, lat)
DX, DY   = mhycom.dx_dy(LON,LAT)

# Read saved monthly bottom U,V:
icc = 0
for YR in range(YR1,YR2+1):
  for MM in range(MM1, MM2+1):
    floutp  = f"gofs31_53X_btmuv_westcoast_{YR}{MM:02d}.pkl"
    dfloutp = os.path.join(pthoutp,floutp)
    print(f"Loading mean UV --> {dfloutp}")
    with open(dfloutp,'rb') as fid:
      UB0,VB0 = pickle.load(fid)

    icc += 1
    if icc == 1:
      VBTM = VB0
      UBTM = UB0
    else:
      VBTM = VBTM+VB0
      UBTM = UBTM+UB0

UBTM = UBTM/icc
VBTM = VBTM/icc

# Inland contour:
#import mod_anls_gofs as manlsgofs
# Search coast line then go off shore until z=z0
import mod_misc1 as mmisc
import mod_anls_gofs as mag
z0   = -400.
INDX, JNDX, DL = mag.find_NWP_isobath(HH, LON, LAT, DX, z0)
dI       = INDX[:,1] - INDX[:,0] + 1
nIX      = INDX.shape[0]
ndi      = np.max(dI)

# Meridional bottom v on the shelf:
Vsh  = np.zeros((nIX, ndi))-1.e30
Ush  = np.zeros((nIX, ndi))-1.e30
Dst  = np.zeros((nIX, ndi))        # offshore distance
Ylat = np.zeros((nIX, ndi))        # offshore distance

for ikk in range(nIX):
  i1 = INDX[ikk,0]
  i2 = INDX[ikk,1]
  j1 = JNDX[ikk]
  aa = VBTM[j1,i1:i2+1]
  uu = UBTM[j1,i1:i2+1]
  dl = -(np.arange(ndi)+1)*DL[ikk]
  dl = np.flip(dl)

  naa = len(aa) 
  Vsh[ikk,-naa:]  = aa
  Ush[ikk,-naa:]  = uu
  Dst[ikk,:]  = dl*1e-3
  Ylat[ikk,:] = LAT[j1,i2-ndi:i2]



cmpr = mcmp.colormap_uv(nclrs=200)
rmin = -0.1
rmax = 0.1
cmpr.set_bad(color=[0.1, 0.1, 0.1])
cmpr.set_under(color=[0.7, 0.7, 0.7])

if pfld == 'ubtm':
  A2d  = Ush
  cfld = 'U'
elif pfld == 'vbtm':
  A2d  = Vsh
  cfld = 'V'

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.5, 0.8])
im1 = ax1.pcolormesh(Dst, Ylat, A2d, \
                 cmap=cmpr,\
                 vmin=rmin, \
                 vmax=rmax)
#ax1.axis('scaled')
#ax1.set_xlim([600, is2])
#ax1.set_ylim([js1, js2])

stl = f"GOFS3.1-53.X Reanalysis {YR1}-{YR2} M:{MM1}-{MM2} bottom {cfld} m/s"
ax1.set_title(stl)
ax1.set_xlabel('Offshore distance, km')
ax1.set_ylabel('Latitude')

ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2)
#clb = plt.colorbar(im1, cax=ax2, extend='max')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=5)

btx = 'plot_NWPbottomU.py'
bottom_text(btx,pos=[0.02, 0.03])




 
