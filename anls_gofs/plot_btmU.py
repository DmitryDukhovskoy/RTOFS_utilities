"""
  Plot
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

YR = 1994
MM = 1

# Western coast for GLBv grid:
is1 = 650
is2 = 886
js1 = 1775
js2 = 2200


pthoutp    = '/work/Dmitry.Dukhovskoy/data/gofs31_btmu_westcoast/'
#ftopo_regn = "gofs31_GLBv008_topo11_westcoast.pkl"
ftopo_regn = "gofs31_GLBv008_11_westcoast.pkl"
dftopo_regn = os.path.join(pthoutp,ftopo_regn)
print(f"Loading region topo, grid --> {dftopo_regn}")
with open(dftopo_regn,'rb') as fid:
  lon, lat, HH = pickle.load(fid)

LON, LAT = np.meshgrid(lon, lat)

floutp  = f"gofs31_53X_btmuv_westcoast_{YR}{MM:02d}.pkl"
dfloutp = os.path.join(pthoutp,floutp)
print(f"Loading mean UV --> {dfloutp}")
with open(dfloutp,'rb') as fid:
  UBTM,VBTM = pickle.load(fid)



cmpr = mutil.colormap_ssh(nclrs=200)
rmin = -0.2
rmax = 0.2
cmpr.set_bad(color=[0.1, 0.1, 0.1])

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
im1 = ax1.pcolormesh(VBTM, \
                 cmap=cmpr,\
                 vmin=rmin, \
                 vmax=rmax)

ax1.axis('scaled')
#ax1.set_xlim([600, is2])
#ax1.set_ylim([js1, js2])


 
