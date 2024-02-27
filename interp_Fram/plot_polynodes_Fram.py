"""
  Plot best nodes for polynomial
  Fram strait
  Data saved in obsloc_Fram_FWyr.py
"""
import os
import numpy as np
from copy import copy
import importlib
import matplotlib.pyplot as plt
import sys
import pickle

sys.path.append('/home/ddmitry/codes/MyPython/hycom_utils')
sys.path.append('/home/ddmitry/codes/MyPython/draw_map')
sys.path.append('/home/ddmitry/codes/MyPython')

plt.close('all')

from mod_utils_fig import bottom_text
plt.ion()  # enables interactive mode

import mod_polynom as mpol
importlib.reload(mpol)

Mplnm = 'MMXI' # method to find best polynom nodes
strnm = 'FramStr'
#FlxNm = 'FWFlux' 
FlxNm = 'VolFlux'
Nnds = 17 # # of nodes

pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/'
pthtopo  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/'


btx = 'plot_polynodes_Fram.py'

#foutp = pthout + 'polynodes_' + strnm + '_' + Mplnm + '.pkl'
foutp = pthout + 'polynodes_' + strnm + '_' + FlxNm + '_' + Mplnm + '.pkl'
print('Loading ---> ' + foutp)
with open(foutp,'rb') as fid:
  PNOD = pickle.load(fid)

# Read xsection coordinates:
flcoord = 'arc08_sect_coord_FramStr.dat'
dL, Lonx, Latx, Indx, Jndx = mpol.read_coord(pthout, flcoord)

import mod_read_hycom as mhycom
importlib.reload(mhycom)

flgrid = 'depth_ARCc0.08_11.nc'
LON, LAT, HH = mhycom.readnc_grid_topo(pthtopo,flgrid)

#
# Extract polynom indices:
# For different years
mnIX, stdIX = mpol.extract_indx_nodes(PNOD,Nnds)
mnIX = np.round(mnIX).astype(int)


from matplotlib import cm
clrs = cm.get_cmap('winter',200)
clrs.set_over(color=[0., 0., 0.])
rmin = -5000.
rmax = 0.

fig1 = plt.figure(1, figsize=(9,9))
plt.clf()

ax1 = plt.axes([0.1, 0.2, 0.7, 0.7])
im1 = ax1.pcolormesh(HH, shading='flat', 
                     cmap = clrs,
                     vmin = rmin,
                     vmax = rmax)

xlim1 = 870
xlim2 = 1100
ylim1 = 850
ylim2 = 1040
#  ax1.axis('equal')
ax1.axis('scaled')
ax1.set_xlim([xlim1,xlim2])
ax1.set_ylim([ylim1,ylim2])

# Plot section:
ax1.plot(Indx,Jndx,'-',color=[0.7, 0.7, 0.7])

# Plot preferred mooring locations:
nNd = mnIX.shape[0]
for ik in range(nNd):
  ii0 = mnIX[ik]
  xx0 = Indx[ii0]
  yy0 = Jndx[ii0]

  ax1.scatter(xx0,yy0,color=[0.8,0.3,0], s=20)

stl = 'Optimal Mooring Locations, n={0}, {1}'.format(nNd,FlxNm)
ax1.set_title(stl)

ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.set_xticks([])
ax1.set_yticks([])

btx = 'plot_polynodes_Fram.py'
bottom_text(btx)



