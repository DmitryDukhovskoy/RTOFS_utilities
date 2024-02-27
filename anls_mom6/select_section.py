"""
  Select sections for vertical T/S/U
  Pick vertices for multi-segment sections
  from an interactive map

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pdb
import importlib
import time
#import struct
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/mom6_utils')

import mod_time as mtime
from mod_utils_fig import bottom_text

import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)

# Function to print mouse click event coordinates
def onclick(event):
   print([event.xdata, event.ydata])

# Check orientation of the line/norms 
dday  = 5       # time stepping for data processing/analysis
dnmb1 = mtime.datenum([2021,1,1])
dnmb2 = mtime.datenum([2021,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)

expt  = '003'
hg    = 1.e15

pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'

import mod_misc1 as mmisc
import mod_mom6 as mom6util
importlib.reload(mom6util)

pthgrid   = pthrun + 'INPUT/'
fgrd_mom  = pthgrid + 'regional.mom6.nc'
ftopo_mom = pthgrid + 'ocean_topog.nc'
LON, LAT  = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
HH        = mom6util.read_mom6depth(ftopo_mom)
Lmsk      = mom6util.read_mom6lmask(ftopo_mom)

# Plot section
plt.ion()
fig1 = plt.figure(1,figsize=(9,8), constrained_layout=False)
fig1.clf()
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8],)
ax1.contour(HH, [0], linestyles='solid',
            colors=[(0,0,0)], linewidths=1)

ax1.contour(HH, [-500], linestyles='solid',
            colors=[(0.7,0.7,0.7)], linewidths=1)
ax1.contour(HH, [-1000], linestyles='solid',
            colors=[(0.8,0.8,0.8)], linewidths=1)

ax1.axis('scaled')
#ax1.set_xlim([600, 1000])
#ax1.set_ylim([800,1100])

# Bind the button_press_event with the onclick() method
fig1.canvas.mpl_connect('button_press_event', onclick)





