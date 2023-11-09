"""
  Pick location on H map
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
from copy import copy
import matplotlib.colors as colors
from matplotlib.patches import Polygon

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_cice6_utils as mc6util
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil

import mod_mom6 as mom6util
importlib.reload(mom6util)
pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_003/'
pthgrid  = pthrun + 'INPUT/'
fgrd_mom = pthgrid + 'regional.mom6.nc'
ftopo_mom= pthgrid + 'ocean_topog.nc'
LON, LAT = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
HH       = mom6util.read_mom6depth(ftopo_mom)
Lmsk     = mom6util.read_mom6lmask(ftopo_mom)

# Function to print mouse click event coordinates
def onclick(event):
   print([event.xdata, event.ydata])

# Plot section
plt.ion()
fig1 = plt.figure(1,figsize=(9,8), constrained_layout=False)
fig1.clf()
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8],)
ax1.contour(HH, [0], linestyles='solid',
            colors=[(0,0,0)], linewidths=1)

ax1.set_xlim([3000, 3800])
ax1.set_ylim([2600,3296])

# Bind the button_press_event with the onclick() method
fig1.canvas.mpl_connect('button_press_event', onclick)




