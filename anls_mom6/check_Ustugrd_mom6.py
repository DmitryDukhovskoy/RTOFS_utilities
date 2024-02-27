"""
  Plot sections with vertical distribution of T or S

  MOM6 grid

  * ------------ V(i,j) -------------- *
  |                                    |  
  |                                    |  
  |                                    |  
  |                                    |  
  |                                    |  
 U(i-1,j)          * h(i,j)          U(i,j)
  |                                    |  
  |                                    |  
  |                                    |  
  |                                    |  
  |                                    |  
  |                                    |  
  * ------------ V(i,j-1) ------------ *


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


expt    = '003'
YR      = 2020
MM      = 1
DD      = 1
jday    = int(mtime.date2jday([YR,MM,DD]))
HR      = 12
hg      = 1.e15
fld     = 'u'  # 

XSCT = mutil.rtofs_sections()
SCTnames = list(XSCT.keys())


pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'

dnmb = mtime.jday2dnmb(YR,jday)
DV   = mtime.datevec(dnmb)
if not MM:
  MM   = DV[1]
  DD   = DV[2]


pthbin = pthrun + 'tarmom_{0}{1:02d}/'.format(YR,MM)
flmom  = 'ocnm_{0}_{1:03d}_{2}.nc'.format(YR,jday,HR)
flin   = pthbin + flmom


import mod_mom6 as mom6util
importlib.reload(mom6util)
pthgrid  = pthrun + 'INPUT/'
fgrd_mom = pthgrid + 'regional.mom6.nc'
ftopo_mom= pthgrid + 'ocean_topog.nc'
LON, LAT = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
HH       = mom6util.read_mom6depth(ftopo_mom)
Lmsk     = mom6util.read_mom6lmask(ftopo_mom)

idm, jdm, kdm = mom6util.mom_dim(flin)
ijdm          = idm*jdm

# the height of the water column = D + ssh
# hbt = sum (dH)
# HH = hbt - ssh
#hbt = mom6util.read_mom6(flin, 'hbt')
# Read layer thicknesses:
dH   = np.zeros((kdm,jdm,idm))
rfld = 'h' 
for kk in range(1,kdm+1):
  A2D = mom6util.read_mom6(flin, rfld, rLayer=kk)
  dH[kk-1,:,:] = A2D

ssh    = mom6util.read_mom6(flin, 'SSH')
ZZ, ZM = mom6util.zz_zm_fromDP(dH, ssh, f_btm=False)
#ZZ, ZM = mom6util.get_zz_zm(flin)

# Read 3D fields:
U2D = mom6util.read_mom6(flin, 'u', rLayer=1)
V2D = mom6util.read_mom6(flin, 'v', rLayer=1)

# 
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

ax1.set_xlim([600, 1000])
ax1.set_ylim([800,1100])

f_setrmu = False
if f_setrmu:
# Bind the button_press_event with the onclick() method
  fig1.canvas.mpl_connect('button_press_event', onclick)

# Pick the near-coast ocean point, coast to the right and above:
# u-component should be nan, v - nan, depth < 0
i=820
j=1002

print('U={0}'.format(U2D[j,i]))
print('V={0}'.format(V2D[j,i]))
print('HH={0}'.format(HH[j,i]))

# Coast to the right only, u - nan, v - not nan
i=820
j=1001
print('U={0}'.format(U2D[j,i]))
print('V={0}'.format(V2D[j,i]))
print('HH={0}'.format(HH[j,i]))



