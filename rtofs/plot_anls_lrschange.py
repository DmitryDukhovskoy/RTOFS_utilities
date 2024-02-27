"""
  Plot analyze layer changes 
  after incremental update
  to identify problem with HYCOM blow up
  data extracted in anls_lrs_change.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import pickle
import datetime
import matplotlib.colors as colors
import matplotlib.mlab as mlab

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text

import mod_utils as mutil
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_read_ncoda as rncoda
importlib.reload(rncoda)
#import mod_extrargo as exargo
#importlib.reload(exargo)
from mod_utils import tsarray

np.set_printoptions(precision=3)

# Cycle over N days comparing layer characteristics
# how they change from day1 to day N
rdateS = '20230303'  # day to start
rdateE = '20230331'  # day to end
expt   = 'paraD'
sfx    = 'n-24'
f_save = True      # save statistics


# check location:
iPrf = 575
jPrf = 1607
#
# Region : Sulu Sea, deep basin
II = np.array([539, 567, 586, 608, 613, 616, 611, 602, 584, 576, \
               569, 561, 544, 536])
JJ = np.array([1610, 1638, 1641, 1642, 1632, 1612, 1605, 1601, 1584, 1578, \
               1576, 1576, 1588, 1605])


plt.ion()

#importlib.reload(mmisc)
import mod_time as mtime
import mod_read_hycom as mhycom
importlib.reload(mhycom)
from mod_utils_fig import bottom_text

pthdump = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_{0}/data_anls'.\
          format(expt)

pthgrid= '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
lonPrf = LON[jPrf,iPrf]
latPrf = LAT[jPrf,iPrf]


# Function to print mouse click event coordinates
def onclick(event):
   print([event.xdata, event.ydata])

f_setrgn = False
if f_setrgn:
# Bind the button_press_event with the onclick() method
  fig1.canvas.mpl_connect('button_press_event', onclick)

dnmbS = int(mtime.rdate2datenum(rdateS))
dnmbE = int(mtime.rdate2datenum(rdateE))


#fldump = pthdump + '/lrthknss_anls.pkl'
fldump  = pthdump + '/lrthknss_anls'+rdateS+'.pkl'

print('Loading analyzed data to '+fldump)
with open(fldump,'rb') as fid:
  RMMX = pickle.load(fid)

  RMIN = RMMX.rmin
  RMAX = RMMX.rmax
  IMIN = RMMX.imin
  JMIN = RMMX.jmin
  IMAX = RMMX.imax
  JMAX = RMMX.jmax
  TM   = RMMX.time

# Ignore 1st empty column
if TM[0] == 0:
  RMIN = RMIN[:,1:]
  RMAX = RMAX[:,1:]
  IMIN = IMIN[:,1:]
  JMIN = JMIN[:,1:]
  IMAX = IMAX[:,1:]
  JMAX = JMAX[:,1:]
  TM   = TM[1:]
   
a1 = RMAX.shape[0]
a2 = RMAX.shape[1]
dtime = TM-TM[0]
xl1 = dtime[0]
xl2 = dtime[-1]
hlrs = np.arange(0,a1)

from copy import copy
#clrmp = copy(plt.cm.Spectral_r)
clrmp = copy(plt.cm.gist_stern_r)
clrmp.set_bad(color=[0.7,0.7,0.7])

fgnmb = 1
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.3, 0.8, 0.55])
im1 = ax1.pcolormesh(dtime,hlrs,RMAX, vmin=0., vmax=1., cmap=clrmp)
#ax1.axis('scaled')
ax1.set_xlim([xl1,xl2])
ax1.set_ylim([0, 25])
ax1.invert_yaxis()
ax1.set_xlabel('Days')
ax1.set_ylabel('HYCOM Layers')


# Colorbar
ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='max')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

dnmb0 = TM[0]
dv0 = mtime.datevec(dnmb0)
ctl = 'layer thkss change (dp(new)-dp(old))/H \n'
ctl = ctl + 'Day0={0}/{1:02d}/{2:02d}'.format(dv0[0],dv0[1],dv0[2])
ax1.set_title(ctl)

btx = 'plot_anls_lrschange.py'
bottom_text(btx)






