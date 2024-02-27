"""
  Plot transports computed in calc_transp_line.py
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

expt    = '003'
hg      = 1.e15
sctnm   = 'Fram79'
dnmb1 = mtime.datenum([2021,1,1])
dnmb2 = mtime.datenum([2021,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)


pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'
pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/' + \
          'MOM6_CICE6/expt{0}/'.format(expt)

floutp = 'mom6-{4}_volTrt_{0}{1:02d}-{2}{3:02d}_{5}.pkl'.\
         format(dv1[0], dv1[1], dv2[0], dv2[1], expt, sctnm)
ffout = pthoutp + floutp
print('Loading '+ffout)
with open(ffout,'rb') as fid:
  VTRP = pickle.load(fid)

# dir(VTRP)
# VTRP.__dict__

TM   = VTRP.TM
VFlx = VTRP.trnsp*1.e-6   # Sv
Hbtm = VTRP.Hbtm
XX   = VTRP.xcrd
YY   = VTRP.ycrd
Tday = TM-TM[0]

# Mean flux:
alf  = 10.
mVF  = np.nanmean(VFlx, axis=0)
pL   = np.percentile(VFlx, alf, axis=0)
pU   = np.percentile(VFlx, (100.-alf), axis=0)
#
# Total transport:
VFtot = np.nansum(VFlx, axis=1)
mnVF  = np.mean(VFtot)
stdVF = np.std(VFtot)


plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
ax1.plot(XX, mVF)
ax1.plot(XX, pL, color=[0.8, 0.85, 1.])
ax1.plot(XX, pU, color=[0.8, 0.85, 1.])

xl1 = np.min(XX)
xl2 = np.max(XX)

ax1.set_xlim([xl1, xl2])
ax1.set_xticks(np.arange(np.floor(xl1),np.ceil(xl2),2.))
ax1.grid(True)

ctl = ('0.08 MOM6-CICE6-{0} {1}, depth.intgr. VFlux (Sv),  \n'.format(expt, sctnm)\
        + 'mean & IDR {0}/{1}/{2} - {3}/{4}/{5}'.\
       format(dv1[0], dv1[1], dv1[2], dv2[0], dv2[1], dv2[2]))
ax1.set_title(ctl)

# Plot bottom:
verts = [(xl1-1.,-8000),(xl1-1., Hbtm[0]),*zip(XX,Hbtm),\
         (xl2+1.,Hbtm[-1]),(xl2+1.,-8000)]
poly = Polygon(verts, facecolor='0.', edgecolor='none')

ax2 = plt.axes([0.1, 0.36, 0.8, 0.1])
ax2.cla()
#ax2.plot(XX, Hbtm)
ax2.add_patch(poly)
ax2.set_xlim([xl1, xl2])
ax2.set_yticks(np.arange(-5000.,0.,1000.))
ax2.set_ylim([np.floor(np.min(Hbtm)), 0])
ax2.set_xticks(np.arange(np.floor(xl1),np.ceil(xl2),2.))
ax2.grid(True)

ax3 = plt.axes([0.1, 0.1, 0.6, 0.2])
ax3.plot(Tday, VFtot)
ax3.set_xlim([Tday[0], Tday[-1]])
ax3.set_xticks(np.arange(0,Tday[-1],30))
ax3.grid(True)
ax3.set_title('VolTransport, Sv')
ax3.set_xlabel('Days, {0}'.format(dv1[0]))

ax4 = plt.axes([0.72, 0.15, 0.15, 0.2])
sinfo = 'mean Flux = {0:7.2f} Sv\n'.format(mnVF)
sinfo = sinfo + 'sgm = {0:6.3f} Sv'.format(stdVF)
ax4.text(0., 0., sinfo)
ax4.axis('off')

btx = 'plot_transp_line.py'
bottom_text(btx,pos=[0.02, 0.03])









