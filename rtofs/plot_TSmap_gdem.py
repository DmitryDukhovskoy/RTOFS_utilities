"""
  Plot GDEM from binary 
  see code:
  /scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/ncoda/sorc/ncoda_qc/libsrc/prfobs
  gdem_mod.f
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import pickle
import importlib
import datetime
import matplotlib.colors as colors
import matplotlib.mlab as mlab

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text

import mod_time as mtime
import mod_read_hycom as mhycom
import mod_utils as mutil
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_read_ncoda as rncoda
#importlib.reload(rncoda)
import mod_gdem as mgdem
importlib.reload(mgdem)
import mod_utils_fig as mufig


imo = 10     # Month
fld = 'salt' # salt or temp
z0  = -150.  # Depth to plot

pthgdem = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/ncoda/fix/gdem/'

LON, LAT, ZM = mgdem.gdem_grid()
NI = LON.shape[0]
NJ = LAT.shape[0]
NK = ZM.shape[0]

LON = np.where(LON > 180., LON-360., LON)

A3d = mgdem.read_gdem3d(pthgdem, imo, fld, NI, NJ, NK)

# Load GDEM land mask:
pthout   =  '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/GDEM/'
gdem_out = pthout + 'GDEM_sealand_mask.pkl'
print('Loading GDEM land mask --> ' + gdem_out)

with open(gdem_out,'rb') as fid:
  Lmsk = pickle.load(fid)

#
# Find top/btm layers for interpolation if needed
Itop = max(np.where(ZM >= z0)[0])
Ibtm = Itop + 1

if abs(Itop - z0) < 1.e-3:
  Aint = np.squeeze(A3d[Itop,:,:])
else:
  Atop = np.squeeze(A3d[Itop,:,:])
  Abtm = np.squeeze(A3d[Ibtm,:,:])
  Ztop = ZM[Itop]
  Zbtm = ZM[Ibtm]

# linear interp
# Lagrange polynomial:
  Aint = Atop*((z0-Zbtm)/(Ztop-Zbtm)) + \
         Abtm*((z0-Ztop)/(Zbtm-Ztop))

Aint = np.where(Lmsk == 0, np.nan, Aint)

def find_indx(XX, x1):
  dd = np.abs(XX-x1)
  ix = np.where(dd==np.min(dd))[0][0]
  return ix

import mod_utils as mutil
importlib.reload(mutil)

if fld == 'salt':
  cmpr = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  cmpr.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 35.
  rmax = 37.
  fld1 = 'mapS'
elif fld == 'temp':
  cmpr = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  cmpr.set_bad(color=[0.2,0.2,0.2])
  rmin = -2.
  rmax = 28.
  fld1 = 'mapT'


# Select region:
lon1  = -101.84
lon2  = -29.84
lat1  = 3.7
lat2  = 42.97
xlim1 = find_indx(LON, lon1)
xlim2 = find_indx(LON, lon2)
ylim1 = find_indx(LAT, lat1)
ylim2 = find_indx(LAT, lat2)

Asub = Aint[ylim1:ylim2+1,xlim1:xlim2+1]
alf=2.
#rmin, rmax = mutil.minmax_clrmap(Asub, pmin=alf, pmax=100-alf, cpnt=0.01)
rmin = 34.79
rmax = 37.03

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))

plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
im1 = ax1.pcolormesh(Aint, \
               cmap=cmpr,\
               vmin=rmin, \
               vmax=rmax)
#ax1.contour(HH, [0.0], colors=[(0,0,0)], linewidths=1)
ax1.axis('scaled')
ax1.set_xlim([xlim1,xlim2])
ax1.set_ylim([ylim1,ylim2])

ctl = 'GDEM {0}, {1}, Mo={2}'.format(fld, z0, imo)
ax1.set_title(ctl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='both')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.tick_params(direction='in', length=12)

# Vert profiles:
ax3  = plt.axes([0.1, 0.08, 0.3, 0.35])
i0   = 1125
j0   = 377
zlim = -500.

S1  = np.squeeze(A3d[:,j0,i0])
zzm = ZM
si = Aint[j0,i0]

if min(zzm) <= zlim:
  izlim = np.min(np.where(zzm <= zlim )) + 1
else:
  izlim = np.max(np.where(zzm == np.nanmin(zzm)))

xl1 = np.floor(np.nanmin(S1[:izlim])*10.)/10.
xl2 = np.ceil(np.nanmax(S1[:izlim])*10.)/10.

ax3.plot(S1,zzm,'.-')
ax3.plot(si,z0,'r.')
ax3.set_xlim(xl1,xl2)
ax3.set_ylim(-500,0)
ax3.plot([xl1, xl2],[z0, z0],'r--')
ax3.grid(True)
ax3.set_title('S i,j={0}, {1}'.format(i0,j0))

btx = 'plot_TSmap_gdem.py'
mufig.bottom_text(btx, pos=[0.1, 0.03])






