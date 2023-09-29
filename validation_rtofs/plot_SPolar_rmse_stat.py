# Plot ssh stat for Southern Ocean using polar projection
# steps: 
# (1) get production & paraXX runs from WCOSS2 (scripts/rtofs)
# (2) Use interpolation indices to do bilinear interpolation 
#     of HYCOM --> CMCES 0.25 grid
#     interp_hycom2cmems.py
# (3) Calculate RMSE (mean 2D and time series):
#     calc_rmse.py for paraD5 and production
#
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
from netCDF4 import Dataset as ncFile
from mpl_toolkits.basemap import Basemap, cm

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text

import mod_utils as mutil
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_time as mtime
import mod_colormaps as mclrs


# Compare parallel vs production (using NRT altimetry SSH as reference)
expt    = 'paraD5'  # paraXX - parallel runs, production is loaded everytime
regn_nm = 'SOcean'
rdateS  = '20230303'
rdateE  = '20230508'
sfx     = 'n-24'

pthsave  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/data_outp/ssh_intrp2cmems/'
frmseout = 'ssh_rmse_' + expt + sfx + '_' + regn_nm + '.pkl'
flrmse   = pthsave + frmseout

# Define region:
import mod_valid_utils as mvutil
importlib.reload(mvutil)
REGNS = mvutil.ssh_regions()
REGNMS = list(REGNS.keys())
xl1    = REGNS[regn_nm]["xl1"]
xl2    = REGNS[regn_nm]["xl2"]
yl1    = REGNS[regn_nm]["yl1"]
yl2    = REGNS[regn_nm]["yl2"]
Ip     = np.array(REGNS[regn_nm]["Ip"])
Jp     = np.array(REGNS[regn_nm]["Jp"])

# Read lon/lat 
pthssh = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/OBS/ssh_altimetry/'
flnm   = 'ssh_duacs_nrt_global_20230301_20230430.nc'
flssh  = pthssh + flnm
nc     = ncFile(flssh,'r')
XA     = nc.variables['longitude'][:].squeeze().data
YA     = nc.variables['latitude'][:].squeeze().data

# RMSE time series
print('Loading  ' + flrmse )
with open(flrmse,'rb') as fid:
  [RMSE_TS, TM] = pickle.load(fid)

# Mean RMSE:
# Parallel run paraD5
frmsemn = 'ssh_rmsemn2D_' + expt + sfx + '_'  + regn_nm + '.pkl'
flrmsemn  = pthsave + frmsemn
print('Loading  ' + flrmsemn )
with open(flrmsemn,'rb') as fid:
  [RMSE_mn, Lmsk] = pickle.load(fid)

# Points inside the region
jdma = RMSE_mn.shape[0]
idma = RMSE_mn.shape[1]
X, Y = np.meshgrid(np.arange(idma), np.arange(jdma))
#Rmsk, IRg, JRg = mmisc.inpolygon_v2(X,Y,Ip,Jp)  # region

# Remove RMSE outside the region:
#RMSE_mn = np.where(Rmsk==1, RMSE_mn, np.nan)
RMSE_mn = np.where(Lmsk==1, RMSE_mn, np.nan)

#
# Load production
frmseout = 'ssh_rmse_product' + sfx + '_' + regn_nm + '.pkl'
flrmse   = pthsave + frmseout
print('Loading  ' + flrmse)
with open(flrmse,'rb') as fid:
  [RMSEp_TS, TMp] = pickle.load(fid)

#frmsemn = 'ssh_rmsemn2D' + sfx + '_product' + '.pkl'
frmsemn = 'ssh_rmsemn2D_product' + sfx + '_'  + regn_nm + '.pkl'
flrmsemn  = pthsave + frmsemn
print('Loading  ' + flrmsemn )
with open(flrmsemn,'rb') as fid:
  [RMSEp_mn, Lmsk] = pickle.load(fid)

# Remove RMSE outside the region:
#RMSEp_mn = np.where(Rmsk==1, RMSEp_mn, np.nan)
RMSEp_mn = np.where(Lmsk==1, RMSEp_mn, np.nan)

#
# Colorbar
from copy import copy
import mod_colormaps as mclrs
CLR = [[1,  1,  1],
       [1,  0.9, 0.8],
       [1,  0.5, 0.],
       [1,  0.,  0.],
       [0.4, 0., 0.]]

CLR = np.array(CLR)
clrmp1 = mclrs.create_colormap(CLR, 200)
clrmp1.set_bad(color=[0.8,0.8,0.8])

Lmsk = Lmsk.astype('float')

# ====================
# Subset S. Ocean region:
sR  = RMSE_mn[0:yl2,:]
sRp = RMSEp_mn[0:yl2,:]
sXA = XA.copy()
sYA = YA[0:yl2]
ny  = sR.shape[0]
nx  = sR.shape[1]
sXX, sYY = np.meshgrid(sXA, sYA)

rmin = 0.
rmax = 0.25

print('Plotting ...')
m = Basemap(projection='spstere',boundinglat=-40,lon_0=0,resolution='l')
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly spaced grid.
x, y   = m(lons, lats) # compute map proj coordinates.
xh, yh = m(sXX, sYY)  # interp. CMCMES coordinates on the projections



# Plot mean RMSE:
plt.ion()

# Plot region and select points if needed:
fig1 = plt.figure(1,figsize=(9,9))
fig1.clf()

ax11 = fig1.add_axes([0.3, 0.55, 0.4, 0.4])
# draw parallels and meridians
parallels = np.arange(-80.,0.,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
meridians = np.arange(-360,0.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
im11 = ax11.pcolormesh(xh, yh, sR, 
                       cmap=clrmp1,
                       vmin=rmin, vmax=rmax)

ax11.set_title('{0} mean RMSE(m) {1}-{2}'.format(expt,rdateS,rdateE))

# Colorbar
cax = fig1.add_axes([0.8, 0.55, 0.015, 0.4])
clb = fig1.colorbar(im11, cax=cax, extend='max')
cax.set_yticklabels(cax.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=5)


# Plot production:
ax12 = fig1.add_axes([0.3, 0.07, 0.4, 0.4])
parallels = np.arange(-80.,0.,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
meridians = np.arange(-360,0.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
im12 = ax12.pcolormesh(xh, yh, sRp,
                     cmap=clrmp1,
                     vmin=rmin, vmax=rmax)

ax12.set_title('Product mean RMSE(m) {0}-{1}'.format(rdateS,rdateE))

# Colorbar
cax2 = fig1.add_axes([0.8, 0.07, 0.015, 0.4])
clb2 = fig1.colorbar(im12, cax=cax2, extend='max')
cax2.set_yticklabels(cax2.get_yticks())
ticklabs = clb2.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb2.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb2.ax.tick_params(direction='in', length=5)

btx = 'plot_SPolar_rmse_stat.py'
bottom_text(btx, pos=[0.02,0.01])


# Plot RMSE time series
dTM = TM - mtime.rdate2datenum('20230301')+1
xt1 = 1
xt2 = np.max(dTM)

md1  = np.percentile(RMSE_TS,50.)
iqr2 = np.percentile(RMSE_TS,90.)
iqr1 = np.percentile(RMSE_TS,10.)

mdP1  = np.percentile(RMSEp_TS,50.)
iqrP2 = np.percentile(RMSEp_TS,90.)
iqrP1 = np.percentile(RMSEp_TS,10.)

mqmin = round(min([iqr1,iqrP1]),2)
mqmax = round(max([iqr2,iqrP2]),2)
mmin  = round(min([min(RMSE_TS),min(RMSEp_TS)]),2)
mmax  = round(max([max(RMSE_TS),max(RMSEp_TS)]),2)

fig2 = plt.figure(2, figsize=(9,8), constrained_layout=False)
fig2.clf()

clr1 = [0., 0.4, 0.9]
clr2 = [0.9, 0.6, 0]

ax21 = fig2.add_axes([0.1, 0.5, 0.8, 0.4])
ln1, = ax21.plot(dTM, RMSE_TS, color=clr1, linewidth=2, label=expt)
ln2, = ax21.plot(dTM, RMSEp_TS, color=clr2, linewidth=2, label="Product")

ax21.set_xlim([xt1,xt2])
ax21.set_xticks(np.arange(0,xt2,5))
ax21.set_ylim([0, mmax+0.1])
ax21.grid(True)
ax21.set_xlabel('Days since 2023/03/01')
ax21.set_title('{0} & Product RMSE(m) {1}-{2}'.format(expt,rdateS,rdateE))

ax23 = plt.axes([0.7, 0.2, 0.25, 0.25])
lgd = plt.legend(handles=[ln1,ln2], loc='upper left')
ax23.axis('off')

# Plot boxplot:
dxx  = 0.1
xx0  = 1.

ax22 = fig2.add_axes([0.1, 0.08, 0.22, 0.35])
ax22.plot(xx0,md1,'.',markersize=20, color=clr1)
ax22.plot([xx0,xx0], [iqr1,iqr2], '-', linewidth=2, color=clr1)
ax22.plot([xx0-dxx,xx0+dxx], [iqr1,iqr1], '-', linewidth=2, color=clr1)
ax22.plot([xx0-dxx,xx0+dxx], [iqr2,iqr2], '-', linewidth=2, color=clr1)

xx0  = 2.
ax22.plot(xx0,mdP1,'.',markersize=20, color=clr2)
ax22.plot([xx0,xx0], [iqrP1,iqrP2], '-', linewidth=2, color=clr2)
ax22.plot([xx0-dxx,xx0+dxx], [iqrP1,iqrP1], '-', linewidth=2, color=clr2)
ax22.plot([xx0-dxx,xx0+dxx], [iqrP2,iqrP2], '-', linewidth=2, color=clr2)

ax22.set_xlim([0.5, 2.5])
ax22.set_ylim([mqmin-0.1,mqmax+0.1])
ax22.set_xticks([1,2])
ax22.grid(True)
ax22.set_xticklabels([expt,'Product'])

ax22.set_title('Inter-decile range, RMSE')


bottom_text(btx, pos=[0.02,0.01])



