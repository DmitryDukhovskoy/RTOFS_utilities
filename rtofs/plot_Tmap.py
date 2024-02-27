"""
  Plot T 2D map 
"""
import os
import numpy as np
import sys
import importlib
import matplotlib
import matplotlib.pyplot as plt
import datetime

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hausdorff')

import mod_read_hycom as mhycom
import mod_misc1 as mmisc
import mod_hausdorff_distance as mmhd
import mod_utils_fig as mufig
import mod_utils as mutil
import mod_rtofs as mrtofs
importlib.reload(mmisc)
import mod_gulfstream as mgulf
importlib.reload(mgulf)
import mod_time as mtime
importlib.reload(mtime)


importlib.reload(mutil)

expt    = 'paraD5'
rdate0  = '20230513'
sfx     = 'n-24'
isct    = 4     # =4 - Augulhas

REGN = mutil.rtofs_reg2Dmaps()
REGnames = list(REGN.keys())

rgn_name = REGnames[isct]
print(' Plotting T maps ' + rgn_name)


# Figure output directory:
f_intract = True  # False - no figures shown
pthfig  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/' + expt + '/fig/'
#pthscr  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/'
pthscr  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/{0}/'.format(expt) 
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

if not f_intract:
  print('Interactive mode is OFF')
  matplotlib.use('Agg')
  plt.close('all')
  plt.ioff()
else:
  plt.ion()


#pthhcm = pthscr + 'rtofs_{0}/run_diagn/rtofs.{1}/'.format(expt,rdate0)
pthhcm = pthscr + 'rtofs.{0}/'.format(rdate0)
fhcm   = 'rtofs_glo.t00z.' + sfx + '.archv'
fina   = pthhcm + fhcm + '.a'
finb   = pthhcm + fhcm + '.b'

if len(rdate0) > 0:
  YR     = int(rdate0[0:4])
  MM     = int(rdate0[4:6])
  DD     = int(rdate0[6:8])
#  yrday  = mmisc.date_yearday(YR,MM,DD)
  yrday  = mtime.rdate2jday(rdate0)
  dnmb0  = mtime.rdate2datenum(rdate0)

#
# Date of plotted fields:
if sfx == 'n-24':
  dnmbP = dnmb0-1
  hr = 0
elif sfx[0] == 'f':
  hr = int(sfx[1:])
  dnmbP = dnmb0+float(hr)/24.


dvP = mtime.datevec(dnmbP)
YRp = dvP[0]
MMp = dvP[1]
DDp = dvP[2]
HRp = dvP[3]


huge = 1.e20
rg   = 9806.
IDM, JDM, KDM = mhycom.hycom_dim(fina,finb)

ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

# Read layer pressures:
dH = np.zeros((KDM,JDM,IDM))
for kk in range(1,KDM+1):
  F,nn,mm,lmm = mhycom.read_hycom(fina,finb,'thknss',rLayer=kk)

  F = F/rg
  F[np.where(F>huge)] = np.nan
  F[np.where(F<0.001)] = 0.
  dH[kk-1,:,:] = F

ZZ, ZM = mhycom.zz_zm_fromDP(dH)
kdm = ZM.shape[0]
jdm = ZM.shape[1]
idm = ZM.shape[2]

# Read T:
fld  = 'temp'
lvl  = 1
A2d, _, _, _  = mhycom.read_hycom(fina, finb, fld, rLayer=lvl)
A2d = np.where(A2d >= 0.1*huge, np.nan, A2d)


# Function to print mouse click event coordinates
def onclick(event):
   print([event.xdata, event.ydata])


# Define region of interest:
# For the HYCOM  0.08 Global grid !!! 
xlim1 = REGN[rgn_name]["xl1"]
xlim2 = REGN[rgn_name]["xl2"]
ylim1 = REGN[rgn_name]["yl1"]
ylim2 = REGN[rgn_name]["yl2"]
Regn  = REGN[rgn_name]["Reg"]
Rname = REGN[rgn_name]["Name"]


print('=======   Start Plotting   =========')
from matplotlib import cm
#clrs   = cm.get_cmap('viridis',200)
clrs   = cm.get_cmap('rainbow',200)
#  clrs.set_under(color=[0.8,0.8,0.8])
#  clrs.set_under(color=[0.8,0.7,1])
clrs = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
#  cmpr.set_under(color=[1, 1, 1])  # ice
clrs.set_bad(color=[0.2, 0.2, 0.2])

plt.ion()

#rmin = 0.
rmin = -2.
rmax = 28.
tz0  = 12.0


fig1 = plt.figure(1,figsize=(9,9))
plt.clf()

ax1 = plt.axes([0.1, 0.2, 0.7, 0.7])

im1 = ax1.pcolormesh(A2d, shading='flat', \
               cmap=clrs,\
               vmin=rmin, \
               vmax=rmax)
#  im1.set_clim(rmin,rmax)
#  ax1.contour(HH, [0.0], colors=[(0,0,0)], linewidths=1)
#  ax1.contour(T12, [tz0], colors=[(0.,0.4,0.6)], linewidths=1)
lons = np.linspace(-180,180,73)
lats = np.linspace(-90,90,37)
ax1.contour(LON, lons, linestyles='dotted', colors=[(0.8,0.8,0.8)])
ax1.contour(LAT, lats, linestyles='dotted', colors=[(0.8,0.8,0.8)])
#
#  ax1.axis('equal')
ax1.axis('scaled')
ax1.set_xlim([xlim1,xlim2])
ax1.set_ylim([ylim1,ylim2])
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.set_xticks([])
ax1.set_yticks([])

# Actual limits:
Ylim = ax1.get_ylim()
Xlim = ax1.get_xlim()
Yl1  = int(Ylim[0])
Yl2  = int(Ylim[1])
Xl1  = int(Xlim[0])
Xl2  = int(Xlim[1])
# Put lon/lats on axis:
lon1 = LON[Yl1,Xl1]
lon2 = LON[Yl2,Xl2]
lat1 = LAT[Yl1,Xl1]
lat2 = LAT[Yl2,Xl2]

iLN1 = np.min(np.where(lons>=lon1)[0])
iLN2 = np.max(np.where(lons<=lon2)[0])
iLT1 = np.min(np.where(lats>=lat1)[0])
iLT2 = np.max(np.where(lats<=lat2)[0])
dltx = 1
dlty = 1
if iLN2-iLN1 >= 8:
  dltx = 2
if iLT2-iLT1 >= 8:
  dlty = 2

# X-axis labels
for ikk in range(iLN1,iLN2+1,dltx):
  xx0 = lons[ikk]
  yy0 = lat1   # for Mercator part of the grid, lat = const along j=fixed
  ii0, jj0 = mrtofs.find_indx_lonlat(xx0, yy0, LON, LAT)
  jj0 = Yl1-20
  xstl = '{0:3.1f}W'.format(abs(xx0))
  ax1.text(ii0, jj0, xstl, 
           fontsize=12,
           horizontalalignment='center')

# Y-axis labels
for ikk in range(iLT1,iLT2+1,dlty):
  yy0 = lats[ikk]
  xx0 = lon1
  ii0, jj0 = mrtofs.find_indx_lonlat(xx0, yy0, LON, LAT)
  ii0 = Xl1-10
  ystl = '{0:3.1f}N'.format(abs(yy0))
  if jj0 > Yl2:
    continue
  ax1.text(ii0, jj0, ystl, 
           fontsize=12,
           verticalalignment='center',
           horizontalalignment='right')

ctl = 'T level {0}, {1}/{2:02d}/{3:02d}:{4:02d}\n {6}'.\
       format(lvl, YR, MM, DD, hr, fhcm, expt)
ax1.set_title(ctl)

ax1.set_xlim(Xlim)
ax1.set_ylim(Ylim)

# Select coordinate of the region of interest:
f_setrgn = False
if f_setrgn:
# Bind the button_press_event with the onclick() method
  fig1.canvas.mpl_connect('button_press_event', onclick)

ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='both')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.tick_params(direction='in', length=12)

# Legend:
ssinf = 'RTOFS: {0}/{1:02d}/{2:02d}:{3:02d}\n'.format(YR,MM,DD,hr)
ax5 = plt.axes([0.7, 0.03, 0.2, 0.13])
ax5.text(0,0.01,ssinf)
ax5.axis('off')


btx = 'plot_Tmap.py'
mufig.bottom_text(btx, pos=[0.1, 0.03])


  


