"""
  Plot 2D fields from RTOFS runs
  on orthographic projections
"""
import os
import numpy as np
import sys
import importlib
import matplotlib.pyplot as plt
import datetime

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

import mod_read_hycom as mhycom
import mod_misc1 as mmisc
importlib.reload(mmisc)
import mod_utils_fig as mufig
import mod_time as mtime
import mod_utils as mutil

expt    = 'paraD5'
rdate   = '20230513'
#fld     = 'temp'
fld    = 'salin'
pthscr  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/'
pthgrid = pthscr+'hycom_fix/'
pthbin  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/' + expt + '/rtofs.' + rdate + '/'
YR     = int(rdate[0:4])
MM     = int(rdate[4:6])
DD     = int(rdate[6:8])
dnmb0  = mtime.datenum([YR,MM,DD])
yrday  = mtime.dnmb2jday(dnmb0)
hr     = -24
lrplt  = 1   # layer to plot

# F/cast output:
pthhcm = pthbin 
flhcm  = 'rtofs_glo.t00z.'.format(YR,yrday,hr)
# Renamed:
#pthhcm  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+rdate+'/'

if hr <= 0:
  flhcm = 'rtofs_glo.t00z.n-{0:02d}.archv'.format(abs(hr))
else:
  flhcm = 'rtofs_glo.t00z.f{0:02d}.archv'.format(abs(hr))

fina = pthhcm+flhcm+'.a'
finb = pthhcm+flhcm+'.b'

huge = 1.e20
rg   = 9806.
IDM, JDM, KDM = mhycom.hycom_dim(fina,finb)

ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

#idm, jdm, kdm = hycom_dim(fina, finb)

# Read S or T:
A2d, _, _, _  = mhycom.read_hycom(fina,finb,fld,rLayer=lrplt)
A2d = np.where(A2d >= 0.1*huge, np.nan, A2d)

    
from matplotlib import cm
#clrs   = cm.get_cmap('viridis',200)
clrs   = cm.get_cmap('rainbow',200)
clrs.set_bad(color=[0.7,0.7,0.7])

if fld == 'salin':
  cmpr = mutil.colormap_salin(clr_ramp=[1,0.85,1])
#  cmpr.set_bad(color=[0.2, 0.2, 0.2]) # do not use for ortho proj! 
  rmin = 30.5
  rmax = 38.5
  fld1 = 'SSS'
elif fld == 'temp' or fld == 'potT':
  cmpr = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
#  cmpr.set_bad(color=[0.2, 0.2, 0.2])  # This doesn't work with ortho proj!
#  cmpr.set_under(color=[1, 1, 1])  # ice
  rmin = -2.
  rmax = 28.
  fld1 = 'SST'
else:
  raise Exception('field {0} not defined'.format(fld))


btx  = 'plot_rtofs_scalar_orthproj.py'
sttl = 'RTOFS-DA {0} {1} {2}, {3}'.format(expt, rdate, flhcm, fld)   

importlib.reload(mutil)
mutil.plot_orthographic_proj(LON, LAT, A2d, cmpr=cmpr, sttl=sttl, \
        lon0=-50, lat0=40, rmin=rmin, rmax=rmax, btx=btx, fgnmb=1)

mutil.plot_orthographic_proj(LON, LAT, A2d, cmpr=cmpr, sttl=sttl, \
        lon0=-50, lat0=-30, rmin=rmin, rmax=rmax, btx=btx, fgnmb=2)


#plt.ion()
#fig1 = plt.figure(1,figsize=(9,9))
#plt.clf()
#ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
#im1 = ax1.pcolormesh(A2d, \
#               cmap=cmpr,\
#               vmin=rmin, \
#               vmax=rmax)
#ax1.contour(HH, [0.0], colors=[(0,0,0)], linewidths=1)
#ax1.axis('equal')
#ax1.set_xlim([2300,3200])
#ax1.set_ylim([1500,2100])

#ctl = 'paraCd S {0} m, {5} {1}/{2}/{3}:{4}'.format(z0,YR,MM,DD,hr,flhcm)
#ax1.set_title(ctl)
#
#ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
#             ax1.get_position().y0,0.02,
#             ax1.get_position().height])
#clb = plt.colorbar(im1, cax=ax2, extend='both')
#ax2.set_yticklabels(ax2.get_yticks())
#ticklabs = clb.ax.get_yticklabels()
#clb.ax.set_yticklabels(ticklabs,fontsize=10)
#clb.ax.tick_params(direction='in', length=12)
#

#btx = 'plot_rtofs_scalar_orthproj.py'
#mufig.bottom_text(btx, pos=[0.1, 0.03])






