"""
  Create Land mask for gdem
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import pickle
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


# Read HYCOM grid:
expt    = 'production'
rmin    = 34.79  # make it [] to get rmin/rmax 
rmax    = 37.03
rdate   = '20231031'
pthscr  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/'
pthgrid = pthscr+'hycom_fix/'
pthbin  = '{0}data/{1}/rtofs.{2}/'.format(pthscr,expt,rdate)
#pthbin  = pthscr+'data/rtofs.'+rdate+'/ocnqc_logs/profile_qc/'
YR     = int(rdate[0:4])
MM     = int(rdate[4:6])
DD     = int(rdate[6:8])
yrday  = mmisc.date_yearday(YR,MM,DD)
hr     = -24

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

#fld = 'temp'
fld = 'salin'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
Lmsk = np.where(HH >= 0., 0., 1.)
JL, IL = np.where(Lmsk < 0.1)
XL = LON[JL,IL]
YL = LAT[JL,IL]

LONG, LATG, ZMG = mgdem.gdem_grid()
NI = LONG.shape[0]
NJ = LATG.shape[0]
NK = ZMG.shape[0]
LONG = np.where(LONG > 180., LONG-360., LONG)

def find_indx(XX, x1):
  dd = np.abs(XX-x1)
  ix = np.where(dd==np.min(dd))[0][0]
  return ix


LmskG           = np.zeros((NJ,NI)) + 1
isp             = max(np.where(LATG < np.min(LAT))[0])
LmskG[:isp+1,:] = 0

import mod_misc1 as mmisc1
npnts = len(XL)
for ii in range(npnts):
  x0 = XL[ii]
  y0 = YL[ii]
  jx = find_indx(LATG, y0)
  yg = LATG[jx]
  latg = LONG.copy()*0.0 + yg
  dd = mmisc1.dist_sphcrd(latg, LONG, y0, x0)
  ix = np.where(dd==np.min(dd))[0][0]
  xg = LONG[ix]
  LmskG[jx,ix] = 0

  if ii%100000 == 0:
    print(' {0:8.2f}% done ...'.format(float(ii)/float(npnts)*100.))
    print(' HYCOM lon/lat = {0:6.2f}/{1:6.2f}, GDEM lon/lat = {2:6.2f}/{3:6.2f}'.\
          format(x0, y0, xg, yg))

LmskG[:,-1] = LmskG[:,0]

pthout   =  '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/GDEM/'
gdem_out = pthout + 'GDEM_sealand_mask.pkl'
print('Saving --> ' + gdem_out)
with open(gdem_out,'wb') as fid:
  pickle.dump(LmskG, fid)



f_plt = False
if f_plt:
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

