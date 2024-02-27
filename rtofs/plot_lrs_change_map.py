"""
  Plot 2D map of layer change metric
  anls_lrs_change.py
  Analyze layer changes 
  after incremental update
  to identify problem with HYCOM blow up
  For 1 day (need to have 2 days)

  Plot 1 day or several days (for movie)
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
# for plotting 1 day: rdateS = dayE - 1
rdateS = '20230418'  # next day will be plotted, last date plotted if continue
rdateE = '20230419'  # if 1 day, this is day to plot
expt   = 'paraD'
sfx    = 'n-24'
f_figsave = False
icc0   = 0    # last saved figure, = 0 - start from 1, or skip saved

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

pthfig  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_{0}/fig/'.\
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

dnmbE = int(mtime.rdate2datenum(rdateE))
dnmbS = int(mtime.rdate2datenum(rdateS))
#dnmbS = dnmbE-1

huge = 1.e20
rg   = 9806.
icc  = icc0
for dnmb0 in range(dnmbS,dnmbE+1): 
  rdate0 = mtime.dnumb2rdate(dnmb0, ihours=False)
  dv0 = mtime.datevec(dnmb0)
  pthpara = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/wcoss2.paraD1' \
            + '/rtofs.' + rdate0 + '/'
  fhcm   = 'rtofs_glo.t00z.' + sfx + '.archv'

  fina  = pthpara + fhcm + '.a'
  finb  = pthpara + fhcm + '.b'

  if dnmb0 == dnmbS:
    IDM, JDM, KDM = mhycom.hycom_dim(fina,finb)
# Find points inside the polygon
# similar to matlab inpolygon:
    X, Y = np.meshgrid(np.arange(IDM), np.arange(JDM))
    MSKx, IMx, JMx = mmisc.inpolygon_v2(X,Y,II,JJ)  # strong relaxation region
    MSKx[np.where(HH>=-500.)] = -1
    Irgn = np.where(MSKx == 1)[0]
    DSGM = np.zeros((KDM,JDM,IDM))
    

# Actual Date of analyzed fields:
  dnmbP = mtime.dnumbTrue(dnmb0,sfx)
  dvP = mtime.datevec(dnmbP)

  dH = np.zeros((KDM,JDM,IDM))
  for kk in range(1,KDM+1):
    F,nn,mm,lmm = mhycom.read_hycom(fina,finb,'thknss',rLayer=kk)

    F = F/rg
    F[np.where(F>huge)] = np.nan
    F[np.where(F<0.001)] = 0.
    dH[kk-1,:,:] = F
# Read layer pressures:
  ZZ, ZM = mhycom.zz_zm_fromDP(dH)
#  ZZ, ZM = mhycom.get_zz_zm(fina,finb)
  ZZprf = ZZ[:,jPrf,iPrf]
  ZMprf = ZM[:,jPrf,iPrf] 
  dHprf = dH[:,jPrf,iPrf]

  if dnmb0 == dnmbS:
    dH0 = dH.copy()
    ZZprf0 = ZZprf.copy()
    ZMprf0 = ZMprf.copy()
    dHprf0 = dHprf.copy()
    print('First day {0} skipped, continue to next'.format(rdate0))
    continue

  kk1 = 14  # upper 14 layers - fixed z
#  kk1 = 0
  for kk in range(kk1,KDM):
    dh0 = np.squeeze(dH0[kk,:,:])
    dh0[np.where(MSKx<1)] = np.nan
    dh  = np.squeeze(dH[kk,:,:])
    dh[np.where(MSKx<1)] = np.nan
    if np.nanmax(dh0) < 0.1 or np.nanmax(dh) < 0.1:   # bottom
      continue
# Change in layer thicknes normalized by local depth
    dsgmd = (dh-dh0)/np.abs(HH)  
#
    rmin  = np.nanmin(dsgmd)
    jj,ii = np.where(dsgmd == rmin)
    jmin  = jj[0]
    imin  = ii[0]
    rmax = np.nanmax(dsgmd)
    jj,ii = np.where(dsgmd == rmax)
    jmax  = jj[0]
    imax  = ii[0]
# Get old/new layer thicknesses at min/max points
    dh0_min = dh0[jmin,imin]
    dh_min  = dh[jmin,imin]
    dh0_max = dh0[jmin,imin]
    dh_max  = dh[jmin,imin]

    print('Thknss change: Lr={0} min/max {1:9.6f} {2:9.6f}'.\
           format(kk+1,rmin,rmax))

    DSGM[kk,:,:] = dsgmd
 
  dH0 = dH.copy()
  ZZprf0 = ZZprf.copy()
  ZMprf0 = ZMprf.copy()
  dHprf0 = dHprf.copy()

# Plotting
  from copy import copy
  clrmp = copy(plt.cm.Spectral_r)
  clrmp.set_bad(color=[0.7,0.7,0.7])

# Plot 1 layer of max/min over all layers:
# Set lrplt = 0 for all layers
  lrplt = 0     # lr to plot
  if lrplt > 0:
    kk0 = lrplt-1  # lr index to plot
    dsgmd0 = np.squeeze(DSGM[kk0,:,:])
  else:
#    dsgmd0 = np.zeros((JDM,IDM))
    amin = np.zeros((JDM,IDM))
    amax = np.zeros((JDM,IDM))
    for kk0 in range(KDM):
      aa = np.squeeze(DSGM[kk0,:,:])
      amin = np.where(amin > aa, aa, amin)  # update min
      amax = np.where(amax < aa, aa, amax)  # update max

# Combine min/max change
    dsgmd0 = amin.copy()
    dsgmd0 = np.where(abs(dsgmd0) < amax, amax, dsgmd0) # 
    dsgmd0[np.where(MSKx<1)] = np.nan

  fgnmb = 1
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
  #  im1 = ax1.pcolormesh(rdh, vmin=-1., vmax=1., cmap=clrmp)
  im1 = ax1.pcolormesh(dsgmd0, vmin=-0.5, vmax=0.5, cmap=clrmp)
  ax1.axis('scaled')
  ax1.set_xlim([500, 650])
  ax1.set_ylim([1550,1660])
  lvl = np.arange(-5000,0,1000)
  ax1.contour(HH,[0],colors=[(0.,0.,0)])
  ax1.contour(HH,levels=lvl,linestyles='solid',colors=[(0.5,0.5,0.5)])
  # Change colormap limits:
  #im1.set_clim(-0.5,0.5)
  # Colorbar
  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='both')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  #  ctl = 'layer thkss change (dp(k+1)-dp(k))/min(dp(k),dp(k+1)) lr={0}'.format(kk+1)
  ctl = 'layer thkss change (dp(new)-dp(old))/H lr={0}\n'.format(lrplt)
  if lrplt == 0:
    ctl = 'min/max layer thkss change (dp(new)-dp(old))/H over ALL layers\n'
  ctl = ctl + '{0}/{1:02d}/{2:02d}'.format(dv0[0],dv0[1],dv0[2])
  ax1.set_title(ctl)

  btx = 'plot_lrs_change_map.py' 
  bottom_text(btx)

  icc += 1
  if f_figsave:
    if not os.path.exists(pthfig):
      print('Creating ' + pthfig)
      os.makedirs(pthfig)

    fgnm   = '{0}_lrchange_lr{1:02d}_{2:03d}.png'.format(expt,lrplt,icc)
    fpigout = pthfig + fgnm
    print('Saving date: {0}/{1:02d}/{2:02d}'.format(dv0[0],dv0[1],dv0[2]))
    print('Saving figure ---> ' + fpigout)
    plt.savefig(fpigout)




