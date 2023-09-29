"""
  Analyze layer changes 
  after incremental update
  to identify problem with HYCOM blow up
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
rdateS = '20230403'  # day to start
rdateE = '20230430'  # day to end
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

huge = 1.e20
rg   = 9806.

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
    RMMX = mutil.rmmx(KDM)

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
    continue

# Compute change
#  dltH = dH-dH0
  RMIN   = np.zeros((KDM,1))
  RMAX   = np.zeros((KDM,1))
  IMIN   = np.zeros((KDM,1))
  JMIN   = np.zeros((KDM,1))
  IMAX   = np.zeros((KDM,1))
  JMAX   = np.zeros((KDM,1))
  DH0MIN = np.zeros((KDM,1))
  DH0MAX = np.zeros((KDM,1))
  DHMIN  = np.zeros((KDM,1))
  DHMAX  = np.zeros((KDM,1))
  kk1 = 14  # upper 14 layers - fixed z
#  kk1 = 0
  sdp0 = np.squeeze(np.nansum(dH0[0:kk1-1,:,:],axis=0))
  sdp0[np.where(MSKx<1)] = np.nan
  sdp = np.squeeze(np.nansum(dH[0:kk1-1,:,:],axis=0))
  sdp[np.where(MSKx<1)] = np.nan
  for kk in range(kk1,KDM):
    dh0 = np.squeeze(dH0[kk,:,:])
    dh0[np.where(MSKx<1)] = np.nan
    dh  = np.squeeze(dH[kk,:,:])
    dh[np.where(MSKx<1)] = np.nan
    if np.nanmax(dh0) < 0.1 or np.nanmax(dh) < 0.1:   # bottom
      continue
    sdp0 = sdp0 + dh0
    sdp = sdp + dh    
#
# Change in sigma-depths, i.e. intgr(dp)/total_depth:
    sgmd0 = sdp0/np.abs(HH)
    sgmd  = sdp/np.abs(HH)
    dsgmd = (sgmd/sgmd0 - 1.)  # relative sigma-depth change of interface

#    ddh = np.squeeze(dltH[kk,:,:])
#    ddh[np.where(MSKx<1)] = np.nan
# Normalize by min of previous/updated lr thickness:
#    dhdh0 = dh0.copy()
# Take care about appearing/vanishing bottom layers from 0 to ...
#    dhmin = 0.1
#    dhdh0 = np.where(dh0 < dhmin, dh, dh0)
#    dhdh0 = np.where(dhdh0 < dhmin, dhmin, dhdh0)
#    rdh = ddh/dhdh0 

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

    RMIN[kk,0] = rmin
    IMIN[kk,0] = imin
    JMIN[kk,0] = jmin
    RMAX[kk,0] = rmax
    IMAX[kk,0] = imax
    JMAX[kk,0] = jmax
    DH0MIN[kk,0] = dh0_min
    DH0MAX[kk,0] = dh0_max
    DHMIN[kk,0]  = dh_min
    DHMAX[kk,0]  = dh_max
    print('Thknss change: Lr={0} min/max {1:9.6f} {2:9.6f}'.\
           format(kk+1,rmin,rmax))

  RMMX.add_data(RMIN,RMAX,IMIN,JMIN,IMAX,JMAX,DH0MIN,DH0MAX,DHMIN,DHMAX,dnmb0)
 
  dH0 = dH.copy()
  ZZprf0 = ZZprf.copy()
  ZMprf0 = ZMprf.copy()
  dHprf0 = dHprf.copy()


f_plt=False
if f_plt:
  from copy import copy
  clrmp = copy(plt.cm.Spectral_r)
  clrmp.set_bad(color=[0.7,0.7,0.7])

  fgnmb = 1
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
#  im1 = ax1.pcolormesh(rdh, vmin=-1., vmax=1., cmap=clrmp)
  im1 = ax1.pcolormesh(dsgmd, vmin=-1., vmax=1., cmap=clrmp)
  ax1.axis('scaled')
  ax1.set_xlim([450, 700])
  ax1.set_ylim([1525,1700])
  lvl = np.arange(-5000,0,1000)
  ax1.contour(HH,[0],colors=[(0.,0.,0)])
  ax1.contour(HH,levels=lvl,linestyles='solid',colors=[(0.5,0.5,0.5)])

# Colorbar
  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='both')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

#  ctl = 'layer thkss change (dp(k+1)-dp(k))/min(dp(k),dp(k+1)) lr={0}'.format(kk+1)
  ctl = 'layer thkss change (intgr{dp(z)}/intgr{dp(z)old})-1) lr={0}\n'.format(kk+1)
  clt = ctl + '{0}/{1:02d}/{2:02d}'.format(dv0[0],dv0[1],dv0[2])
  ax1.set_title(ctl)

  btx = 'anls_lrs_change.py' 
  bottom_text(btx)

if f_save:
  fldump = pthdump + '/lrthknss_anls.pkl'
  print('Saving analyzed data to '+fldump)
  with open(fldump,'wb') as fid:
    pickle.dump(RMMX,fid)

  





