"""
  Plot 2D fields from RTOFS runs
  Interpolate 2D field from RTOFS *archv* binary
  to a depth z0
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

import mod_time as mtime
import mod_read_hycom as mhycom
import mod_misc1 as mmisc
importlib.reload(mmisc)
import mod_utils_fig as mufig
importlib.reload(mhycom)

# Interpolation depth
z0 = -150.   

expt    = 'GOFS3.1'
rmin    = 34.79  # make it [] to get rmin/rmax 
rmax    = 37.03
rdate   = '20231030'
pthscr  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/'
pthgrid = pthscr+'hycom_fix/'
pthbin  = '{0}data/{1}/'.format(pthscr,expt)
#pthbin  = pthscr+'data/rtofs.'+rdate+'/ocnqc_logs/profile_qc/'
YR     = int(rdate[0:4])
MM     = int(rdate[4:6])
DD     = int(rdate[6:8])
yrday  = mmisc.date_yearday(YR,MM,DD) 
hr     = 0

# F/cast output:
pthhcm = pthbin 
flhcm  = 'archv.{0}{1:02d}'.format(rdate, hr)
# Renamed:
#pthhcm  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+rdate+'/'

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

# Read S or T:
A3d, _, _, _ = mhycom.read_hycom(fina,finb,fld)
A3d = np.where(A3d >= 0.1*huge, np.nan, A3d)

# 2D Arrays with S above and below z0 for interpolation
print('Searching top/btm values for {0}'.format(z0))
Atop = np.zeros((jdm,idm))*np.nan
Abtm = np.zeros((jdm,idm))*np.nan
Ztop = np.zeros((jdm,idm))*np.nan
Zbtm = np.zeros((jdm,idm))*np.nan
for kk in range(kdm-1):
  zm1 = np.squeeze(ZM[kk,:,:])
  zm2 = np.squeeze(ZM[kk+1,:,:])
  if np.nanmin(zm2) > z0:
    continue
  if np.nanmax(zm1) < z0:
    continue

  print('lvl={0}  max(z1)/min(z2) = {1:5.1f}/{2:5.1f}'.\
         format(kk+1,np.nanmax(zm1),np.nanmin(zm2)))

  [J,I] = np.where((zm1 >= z0) & (zm2 < z0)) 
  Atop[J,I] = A3d[kk,J,I]
  Abtm[J,I] = A3d[kk+1,J,I]
  Ztop[J,I] = ZM[kk,J,I]
  Zbtm[J,I] = ZM[kk+1,J,I]   

# linear interp
# Lagrange polynomial:
Aint = Atop*((z0-Zbtm)/(Ztop-Zbtm)) + \
       Abtm*((z0-Ztop)/(Zbtm-Ztop))

    
from matplotlib import cm
#clrs   = cm.get_cmap('viridis',200)
#clrs   = cm.get_cmap('rainbow',200)
#clrs.set_bad(color=[0.7,0.7,0.7])

import mod_utils as mutil
importlib.reload(mutil)

if fld == 'salin':
  cmpr = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  cmpr.set_bad(color=[0.2, 0.2, 0.2])
#  rmin = 35.
#  rmax = 37.
  fld1 = 'mapS'
elif fld == 'temp':
  cmpr = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  cmpr.set_bad(color=[0.2,0.2,0.2])
#  rmin = -2.
#  rmax = 28.
  fld1 = 'mapT'


xlim1 = 2300
xlim2 = 3200
ylim1 = 1550
ylim2 = 2100

lon1  = LON[ylim1, xlim1]
lon2  = LON[ylim1, xlim2]
lat1  = LAT[ylim1, xlim1]
lat2  = LAT[ylim2, xlim1]

Asub = Aint[ylim1:ylim2+1,xlim1:xlim2+1]
alf=2.
if not rmin:
  rmin, rmax = mutil.minmax_clrmap(Asub, pmin=alf, pmax=100-alf, cpnt=0.01)

dnmb  = mtime.rdate2datenum(rdate)
if hr < 0:
  HR = abs(hr)
  dnmbP = (dnmb-1) + float(24-HR)/24.
else:
  dnmbP = dnmb + float(hr)/24.

dvP = mtime.datevec(dnmbP)
YRp = dvP[0]
MMp = dvP[1]
DDp = dvP[2]
HRp = dvP[3]


plt.ion()
fig1 = plt.figure(1,figsize=(9,9))

plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
im1 = ax1.pcolormesh(Aint, \
               cmap=cmpr,\
               vmin=rmin, \
               vmax=rmax)
ax1.contour(HH, [0.0], colors=[(0,0,0)], linewidths=1)
ax1.axis('scaled')
ax1.set_xlim([xlim1,xlim2])
ax1.set_ylim([ylim1,ylim2])

ctl = '{0}, {1} {2:5.1f}m, {3}/{4}/{5}, {6}'.format(expt, fld, z0, YR, MM, DD, flhcm)
ax1.set_title(ctl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='both')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ss2 = '{0} {1} \n'.format(expt,flhcm,rdate)
ax6 = plt.axes([0.45, 0.4, 0.8, 0.06])
ax6.text(0.0,0.0,ss2)
ax6.axis('off')

# Vert profiles:
ax3  = plt.axes([0.1, 0.08, 0.3, 0.35])
i0   = 2675
j0   = 1675
zlim = -500.

S1  = np.squeeze(A3d[:,j0,i0])
zzm = np.squeeze(ZM[:,j0,i0])
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

btx = 'plot_rtofs2z_TSmap.py'
mufig.bottom_text(btx, pos=[0.1, 0.03])






