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

import mod_read_hycom as mhycom
import mod_misc1 as mmisc
importlib.reload(mmisc)
import mod_utils_fig as mufig

# Interpolation depth
z0 = -150.   

expt    = 'rtofs_paraCd'
rdate   = '20230308'
pthscr  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/'
pthgrid = pthscr+'hycom_fix/'
pthbin  = pthscr+'data/rtofs.'+rdate+'/ocnqc_logs/profile_qc/'
YR     = int(rdate[0:4])
MM     = int(rdate[4:6])
DD     = int(rdate[6:8])
yrday  = mmisc.date_yearday(YR,MM,DD) 
hr     = 12

# F/cast output:
pthhcm = pthscr+expt+'/rtofs.'+rdate+'/'
flhcm  = 'rtofs_glo.t00z.'.format(YR,yrday,hr)
# Renamed:
#pthhcm  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+rdate+'/'

if hr <= 0:
  flhcm = 'rtofs_glo.t00z.n{0:02d}.archv'.format(abs(hr))
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
A3d  = np.array([])
for lvl in range (1,kdm+1):
  F,n1,m1,l1 = mhycom.read_hycom(fina,finb,fld,rLayer=lvl)
  F[np.where(F>huge)] = np.nan
  if A3d.size == 0:
    A3d = F.copy()
    A3d = np.expand_dims(A3d, axis=0)
  else:
    F = np.expand_dims(F, axis=0)
    A3d = np.append(A3d, F, axis=0)

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
clrs   = cm.get_cmap('rainbow',200)
clrs.set_bad(color=[0.7,0.7,0.7])


rmin = 35.
rmax = 37.

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()


ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
im1 = ax1.pcolormesh(Aint, shading='flat', \
               cmap=clrs,\
               vmin=rmin, \
               vmax=rmax)
ax1.contour(HH, [0.0], colors=[(0,0,0)], linewidths=1)
ax1.axis('equal')
ax1.set_xlim([2300,3200])
ax1.set_ylim([1500,2100])

ctl = 'paraCd S {0} m, {5} {1}/{2}/{3}:{4}'.format(z0,YR,MM,DD,hr,flhcm)
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
ax3 = plt.axes([0.1, 0.08, 0.3, 0.35])
i0  = 2675
j0  = 1675
S1  = np.squeeze(A3d[:,j0,i0])
zzm = np.squeeze(ZM[:,j0,i0])
si = Aint[j0,i0]

xl1 = 35.
xl2 = 37.
ax3.plot(S1,zzm,'.-')
ax3.plot(si,z0,'r.')
ax3.set_xlim(xl1,xl2)
ax3.set_ylim(-500,0)
ax3.plot([xl1, xl2],[z0, z0],'r--')
ax3.grid(True)
ax3.set_title('S i,j={0}, {1}'.format(i0,j0))

btx = 'plot_rtofs2z_TSmap.py'
mufig.bottom_text(btx, pos=[0.1, 0.03])






