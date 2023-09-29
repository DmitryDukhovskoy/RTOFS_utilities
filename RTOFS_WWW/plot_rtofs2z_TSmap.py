"""
  Plot 2D fields from RTOFS runs
  Interpolate 2D field from RTOFS *archv* binary
  to a depth z0
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

import mod_read_hycom as mhycom
import mod_misc1 as mmisc
importlib.reload(mmisc)
import mod_utils_fig as mufig

import mod_utils as mutil
importlib.reload(mutil)
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_time as mtime

# Interpolation depth
z0 = 0.

expt    = 'paraD'
rdate0  = '20230419'
fld     = 'salin'  # temp or salin
sfx     = 'n-24'
isct    = [0]      # region to plot
fld     = 'salin'  # salin or temp

REGN = mutil.rtofs_reg2Dmaps()
REGnames = list(REGN.keys())

# Figure output directory:
f_figsave = False
f_intract = True  # False - no figures shown
pthfig = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/' + expt + '/fig/'

if not f_intract:
  print('Interactive mode is OFF')
  matplotlib.use('Agg')
  plt.close('all')
  plt.ioff()
else:
  plt.ion()

pthscr = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/'
#pthhcm = '/scratch2/NCEPDEV/marine/Dan.Iredell/wcoss.paraB/rtofs.' + rdate0 + '/'
#pthhcm = '/scratch2/NCEPDEV/marine/Dan.Iredell/paraD.20230501/'
pthhcm = pthscr + 'rtofs_{0}/run_diagn/rtofs.{1}/'.format(expt,rdate0)
pthgrid= '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
fhcm   = 'rtofs_glo.t00z.' + sfx + '.archv'

fina = pthhcm + fhcm + '.a'
finb = pthhcm + fhcm + '.b'

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


#IDM  = 4500
#JDM  = 3298
huge = 1.e20
rg   = 9806.

get_topo = True
ftopo = 'regional.depth'
fgrid = 'regional.grid'


import mod_read_hycom as mhycom
importlib.reload(mhycom)
from mod_utils_fig import bottom_text

print('Processing '+fina)
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

# interpolated depth should be same or deeper than layer 1
zlr1 = np.nanmin(ZM[0,:,:])
if z0 > zlr1:
  print('    ')
  print('Requested interpolate z={0} shallower than 1st lr={1}'.format(z0,zlr1))
  print('Adjusting z0  ---> {0}'.format(zlr1))
  print('    ')
  z0 = zlr1

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

import plot_sect as psct
importlib.reload(psct)

# Loop over all xsections:
# For interactive job - do 1 section at a time
if fld == 'salin':
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

nrgn = len(isct)
for ii in range(nrgn):
  if ii >= len(REGnames):
    break

  isc = isct[ii]
  rgn_name = REGnames[isc]
  print(' Plotting T/S maps ' + fld + '  ' + rgn_name)

  xlim1 = REGN[rgn_name]["xl1"]
  xlim2 = REGN[rgn_name]["xl2"]
  ylim1 = REGN[rgn_name]["yl1"]
  ylim2 = REGN[rgn_name]["yl2"]
  Regn  = REGN[rgn_name]["Reg"]
  Rname = REGN[rgn_name]["Name"]
  prf1  = REGN[rgn_name]["ij1"]
  prf2  = REGN[rgn_name]["ij2"]
  prf3  = REGN[rgn_name]["ij3"]

  PRF = np.array(prf1)
#  PRF = np.expand_dims(PRF, axis=(0))
  PRF = np.append(PRF,np.array(prf2))
  PRF = np.append(PRF,np.array(prf3))
  PRF = np.reshape(PRF,(3,2))
  nprf = PRF.shape[0]

  Asub = Aint[ylim1:ylim2+1,xlim1:xlim2+1]
  alf=2.
  rmin = np.nanpercentile(Asub,alf)
  rmin = round(rmin,2)
  rmax = np.nanpercentile(Asub,100.-alf)
  rmax = round(rmax,2)


#  plt.ion()
  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()

  ax1 = plt.axes([0.08, 0.55, 0.8, 0.4])
  im1 = ax1.pcolormesh(Aint, shading='flat', \
                 cmap=cmpr,\
                 vmin=rmin, \
                 vmax=rmax)
  ax1.contour(HH, [0.0], colors=[(1,1,1)], linewidths=1)

# Show profile locations
  for ip in range(nprf):
    xx = PRF[ip,0]
    yy = PRF[ip,1]
    ax1.plot(xx+5,yy+5,'.',color=(0,0,0))
    ax1.text(xx+5,yy+5,'{0}'.format(ip+1))


  ax1.axis('scaled')
  ax1.set_xlim([xlim1,xlim2])
  ax1.set_ylim([ylim1,ylim2])

  ctl = '{6} {0}m, {7}/{5} {1}/{2}/{3}:{4}'.format(z0,YR,MM,DD,hr,fhcm,fld,expt)
  ax1.set_title(ctl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='both')
  ax2.set_yticklabels(ax2.get_yticks())
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
#  ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  ss2 = '{0} {2}/{1} \n'.format(expt,fhcm,rdate0)
  ss2 = ss2 + 'RTOFS f/cast: {0}/{1}/{2} '.format(YR,MM,DD)
  ss2 = ss2 + 'Plotted output: {0}/{1}/{2} {3:02d}:00 UTC\n'.format(YRp,MMp,DDp,HRp)
  ax6 = plt.axes([0.06, 0.47, 0.8, 0.06])
  ax6.text(0.0,0.0,ss2)
  ax6.axis('off')

  # Vert profiles:
  hax = 0.35
  wax = 0.22
  zlim  = -500.

  for ip in range(nprf):
    if ip == 0:
      ax3 = plt.axes([0.08, 0.08, wax, hax])
      ax0 = ax3
    elif ip == 1:
      ax4 = plt.axes([0.40, 0.08, wax, hax])
      ax0 = ax4
    else:
      ax5 = plt.axes([0.72, 0.08, wax, hax])
      ax0 = ax5

    i0  = PRF[ip,0]
    j0  = PRF[ip,1]
    S1  = np.squeeze(A3d[:,j0,i0])
    zzm = np.squeeze(ZM[:,j0,i0])
    si = Aint[j0,i0]

    if min(zzm) <= zlim:
      izlim = np.min(np.where(zzm <= zlim )) + 1
    else:
      izlim = np.max(np.where(zzm == np.nanmin(zzm)))

    xl1 = np.floor(np.nanmin(S1[:izlim])*10.)/10.
    xl2 = np.ceil(np.nanmax(S1[:izlim])*10.)/10.
    stl = 'Prof {2} i,j= {0}, {1}'.format(i0,j0,ip+1)
    psct.plot_profile(ax0, S1, zzm, z0, stl=stl, zlim=zlim, xl1=xl1, xl2=xl2)

  btx = 'plot_rtofs2z_TSmap.py'
  mufig.bottom_text(btx, pos=[0.1, 0.02])

  if f_figsave:
    fgnm   = '{0}_{1}_{2}_{3}_{4}_{5}.png'.\
          format(fld1, expt, rdate0, sfx, Regn, rgn_name)
    fpigout = pthfig + fgnm
    print('Saving figure ---> ' + fpigout)
    plt.savefig(fpigout)



