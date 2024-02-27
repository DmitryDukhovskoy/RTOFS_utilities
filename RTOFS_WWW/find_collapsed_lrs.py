"""
  Check for collapsed lrs
  rtofs output fields
  only deep ocean region is considered

  November 2022 
  Dmitry Dukhovskoy, NOAA/NWS/NCEP/EMC 

The first analysis is for dtg=2022111700. You are right about the background. The background archive is in my directory, rtofs.20221117/*n00*archv.a.
Remember that the rtofs day is the day when the system is run, and we allow one day for the observations to arrive, then the rtofs directory name day is one more day than the analysis day.

The first incremental update run starts from restarts rtofs.20221117/rtofs*n-06.restart.[a,b], they are not in hera
The end of the incremental update run produces rtofs.20221118/rtofs*n-24.archv.[a,b]

That is:
rtofs.20221117/*n00*archv.a   background, archive time is 2022111700
lyrprs_pre_1o4500x3298_2022111700_0000_analinc    increments
rtofs.20221118/*n-24*archv.a archive after incremental update , archive time is 2022111700
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
from copy import copy
import matplotlib
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
#from mpl_toolkits.basemap import Basemap, shiftgrid

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
import mod_read_ncoda as rncoda

import mod_utils as mutil
importlib.reload(mutil)
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_time as mtime

rdate0 = '20230303' # fcast date
expt   = 'paraD'
sfx    = 'n-24'
Regn   = 'Glb'
map_name = 'thckRthin'

# Figure output directory:
f_figsave = False
f_intract = True  # False - no figures shown
f_restart = False
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

#fhcm   = 'archv.2023_123_12'

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

LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)


from mod_read_hycom import read_grid_topo, read_hycom, \
                           read_topo, read_hycom_restart
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text

print('Processing '+fina)

def print_1D(A,wd=8,prc=2):
  ndim1 = A.shape[0]
  for k in range (ndim1):
    print('{0}: {1:{width}.{precis}f}'.format(k+1,A[k],width=wd,precis=prc))


fld  = 'thknss'
if f_restart:
  fld = 'dp'
  F,nn,mm,ll = read_hycom_restart(fina,finb,fld,IDM,JDM,rLayer=1)
else:
  F,nn,mm,ll = read_hycom(fina,finb,fld,rLayer=1)

F[np.where(F>huge)] = np.nan
F = F/rg
F[np.where(F<0.001)] = 0.

# For example 30S, 30W is  i=3199 j=1112
dH = np.zeros((ll,mm,nn))
dH[0,:,:] = F
for kk in range(2,ll+1):
  if f_restart:
    fld = 'dp'
    F,nn,mm,llm = read_hycom_restart(fina,finb,fld,IDM,JDM,rLayer=kk)
  else:
    F,nn,mm,lmm = read_hycom(fina,finb,fld,rLayer=kk)

  F = F/rg
  F[np.where(F>huge)] = np.nan
  F[np.where(F<0.001)] = 0.
  dH[kk-1,:,:] = F

ZZ, ZM = mhycom.zz_zm_fromDP(dH)
kdm = ZM.shape[0]
jdm = ZM.shape[1]
idm = ZM.shape[2]

dZZ = np.abs(np.diff(ZZ, axis=0))
# Layers to test:
k1 = 26  # if it is deeper than zUlim to avoid z-levels!
k2 = kdm-1
zUlim = -800. 

hdeep = -1300.
#Joc,Ioc   = np.where(HH < hdeep)
#Jlnd,Ilnd = np.where(HH > hdeep)
print(' Finding collapsed layers ...')
Rdp = np.zeros((jdm,idm))+1.e6
for kk in range(k1-1,k2):
  dm1 = dZZ[kk-1,:,:].copy()
  dm0 = dZZ[kk,:,:].copy()
  dp1 = dZZ[kk+1,:,:].copy()
  zzk = ZZ[kk,:,:].copy()
  zzk[np.where(HH > hdeep)] = np.nan
  zzL = np.nanmin(zzk)
  if zzL > zUlim:
    continue

  dm1[np.where(HH > hdeep)] = np.nan
  dm0[np.where(HH > hdeep)] = np.nan
  dp1[np.where(HH > hdeep)] = np.nan
  zz0 = ZZ[kk+2,:,:]  #  interface depth layer below
  dm0 = np.where((zz0 <= HH+0.5) & (zz0 >= HH-0.5),np.nan,dm0) 
  dp1 = np.where((zz0 <= HH+0.5) & (zz0 >= HH-0.5),np.nan,dp1)

  Jnan,Inan = np.where((~np.isnan(dm0)) & \
                       (~np.isnan(dp1)))
 
  if Jnan.size == 0:
    print(' Near-Bottom, all layers are ~0 or nans')
    break
  
  rm1 = dm0/(dm0+dm1)
  rp1 = dm0/(dm0+dp1)
  Rdp = np.where((rm1 < Rdp),rm1,Rdp)
  Rdp = np.where((rp1 < Rdp),rp1,Rdp)
  print('Lr {0} min/max Rdp = {1:10.8f}/{2:10.8f}'.\
         format(kk,np.nanmin(rm1),np.nanmax(rm1)))


# For balanced layers, the ratio of Lr dp2/(dp1+dp2) 
# is somewhere 0.5
# Assume 1/10 is ok, rdp0=0.1
# Smaller than that is thin/thick
Rdp[np.where(Rdp > 1.e3)] = np.nan
rdp0 = 0.1  # thin layers, as rdp-->0 - collapsed layers
 

# For easier interpretation on the figure, invert Rdp
iRdp = rdp0/Rdp # =1 thinnest layer is ~0.1 of lrs above/below
                #    + this layer
                # <1 - OK
                # >> 1 - collapsing layers

# plt.clim(a,b) = limits of colorbar
# Check lr thkn gradient:
import mod_utils as mutil
importlib.reload(mutil)

grdZ = mutil.anls_lrthkn(dZZ,lrintrf=False)

fgnmb1 = 1
print('Plotting thick-thin layers map ...')
clrmp = copy(plt.cm.afmhot_r)
clrmp.set_bad(color=[0.7,0.7,0.7])

#  plt.ion()

fig1 = plt.figure(fgnmb1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.08, 0.2, 0.8, 0.75])
im1 = plt.pcolormesh(iRdp,shading='flat',\
                     vmin=0, vmax=50, cmap=clrmp)
plt.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
ax1.axis('scaled')

ax2 = fig1.add_axes([ax1.get_position().x1+0.02, 
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='max')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)

plt.sca(ax1)
ctl = 'Ratio {0:3.2f}/min(dH(k)/[dH(k)+dH(k+/-1)]), rtofs.{1}'.\
      format(rdp0,rdate0+'/'+fhcm)
ax1.set_title(ctl)

ss1 = 'Depths > {0}m\n'.format(abs(hdeep))
ss1 = ss1 + 'Lrs = {0}/{1}, ZUlim = {2}m\n'.format(k1+1,k2+1,zUlim)
ss1 = ss1 + 'Rdp = min[dh(i)/(dh(i-1)+dh(i)),dh(i)/(dh(i)+dh(i+1))]\n'
ss1 = ss1 + 'Plotted R={0:3.2f}/Rdp\n'.format(rdp0)
ss1 = ss1 + 'Ratio (thin/thick): R(2)=1m/20m, R(4)=1m/40m, ..., R(100)=1m/1000m'
  
ax4 = plt.axes([0.08,0.08,0.7,0.11])
ax4.text(0.0,0.1,ss1)
ax4.axis('off')


btx = 'find_collapsed_lrs.py'
bottom_text(btx)

if f_figsave:
  fgtype = 'map'
  fld1   = 'thkthn'
  Regn   = 'Glb'
  fgnm   = '{0}_{1}_{2}_{3}_{4}_{5}.png'.format(fgtype,expt,rdate0,sfx,Regn,fld1)

  fpigout = pthfig + fgnm
  print('Saving figure ---> ' + fpigout)
  plt.savefig(fpigout)





