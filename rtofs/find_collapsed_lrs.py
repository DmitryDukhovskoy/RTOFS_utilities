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
import pdb
import netCDF4
import importlib
import struct
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
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
#importlib.reload(rncoda)

#rdate   = '20211103' # no NCODA 
#rdate   = '20221103' # 1 year with NCODA
#pthbin  = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/FOR_DMITRY/rtofs.'+rdate+'/'
#rdate   = '20220616'  # From GOFS restart
#rdate   = '20220617'  # after 1 NCODA + incr. update
#rdate   = '20220618'  # after 2 NCODA
#rdate   = '20220623'
#pthbin  = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/rtofs_expts/rtofs_para8.1/hycom/rtofs.'\
#           + rdate+'/'
#rdate   = '20221117'
rdate   = '20220617'
hr      = 18
#pth1    = '/scratch2/NCEPDEV/marine/Dan.Iredell/para8d/rtofs.'
#pth1    = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/rtofs_expts/'+\
#           'rtofs_para8.1b/rtofs.'
#pth1    = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/'+\
#           'rtofs_para7b/hycom/ncoda_archv_inc'
pth1    = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/'+\
           'rtofs_para7b/hycom/rtofs.'

date_incup = rncoda.date_type_output(rdate,"incup")
date_fcst  = rncoda.date_type_output(rdate,"fcast")
date_outp  = date_incup

if pth1[-15:] == 'ncoda_archv_inc':
  pthbin = pth1+'/'
else:
  pthbin  = pth1+rdate+'/'

pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

#flnm = 'rtofs_glo.t00z.n-24.archv'    # 6hr ingested NCODA incremental updates
#flnm = 'restart_r2022111600_930'      # NRL restart with no NCODA 
#flnm = 'rtofs_glo.t00z.n00.archv'    # 24-hr hindcast or background for next day
#if rdate == '20220616':
#  flnm = 'archv.2022061600'
#jday = rncoda.rdate2julian(rdate)
jday = rncoda.rdate2julian(date_outp)
yr, mo, mday, hr = rncoda.parse_rdate(date_outp)
# Background field:
#flnm = 'archv.{0}_{1:03d}_00'.format(yr,jday)   # background field
# Increment:
#flnm = 'archv_1_inc.2022_{0:03d}_00'.format(jday)  # added increments  ncoda_archv_inc
#
# Forecast:
#flnm = 'archv.{0}_{1:03d}_{2:02d}'.format(yr,jday,hr)
flnm = 'rtofs_glo.t00z.n-24.archv'   # renamed after 6hr incr update

fina = pthbin+flnm+'.a'
finb = pthbin+flnm+'.b'

IDM = 4500
JDM = 3298

f_restart = flnm[0:7]=='restart'

get_topo = True
ftopo = 'regional.depth'
fgrid = 'regional.grid'


import mod_read_hycom
#importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo, read_hycom, \
                           read_topo, read_hycom_restart
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text

print('Processing '+fina)

def print_1D(A,wd=8,prc=2):
  ndim1 = A.shape[0]
  for k in range (ndim1):
    print('{0}: {1:{width}.{precis}f}'.format(k+1,A[k],width=wd,precis=prc))


huge = 1.e20
rg   = 9806.
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

# ================================
def find_indx_lonlat(x0,y0):
  """
  Find closest grid point to lon/lat coordinate
  """
  if x0 > 180.:
    x0 = x0-360.
  
  dmm = np.sqrt((LON-x0)**2+(LAT-y0)**2)
  jj0, ii0 = np.where(dmm == np.min(dmm))

  return ii0[0], jj0[0]

# ================================
ZZ, ZM = zz_zm_fromDP(dH)
kdm = ZM.shape[0]
jdm = ZM.shape[1]
idm = ZM.shape[2]

if get_topo:
#  HH = read_topo(pthgrid,ftopo,nn,mm)
  LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)
  get_topo = False  


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

f_pltrdp = True
if f_pltrdp:
  fgnmb1 = 1
  print('Plotting thick-thin layers map ...')
  clrmp = copy(plt.cm.afmhot_r)
  clrmp.set_bad(color=[0.7,0.7,0.7])
  plt.ion()
  fig1 = plt.figure(fgnmb1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im1 = plt.pcolormesh(iRdp,shading='flat',\
                       vmin=0, vmax=10, cmap=clrmp)
  plt.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
 
  ax2 = fig1.add_axes([ax1.get_position().x1+0.02, 
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='max')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(ticklabs,fontsize=10)

  plt.sca(ax1)
  ctl = 'Inverse min(lr_thkns/dH) wrt {0:3.2f}, rtofs.{1}'.\
        format(rdp0,rdate+'/'+flnm)
  ax1.set_title(ctl)
  ss1 = 'Depths > {0}m'.format(abs(hdeep))
  ss2 = 'Lrs = {0}/{1}'.format(k1+1,k2+1)
  ss5 = 'ZUlim = {0}m'.format(zUlim)
  ss3 = 'Rdp = min[dh(i)/(dh(i-1)+dh(i)),dh(i)/(dh(i)+dh(i+1))]'
  ss4 = '{0:3.2f}/Rdp'.format(rdp0)

  txt = ss1+'\n'+ss2+'\n'+ss5+'\n'+ss3+'\n'+ss4
  ax1.text(10,2700,txt) 
 
  btx = 'find_collapsed_lrs.py'
  bottom_text(btx)


# ==========================
# Plot selected W-E profiles  
import plot_sect as psct
importlib.reload(psct)

itst = 2412
jtst = 1798
f_plt = True
if f_plt:
#  xsct = "SeaJpn"
  xsct = "Test1"
#  xsct = "Test2"
#  xsct = "GoM1"
#  xsct = "NAtl1"
# find index to print out layers:
  ipp,jpp,IP0,JP0 = psct.find_indx_lonlat(-37.681,-47.96,LON,LAT,xsct=xsct)

  if xsct == "GoM1":
    jp0 = 72  # GoM corresponds the itest, jtest in ncoda_archv_inc
  elif xsct == "SeaJpn":
    jp0 = 56
  elif xsct == "Test1":
    jp0 = 102
  elif xsct == "NAtl1":
    jp0 = 50
  elif xsct == "Test2":
    jp0 = 60

#  drfnm = "rtofs."+rdate+"/"+flnm
  drfnm = pthbin[-16:]+flnm
  ZZs, Hb, XX, dZZs = psct.plot_Jsection(xsct,ZZ,HH,LON,LAT,drfnm,sct_show=True,\
                      fgnmb=5,dcntr=1, zstart=-1000., jlrs=jp0, btx=btx)
#  print_1D(dZZs[:,jp0])
  psct.print_2D(ZZs[:,jp0],dZZs[:,jp0],kend=kdm)
 






