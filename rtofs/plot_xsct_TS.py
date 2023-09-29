"""
  Plot section vertical distr of T or S
  from RTOFS archv output

  WARNING: the code is used in automated diagnostics of para-experiments
           do not modify
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import struct
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

import mod_utils as mutil
importlib.reload(mutil)
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_time as mtime

# rdate = date of the f/cast
# incup 
#toutp = "incup"   # type of archv: incup - from 6hr incr update, fcast
#hr0   = 12

rdate0 = '20230703'  # fcast date
expt   = 'paraD5d'
sfx    = 'n-24'
isct   = [0]  # specify sections to plot
fld    = 'temp'  # salin or temp

#sct_name = 'GoM1'
XSCT = mutil.rtofs_sections()
SCTnames = list(XSCT.keys())

if max(isct) > len(SCTnames):
  raise Exception("# of requested xsections {0} > existing sections {1}".\
                   format(max(isct),len(SCTnames)))

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

#if len(rdate0) == 0:
#  yrday  = 98
#  pthhcm = '/scratch2/NCEPDEV/marine/Dan.Iredell/wcoss.20230407.run2/'
#  fhcm   = 'archv.2023_{0:03d}_12'.format(yrday)
#else:
#  pthhcm = '/scratch2/NCEPDEV/marine/Dan.Iredell/wcoss.paraB/rtofs.' + \
#           rdate0 + '/'
#  fhcm   = 'rtofs_glo.t00z.n-24.archv'
#

#pthscr = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/'
pthscr  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/{0}/'.format(expt)
#pthhcm = '/scratch2/NCEPDEV/marine/Dan.Iredell/wcoss.paraB/rtofs.' + rdate0 + '/'
#pthhcm = '/scratch2/NCEPDEV/marine/Dan.Iredell/paraD.20230501/'
#pthhcm = pthscr + 'rtofs_{0}/run_diagn/rtofs.{1}/'.format(expt,rdate0)
pthhcm = pthscr + 'rtofs.{0}/'.format(rdate0)
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

def print_1D(A,wd=8,prc=2):
  ndim1 = A.shape[0]
  for k in range (ndim1):
    print('{0}: {1:{width}.{precis}f}'.format(k+1,A[k],width=wd,precis=prc))


LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

# Read layer pressures:
dH = np.zeros((KDM,JDM,IDM))
for kk in range(1,KDM+1):
  F,nn,mm,lmm = mhycom.read_hycom(fina,finb,'thknss',rLayer=kk)

  F = F/rg
  F[np.where(F>huge)] = np.nan
  F[np.where(F<0.001)] = 1.e-6
  dH[kk-1,:,:] = F

ZZ, ZM = mhycom.zz_zm_fromDP(dH, f_btm=False)

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


btx = 'plot_xsct_TS.py'


# ==========================
# Plot selected W-E sections
import plot_sect as psct
importlib.reload(psct)


# Loop over all xsections:
# For interactive job - do 1 section at a time
if fld == 'salin':
  cmpr = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  cmpr.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 30.
  rmax = 40.
  fld1 = 'SctS'
elif fld == 'temp':
  cmpr = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  cmpr.set_bad(color=[1,1,1])
  rmin = -2.
  rmax = 28.
  fld1 = 'SctT'


nsct=len(isct)
for ii in range(nsct):
  if ii >= len(SCTnames):
    break

  isc = isct[ii]
  sct_name = SCTnames[isc]
# Specify location of the section:
# xsection in mod_utils <--- get di, dj, ix0, jx0 from XSCT
  ix1  = XSCT[sct_name]["x1"]
  ix2  = XSCT[sct_name]["x2"]
  jx0  = XSCT[sct_name]["y0"]
  jx1  = XSCT[sct_name]["y1"]
  jx2  = XSCT[sct_name]["y2"]
  ix0  = XSCT[sct_name]["x0"]
  Regn = XSCT[sct_name]["Reg"]
  Sname= XSCT[sct_name]["Name"]

  if f_figsave and (not os.path.exists(pthfig)):
    print('Creating ' + pthfig)
    os.makedirs(pthfig)

  # find index to print out layers:
  #  ipp,jpp,IP0,JP0 = psct.find_indx_lonlat(-37.681,-47.96,LON,LAT,xsct=xsct)
  lat0 = LAT[jx0,ix0]
  lon0 = LON[jx0,ix0]
  lon1 = LON[jx0,ix1]
  lon2 = LON[jx0,ix2]
  lat1 = LAT[jx1,ix0]
  lat2 = LAT[jx2,ix0]

  ss1 = '{0} {2}/{1} \n'.format(expt,fhcm,rdate0)
  ss1 = ss1 + 'RTOFS f/cast: {0}/{1}/{2}\n'.format(YR,MM,DD)
  ss1 = ss1 + 'Plotted output: {0}/{1}/{2} {3:02d}:00 UTC\n'.format(YRp,MMp,DDp,HRp) 
  ss1 = ss1 + 'Fort i/i, j/j: {0}/{1}, {2}/{3}\n'.format(ix1+1,ix2+1,jx0+1,jx0+1)
  ss1 = ss1 + 'Lon, lat: {0:7.2f}/{1:7.2f}E, {2:7.2f}/{3:7.2f}N'.\
         format(lon1,lon2,lat0,lat0)

  stl  = 'EW section {2}, {0} {1}'.format(Sname, sct_name, fld)
  print('Plotting '+stl)
  ZZs, Hb, XX, dZZs = psct.plot_EWsectTS(ix1, ix2, ix0, jx0, ZZ, A3d, HH, LON, LAT,\
                      cmpr, lrstart=42, stl=stl, fgnmb=4, btx=btx, \
                      sinfo=ss1, f_intract=f_intract)

  if f_figsave:
    fgnm   = '{5}EW_{0}_{1}_{2}_{4}_{3}.png'.\
          format(expt, rdate0, sfx, sct_name, Regn,fld1)
    fpigout = pthfig + fgnm
    print('Saving figure ---> ' + fpigout)
    plt.savefig(fpigout)

  # plt.show()
  # S- N section:
  ss2 = '{0} {2}/{1} \n'.format(expt,fhcm,rdate0)
  ss2 = ss2 + 'RTOFS f/cast: {0}/{1}/{2}\n'.format(YR,MM,DD)
  ss2 = ss2 + 'Plotted output: {0}/{1}/{2} {3:02d}:00 UTC\n'.format(YRp,MMp,DDp,HRp)
  ss2 = ss2 + 'Fort i/i, j/j: {0}/{1}, {2}/{3}\n'.format(ix0+1,ix0+1,jx1+1,jx2+1)
  ss2 = ss2 + 'Lon, lat: {0:7.2f}/{1:7.2f}E, {2:7.2f}/{3:7.2f}N'.\
         format(lon0,lon0,lat1,lat2)

  stl  = 'SN section {2}, {0} {1}'.format(Sname, sct_name, fld)
  print('Plotting '+stl)
  ZZn,Hbn,XXn,dZZn = psct.plot_SNsectTS(jx1, jx2, ix0, jx0, ZZ, A3d, HH, LON, LAT,\
                      cmpr, lrstart=42, stl=stl, fgnmb=5, btx=btx, \
                      sinfo=ss2, f_intract=f_intract)

  if f_figsave:
    fgnm   = '{5}SN_{0}_{1}_{2}_{4}_{3}.png'.\
          format(expt, rdate0, sfx, sct_name, Regn, fld1)
    fpigout = pthfig + fgnm
    print('Saving figure ---> ' + fpigout)
    plt.savefig(fpigout)


