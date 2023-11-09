"""
  Plot sections with vertical distribution of T or S
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
from copy import copy
import matplotlib.colors as colors
from matplotlib.patches import Polygon

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_cice6_utils as mc6util
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil


expt    = '003'
YR      = 2020
MM      = 6
DD      = 15
jday    = int(mtime.date2jday([YR,MM,DD]))
HR      = 12
hg      = 1.e15
isct    = [4]  # specify sections to plot
#fld     = 'potT'  # "salt" or "potT"
fld     = 'salt'  # salt or potT

XSCT = mutil.rtofs_sections()
SCTnames = list(XSCT.keys())


pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'

dnmb = mtime.jday2dnmb(YR,jday)
DV   = mtime.datevec(dnmb)
if not MM:
  MM   = DV[1]
  DD   = DV[2]


pthbin = pthrun + 'tarmom_{0}{1:02d}/'.format(YR,MM)
flmom  = 'ocnm_{0}_{1:03d}_{2}.nc'.format(YR,jday,HR)
flin   = pthbin + flmom


import mod_mom6 as mom6util
importlib.reload(mom6util)
pthgrid  = pthrun + 'INPUT/'
fgrd_mom = pthgrid + 'regional.mom6.nc'
ftopo_mom= pthgrid + 'ocean_topog.nc'
LON, LAT = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
HH       = mom6util.read_mom6depth(ftopo_mom)
Lmsk     = mom6util.read_mom6lmask(ftopo_mom)

idm, jdm, kdm = mom6util.mom_dim(flin)
ijdm          = idm*jdm

# the height of the water column = D + ssh
# hbt = sum (dH)
# HH = hbt - ssh
#hbt = mom6util.read_mom6(flin, 'hbt')
# Read layer thicknesses:
dH   = np.zeros((kdm,jdm,idm))
rfld = 'h' 
for kk in range(1,kdm+1):
  A2D = mom6util.read_mom6(flin, rfld, rLayer=kk)
  dH[kk-1,:,:] = A2D

ssh    = mom6util.read_mom6(flin, 'SSH')
ZZ, ZM = mom6util.zz_zm_fromDP(dH, ssh, f_btm=False)
#ZZ, ZM = mom6util.get_zz_zm(flin)

#j=1250
#i=1000
#mmisc.print_1col(ZZ[:,j,i])
#mmisc.print_2col(ZZ[:,j,i],dH[:,j,i])

# Read S/T fields:
A3d  = np.zeros([kdm, jdm, idm])
for kk in range (1,kdm+1):
  A2D = mom6util.read_mom6(flin, fld, rLayer=kk)
  A3d[kk-1,:,:] = A2D

# 
# Plot section
btx = 'plot_xsctTS_mom6.py'
plt.ion()

import plot_sect as psct
importlib.reload(psct)

if fld == 'salt':
  cmpr = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  cmpr.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 30.
  rmax = 40.
  fld1 = 'SctS'
elif fld == 'temp' or fld == 'potT':
  cmpr = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  cmpr.set_bad(color=[1,1,1])
  rmin = -2.
  rmax = 28.
  fld1 = 'SctT'
else:
  raise Exception('field {0} not defined'.format(fld))

smin = []
smax = []
tmin = []
tmax = []
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
  try: 
    smin = XSCT[sct_name]["smin"]
    smax = XSCT[sct_name]["smax"]
  except:
    smin = []
    smax = []
  try:
    tmin = XSCT[sct_name]["tmin"]
    tmax = XSCT[sct_name]["tmax"]
  except:
    tmin = []
    tmax = []

  if fld == 'salt' and smin:
    rmin = smin
    rmax = smax
  elif fld == 'potT' and tmax:
    rmin = tmin
    rmax = tmax
  else:
    rmin = []
    rmax = []

  # find index to print out layers:
  #  ipp,jpp,IP0,JP0 = psct.find_indx_lonlat(-37.681,-47.96,LON,LAT,xsct=xsct)
  lat0 = LAT[jx0,ix0]
  lon0 = LON[jx0,ix0]
  lon1 = LON[jx0,ix1]
  lon2 = LON[jx0,ix2]
  lat1 = LAT[jx1,ix0]
  lat2 = LAT[jx2,ix0]

  
  ss1 = '{0} {1} \n'.format(expt,flmom)
  ss1 = ss1 + '0.08 MOM6-CICE6 {0}/{1}/{2}\n'.format(YR,MM,DD)
  ss1 = ss1 + 'Fort i/i, j/j: {0}/{1}, {2}/{3}\n'.format(ix1+1,ix2+1,jx0+1,jx0+1)
  ss1 = ss1 + 'Lon, lat: {0:7.2f}/{1:7.2f}E, {2:7.2f}/{3:7.2f}N'.\
         format(lon1,lon2,lat0,lat0)

  stl  = 'EW section {2}, {0} {1}, {3}/{4}/{5}'.\
         format(Sname, sct_name, fld, YR, MM, DD)
  print('Plotting '+stl)
  ZZs, Hb, XX, dZZs = psct.plot_EWsectTS(ix1, ix2, ix0, jx0, ZZ, A3d, HH, LON, LAT,\
                      cmpr, lrstart=20, stl=stl, fgnmb=1, btx=btx, \
                      rmin=rmin, rmax=rmax, sinfo=ss1)

  # S- N section:
  ss2 = '{0} {1} \n'.format(expt,flmom)
  ss2 = ss2 + '0.08 MOM6-CICE6 {0}/{1}/{2}\n'.format(YR,MM,DD)
  ss2 = ss2 + 'Fort i/i, j/j: {0}/{1}, {2}/{3}\n'.format(ix0+1,ix0+1,jx1+1,jx2+1)
  ss2 = ss2 + 'Lon, lat: {0:7.2f}/{1:7.2f}E, {2:7.2f}/{3:7.2f}N'.\
         format(lon0,lon0,lat1,lat2)

  stl  = 'SN section {2}, {0} {1}, {3}/{4}/{5}'.\
         format(Sname, sct_name, fld, YR, MM, DD)
  print('Plotting '+stl)
  ZZn, Hbn, XXn, dZZn = psct.plot_SNsectTS(jx1, jx2, ix0, jx0, ZZ, A3d, \
                        HH, LON, LAT,\
                        cmpr, lrstart=20, stl=stl, fgnmb=2, btx=btx, \
                        rmin=rmin, rmax=rmax, sinfo=ss2)




