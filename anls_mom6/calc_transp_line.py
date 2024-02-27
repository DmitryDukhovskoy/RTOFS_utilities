"""
  Calculate volume transport
  across sections/straits
  connected with a straight line
  EW or SN 

  zigzagging line - needs different code

  Note variables are on staggered grid in the output files
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pdb
import importlib
import time
#import struct
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/mom6_utils')

import mod_time as mtime
from mod_utils_fig import bottom_text

import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)

expt    = '003'
hg      = 1.e15
sctnm   = 'Fram79'
dnmb1 = mtime.datenum([2021,1,1])
dnmb2 = mtime.datenum([2021,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)


pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'
pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/' + \
          'MOM6_CICE6/expt{0}/'.format(expt)

STR = mom6vld.ocean_straits()
i1  = STR[sctnm]["xl1"]
i2  = STR[sctnm]["xl2"]
j1  = STR[sctnm]["yl1"]
j2  = STR[sctnm]["yl2"]
ni  = i2+1
nj  = j2+1

if j1 == j1:
  xsct_EW = True
else:
  xsct_EW = False


import mod_misc1 as mmisc
import mod_mom6 as mom6util
importlib.reload(mom6util)

dday=1
TPLT = mom6util.create_time_array(dnmb1, dnmb2, dday, date_mat=True)
nrec = len(TPLT)
for irec in range(nrec):
  rtimeS = time.time()
  DV   = TPLT[irec][1]
  YR   = DV[0]
  MM   = DV[1]
  DD   = DV[2]
  jday = int(TPLT[irec][2])
  dnmb = TPLT[irec][0]  
#  jday    = int(mtime.date2jday([YR,MM,DD]))
#  dnmb = mtime.jday2dnmb(YR,jday)
  HR   = 12

  pthbin = pthrun + 'tarmom_{0}{1:02d}/'.format(YR,MM)
  flmom  = 'ocnm_{0}_{1:03d}_{2}.nc'.format(YR,jday,HR)
  flin   = pthbin + flmom

  pthgrid   = pthrun + 'INPUT/'
  fgrd_mom  = pthgrid + 'regional.mom6.nc'
  ftopo_mom = pthgrid + 'ocean_topog.nc'
  if irec == 0:
    LON, LAT  = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
    HH        = mom6util.read_mom6depth(ftopo_mom)
    Lmsk      = mom6util.read_mom6lmask(ftopo_mom)
    idm, jdm, kdm = mom6util.mom_dim(flin)
    ijdm          = idm*jdm

  # Read layer thicknesses:
  rfld = 'h'
  dH = mom6util.read_mom6(flin, rfld, finfo=False)

  ssh    = mom6util.read_mom6(flin, 'SSH', finfo=False)
  ZZ, ZM = mom6util.zz_zm_fromDP(dH, ssh, f_btm=False, finfo=False)
  #ZZ, ZM = mom6util.get_zz_zm(flin)

  # Calc segment lengths:
  if irec==0:
    DX, DY = mom6util.dx_dy(LON,LAT)

  #j=1250
  #i=1000
  #mmisc.print_1col(ZZ[:,j,i])
  #mmisc.print_2col(ZZ[:,j,i],dH[:,j,i])

  # Read U/V only normal velocity needed:
  if xsct_EW:
    fld='v'
    XX = LON[j1,i1:ni].squeeze()
    YY = LAT[j1,i1:ni].squeeze()
    if XX[-1] < XX[0] and XX[-1] < 0.:
      XX = np.where(XX<0, XX+360., XX)
    sttl = 'V'
    LSgm = DX[j1,i1:ni].squeeze()
  else:
    fld='u'
    XX = LAT[j1:nj,i1].squeeze()
    YY = LAT[j1:nj,i1].squeeze()
    LSgm = DY[j1:nj,i1].squeeze()

    sttl = 'U'

  A3d = mom6util.read_mom6(flin, fld, finfo=False)

  # Interpolated onto h-pnt:
  if xsct_EW:
  # Collocate V components to h pnt
    Uj   = A3d[:,j1,i1:ni].squeeze()
    Ujm1 = A3d[:,j1-1,i1:ni].squeeze() 
    Uj   = np.where(np.isnan(Uj), 0., Uj)
    Ujm1 = np.where(np.isnan(Ujm1), 0., Ujm1)
    Uath = 0.5*(Uj + Ujm1)   # collocated with h
  else:
  # Collocate U components to h pnt
    Ui   = A3d[:,j1:nj,i1].squeeze()
    Uim1 = A3d[:,j1:nj,i1-1].squeeze()
    Ui   = np.where(np.isnan(Ui), 0., Ui)
    Uim1 = np.where(np.isnan(Uim1), 0., Uim1)
    Uath = 0.5*(Ui + Uim1)   # collocated with h

  dH2d  = dH[:,j1:nj,i1:ni].squeeze()
  Hb    = HH[j1:nj,i1:ni].squeeze()
  ZZ2d  = ZZ[:,j1:nj,i1:ni].squeeze()

  dZ = abs(np.diff(ZZ2d, axis=0))
  # 
  # Plot section
  f_plt = False
  if f_plt:
    btx = 'calc_transp_line.py'
    plt.ion()

    stl = ('0.08 MOM6-CICE6-{0} {1} {2}, {3}/{4}/{5}'.\
            format(expt, sttl, sctnm, YR, MM, DD))
    importlib.reload(mom6vld)
    mom6vld.plot_xsect(XX, Hb, ZZ2d, Uath, HH, stl=stl, ijsct=[i1,i2,j1,j2], btx=btx)


  # Vol transport:
  VFlx = mom6util.vol_transp_2Dsection(LSgm, ZZ2d, Uath)

  if irec == 0:
    VTRP = mom6vld.TRANSP(dnmb, VFlx, XX, YY, Hb)
  else:
    VTRP.add_array(dnmb, VFlx)

  rtime0 = time.time()
  dtime  = (rtime0-rtimeS)/60.
  print('{0:4d}: {1}/{2:02d}/{3:02d} {4} VolTrt = {5:7.3f} Sv'.\
        format(irec, YR, MM, DD, sctnm, np.nansum(VFlx)*1.e-6))
  print('  1 day processed, {0:8.3f} min ...'.format(dtime))

# Save:
floutp = 'mom6-{4}_volTrt_{0}{1:02d}-{2}{3:02d}_{5}.pkl'.\
         format(dv1[0], dv1[1], dv2[0], dv2[1], expt, sctnm)
ffout = pthoutp + floutp
print('Saving to '+ffout)
with open(ffout,'wb') as fid:
  pickle.dump(VTRP,fid)

print(' All done ' )

f_chck=False
if f_chck:
  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
  ax1.plot(XX,VFlx*1e-6)

