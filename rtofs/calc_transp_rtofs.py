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
sys.path.append('/home/Dmitry.Dukhovskoy/python/anls_mom6')

import mod_read_hycom as mhycom
import mod_time as mtime
from mod_utils_fig import bottom_text

import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)

f_cont = False      # True - load saved and start from last saved
expt   = 'product'
sfx    = 'n-24'
YR1    = 2021
sctnm   = 'Fram79'
huge   = 1.e15
rg     = 9806.

dnmb1 = mtime.datenum([2021,1,1])
dnmb2 = mtime.datenum([2021,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)

pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/wcoss2.prod/'
pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/RTOFS_production/'
floutp = f"rtofs-{expt}_volTrt_{dv1[0]}" + \
         f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
ffout = pthoutp + floutp

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

# Time array
dday = 7
irc  = -1
nrec = 48
TPLT = np.zeros(nrec, dtype=[('dnmb', float),
                              ('date', int, (4,)),
                              ('yrday', float)])
for imo in range(1, 13):
  for mday in range(2, 29, dday):
    irc += 1
    dnmb = mtime.datenum([YR1,imo,mday])
    DV   = mtime.datevec(dnmb)
    _, jday = mtime.dnmb2jday(dnmb)
    TPLT['dnmb'][irc]   = dnmb
    TPLT['yrday'][irc]  = jday
    TPLT['date'][irc,:] = DV[0:4]

dnmbL = 0
if f_cont:
  print('Start from last saved record')
  print('Loading ' + ffout)
  try:
    with open(ffout,'rb') as fid:
      VTRP = pickle.load(fid)

    TM    = VTRP.TM
    nlast = len(TM)
    dnmbL = TM[-1]
    dvL   = mtime.datevec(dnmbL)
    print(f"{nlast} saved records, last {dvL[0]}/{dvL[1]}/{dvL[2]}")
  except:
#    raise Exception('Saved output not found ...')
    print(f"Saved output not found {ffout}")
    print("Starting from 0")


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

  print(f"Processing rec {irec+1} {YR}/{MM}/{DD}")
  pthbin  = pthrun + 'tarmom_{0}{1:02d}/'.format(YR,MM)
  pthbin  = pthrun + f"rtofs.{YR}{MM:02d}{DD:02d}/"
  flhycom = f"rtofs_glo.t00z.{sfx}.archv"
  fina    = pthbin + flhycom + '.a'
  finb    = pthbin + flhycom + '.b'

  pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
  ftopo = 'regional.depth'
  fgrid = 'regional.grid'
  if irec == 0:
    LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

  # Calc segment lengths:
  if irec==0:
    DX, DY = mhycom.dx_dy(LON,LAT)

  if dnmb <= dnmbL:
    print(f"Skipping rec {irec+1} {YR}/{MM}/{DD} --->")
    continue

  # Read layer thicknesses:
  dH, idm, jdm, kdm = mhycom.read_hycom(fina,finb,'thknss')
  dH = dH/rg
  dH = np.where(dH>huge, np.nan, dH)
  dH = np.where(dH<0.001, 0., dH)

  ZZ, ZM = mhycom.zz_zm_fromDP(dH, f_btm=False)

  if xsct_EW:
    fld = 'v-vel.'
# Read barotropic V
# Only for archv output !
    F, _, _, _ = mhycom.read_hycom(fina,finb,'v_btrop')
    F = np.where(F>=huge, 0., F)
    Ubtrop = F.copy()
  else:
    fld = 'u-vel.'
# Read barotropic U
    F, _, _, _ = mhycom.read_hycom(fina,finb,'u_btrop')
    F = np.where(F>=huge, 0., F)
    Ubtrop = F.copy()

  # Read U/V only normal velocity needed:
  if xsct_EW:
    XX = LON[j1,i1:ni].squeeze()
    YY = LAT[j1,i1:ni].squeeze()
    if XX[-1] < XX[0] and XX[-1] < 0.:
      XX = np.where(XX<0, XX+360., XX)
    sttl = 'V'
    LSgm = DX[j1,i1:ni].squeeze()
  else:
    XX = LAT[j1:nj,i1].squeeze()
    YY = LAT[j1:nj,i1].squeeze()
    LSgm = DY[j1:nj,i1].squeeze()

    sttl = 'U'

  A3d, _, _, _ = mhycom.read_hycom(fina, finb, fld)
  A3d = np.where(A3d >= huge, np.nan, A3d)
  # Interpolated onto h-pnt:
  if xsct_EW:
# Add barotropic U:
    for kk in range(kdm):
      A3d[kk,:,:] = A3d[kk,:,:] + Ubtrop
# Collocate V components to h pnt
    Uj   = A3d[:,j1,i1:ni].squeeze()
    Ujm1 = A3d[:,j1-1,i1:ni].squeeze()
    Uj   = np.where(np.isnan(Uj), 0., Uj)
    Ujm1 = np.where(np.isnan(Ujm1), 0., Ujm1)
    A2d = 0.5*(Uj + Ujm1)   # collocated with h
  else:
# Add barotropic U:
    for kk in range(kdm):
      A3d[kk,:,:] = A3d[kk,:,:] + Ubtrop
# Collocate U components to h pnt
    Ui   = A3d[:,j1:nj,i1].squeeze()
    Uip1 = A3d[:,j1:nj,i1+1].squeeze()
    Ui   = np.where(np.isnan(Ui), 0., Ui)
    Uip1 = np.where(np.isnan(Uip1), 0., Uip1)
    A2d = 0.5*(Ui + Uip1)   # collocated with h

  dH2d  = dH[:,j1:nj,i1:ni].squeeze()
  Hb    = HH[j1:nj,i1:ni].squeeze()
  ZZ2d  = ZZ[:,j1:nj,i1:ni].squeeze().copy()
  ZM2d  = ZM[:,j1:nj,i1:ni].squeeze().copy()
  dZ    = abs(np.diff(ZZ2d, axis=0))

  # 
  # Plot section
  f_plt = False
  if f_plt:
    btx = 'calc_transp_rtofs.py'
    plt.ion()

    stl = ('RTOFS-{0} {1} {2}, {3}/{4}/{5}'.\
            format(expt, sttl, sctnm, YR, MM, DD))
    importlib.reload(mom6vld)
    mom6vld.plot_xsect(XX, Hb, ZZ2d, A2d, HH, stl=stl, ijsct=[i1,i2,j1,j2], btx=btx)


  # Vol transport:
  VFlx = mom6util.vol_transp_2Dsection(LSgm, ZZ2d, A2d)

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
  if irec%10 == 0 or irec == nrec-1:
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

