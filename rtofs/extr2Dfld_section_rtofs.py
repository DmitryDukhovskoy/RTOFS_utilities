"""
  RTOFS - production
  scripts to bring rtofs files: scripts/rtofs/get_production_hycom_Ndays.sh

  Extract 2D field along a section
  Note variables are on staggered hycom grid in the output archv files
  U,V (i,j) is different from MOM6

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

f_cont = False       # True - load saved and start from last saved
expt   = 'product'
sfx    = 'n-24'
YR1    = 2021
#sctnm  = 'Fram79'
sctnm = 'DavisStr'
fld2d  = 'Unrm'
#fld2d = 'salin'
#fld2d = 'temp'
huge   = 1.e15
rg     = 9806.

dnmb1 = mtime.datenum([YR1,1,1])
dnmb2 = mtime.datenum([YR1,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)

pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/wcoss2.prod/'
pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/RTOFS_production/'
floutp = f"rtofs-{expt}_{fld2d}xsct_{dv1[0]}" + \
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


# Interpolate onto fixed z-levels:
ZZi = np.concatenate((np.arange(0.,     -21,    -1),
                      np.arange(-22.,   -202.,  -2),
                      np.arange(-205.,  -250.,  -5),
                      np.arange(-250.,  -510.,  -10),
                      np.arange(-525.,  -1000., -25),
                      np.arange(-1000., -2000., -50),
                      np.arange(-2000., -5001., -1000)))

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
      F2D = pickle.load(fid)

    TM    = F2D.TM
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

#j=1250
#i=1000
#mmisc.print_1col(ZZ[:,j,i])
#mmisc.print_2col(ZZ[:,j,i],dH[:,j,i])

  # Read U/V only normal velocity needed:
  if fld2d == 'Unrm':
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
  else:
     fld = fld2d

  if xsct_EW:
    XX = LON[j1,i1:ni].squeeze()
    YY = LAT[j1,i1:ni].squeeze()
    if XX[-1] < XX[0] and XX[-1] < 0.:
      XX = np.where(XX<0, XX+360., XX)
    LSgm = DX[j1,i1:ni].squeeze()
#    dHsgm = dH[:,j1,i1:ni].squeeze()
  else:
    XX = LAT[j1:nj,i1].squeeze()
    YY = LAT[j1:nj,i1].squeeze()
    LSgm = DY[j1:nj,i1].squeeze()

  sttl = fld2d

  A3d, _, _, _ = mhycom.read_hycom(fina, finb, fld)
  A3d = np.where(A3d >= huge, np.nan, A3d)
  # Interpolated onto h-pnt:
  if xsct_EW:
    if fld2d == 'Unrm':
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
      A2d = A3d[:,j1,i1:ni].squeeze()
  else:
    if fld2d == 'Unrm':
# Add barotropic U:
      for kk in range(kdm):
        A3d[kk,:,:] = A3d[kk,:,:] + Ubtrop
# Collocate U components to h pnt
      Ui   = A3d[:,j1:nj,i1].squeeze()
      Uip1 = A3d[:,j1:nj,i1+1].squeeze()
      Ui   = np.where(np.isnan(Ui), 0., Ui)
      Uip1 = np.where(np.isnan(Uip1), 0., Uip1)
      A2d = 0.5*(Ui + Uip1)   # collocated with h
    else:
      A2d = A3d[:,j1:nj,i1].squeeze()

  dH2d  = dH[:,j1:nj,i1:ni].squeeze()
  Hb    = HH[j1:nj,i1:ni].squeeze()
  ZZ2d  = ZZ[:,j1:nj,i1:ni].squeeze().copy()
  ZM2d  = ZM[:,j1:nj,i1:ni].squeeze().copy()
  dZ    = abs(np.diff(ZZ2d, axis=0))

  
#
# Interpolate:
  A2di = mom6vld.interp_2Dsect(A2d, ZZi, ZM2d, Hb)

  # 
  # Plot section
  f_plt = False
  if f_plt:
    btx = 'extr2Dfld_section_rtofs.py'
    plt.ion()

    importlib.reload(mom6vld)
    import mod_utils as mutil
    import plot_sect as psct

#    cmpr = mutil.colormap_ssh(nclrs=100)
    if fld2d == 'salt':
      cmpr = mutil.colormap_salin(clr_ramp=[0.94,0.95,1])
      rmin = 32.5
      rmax = 35.0
    elif fld2d == 'potT':
      cmpr = mutil.colormap_temp()
      rmin = -1.5
      rmax = 6.5
    elif fld2d == 'Unrm':
      cmpr = mutil.colormap_ssh(nclrs=100)
      rmin = -0.1
      rmax = 0.1

#    stl = ('0.08 MOM6-CICE6-{0} {1} {2}, {3}/{4}/{5}'.\
#            format(expt, sttl, sctnm, YR, MM, DD))
    stl = f"RTOFS-{expt} {sttl} {sctnm}, {YR}/{MM:02d}/{DD:02d}"

    mom6vld.plot_xsect(XX, Hb, ZZ2d, A2d, HH, stl=stl,\
                       rmin = rmin, rmax = rmax, clrmp=cmpr,\
                       ijsct=[i1,i2,j1,j2], btx=btx)

    stl = f"RTOFS-{expt} {sttl} {sctnm} Interpelodated, {YR}/{MM:02d}/{DD:02d}"
    II = np.arange(i1,ni)
    mom6vld.plot_xsect(XX, Hb, ZZi, A2di, HH, stl=stl, fgnmb=2,\
                       rmin = rmin, rmax = rmax, clrmp=cmpr,\
                       ijsct=[i1,i2,j1,j2], btx=btx)

#    AB = CC

  if irec == 0:
    F2D = mom6vld.FIELD2D(dnmb, XX, YY, LSgm, ZZi, Hb, A2di)
  else:
    F2D.add_array(dnmb, A2di)

  rtime0 = time.time()
  dtime  = (rtime0-rtimeS)/60.
  print('{0:4d}: {1}/{2:02d}/{3:02d} '.format(irec, YR, MM, DD) + \
        'min/max {4} = {0:5.2f}/{1:5.2f}, {4}_intrp= {2:5.2f}/{3:5.2f}'.\
        format(np.nanmin(A2d), np.nanmax(A2d),\
               np.nanmin(A2di), np.nanmax(A2di), fld2d))
  print('  1 day processed, {0:8.3f} min ...'.format(dtime))

# Save:
  if irec%10 == 0 or irec == nrec-1:
    print('Saving to '+ffout)
    with open(ffout,'wb') as fid:
      pickle.dump(F2D,fid)

print(' All done ' )

f_chck=False
if f_chck:
  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
  ax1.plot(XX,VFlx*1e-6)


