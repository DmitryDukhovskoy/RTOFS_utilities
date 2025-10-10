"""
  Extract 2D field along a section 1D - straight line
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

PPTHN = '/home/Dmitry.Dukhovskoy/python'
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

import mod_time as mtime
from mod_utils_fig import bottom_text

import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)

expt  = '003'
hg    = 1.e15
#sctnm = 'Fram79'
sctnm = 'DavisStr2'
#sctnm = 'DavisStr'
fld2d = 'Unrm'
#fld2d = 'salt'
#fld2d = 'potT'

f_plt = True

dnmb1 = mtime.datenum([2021,1,1])
dnmb2 = mtime.datenum([2021,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)

pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'
pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/' + \
          'MOM6_CICE6/expt{0}/'.format(expt)
#floutp = 'mom6-{4}_u2dsect_{0}{1:02d}-{2}{3:02d}_{5}.pkl'.\
#         format(dv1[0], dv1[1], dv2[0], dv2[1], expt, sctnm)
floutp = f"mom6-{expt}_{fld2d}xsct_{dv1[0]}" + \
         f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
ffout = pthoutp + floutp

STR = mom6vld.ocean_straits()
i1  = STR[sctnm]["xl1"]
i2  = STR[sctnm]["xl2"]
j1  = STR[sctnm]["yl1"]
j2  = STR[sctnm]["yl2"]
ni  = i2+1
nj  = j2+1

if j1 == j2:
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

dday=5
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
  ZZ, ZM = mom6util.zz_zm_fromDP(dH, ssh, f_intrp=True, finfo=False)
  #ZZ, ZM = mom6util.get_zz_zm(flin)

  # Calc segment lengths:
  if irec==0:
    DX, DY = mom6util.dx_dy(LON,LAT)

  #j=1250
  #i=1000
  #mmisc.print_1col(ZZ[:,j,i])
  #mmisc.print_2col(ZZ[:,j,i],dH[:,j,i])

  # Read U/V only normal velocity needed:
  if fld2d == 'Unrm':
    if xsct_EW:
      fld = 'v'
    else:
      fld = 'u'
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

  A3d = mom6util.read_mom6(flin, fld, finfo=False)

  # Interpolated onto h-pnt:
  if xsct_EW:
    if fld2d == 'Unrm':
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
  # Collocate U components to h pnt
      Ui   = A3d[:,j1:nj,i1].squeeze()
      Uim1 = A3d[:,j1:nj,i1-1].squeeze()
      Ui   = np.where(np.isnan(Ui), 0., Ui)
      Uim1 = np.where(np.isnan(Uim1), 0., Uim1)
      A2d = 0.5*(Ui + Uim1)   # collocated with h
    else:
      A2d = A3d[:,j1:nj,i1].squeeze()

  dH2d  = dH[:,j1:nj,i1:ni].squeeze()
  Hb    = HH[j1:nj,i1:ni].squeeze()
  ZZ2d  = ZZ[:,j1:nj,i1:ni].squeeze()
  ZM2d  = ZM[:,j1:nj,i1:ni].squeeze()
  dZ    = abs(np.diff(ZZ2d, axis=0))
#
# Interpolate:
  A2di = mom6vld.interp_2Dsect(A2d, ZZi, ZM2d, Hb)

  # 
  # Plot section
  if f_plt:
    btx = 'extr2Dfld_section.py'
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
    stl = f"0.08 MOM6-CICE6-{expt} {sttl} {sctnm}, {YR}/{MM:02d}/{DD:02d}"

    mom6vld.plot_xsect(XX, Hb, ZZ2d, A2d, HH, stl=stl,\
                       rmin = rmin, rmax = rmax, clrmp=cmpr,\
                       ijsct=[i1,i2,j1,j2], btx=btx)

    stl = ('0.08 MOM6-CICE6-{0} Interpolated  {1} {2}, {3}/{4}/{5}'.\
            format(expt, sttl, sctnm, YR, MM, DD))

    II = np.arange(i1,ni)
    mom6vld.plot_xsect(XX, Hb, ZZi, A2di, HH, stl=stl, fgnmb=2,\
                       rmin = rmin, rmax = rmax, clrmp=cmpr,\
                       ijsct=[i1,i2,j1,j2], btx=btx)

    A = STOP
#    i0 = 3513-i1
#    mmisc.print_1col(A2di[:,i0], prc=4)
#    mmisc.print_1col(A2d[:,i0], prc=4)
#
#    fig2 = plt.figure(2,figsize=(9,9))
#    plt.clf()
#    ax1 = plt.axes([0.1, 0.1, 0.35, 0.8])
#    psct.plot_profile(ax1, Aint, ZZi, 0., zlim=-1700., xl1=-0.1, xl2=0.6) 
#    psct.plot_profile(ax1, A2d[:,i0], ZM2d[:,i0], 0., zlim=-1700., \
#                      xl1=-0.1, xl2=0.6, \
#                      lclr=(0.9, 0.2, 0)) 
#

  # 2D field
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
  Is = i1
  Ie = i2
  Js = j1
  Je = j1+150
  IIv=[Is,Ie]
  JJv=[Js,Je]
  II, JJ = mmisc.xsect_indx(IIv,JJv)

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
  ax1.contour(HH,[0], colors=[(0,0.,0)])
  ax1.contour(HH,[-1000,-500], colors=[(0.8,0.8,0.8)], linestyles='solid')
  ax1.set_xlim([Is-50, Ie+50])
  ax1.set_ylim([Js-50, Je+50])
  ax1.plot(IIv,JJv,'r')
  ax1.plot(II, JJ, '.-')


#  ax1.plot(XX,VFlx*1e-6)

