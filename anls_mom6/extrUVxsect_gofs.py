"""
  RTOFS - production
  scripts to bring rtofs files: scripts/rtofs/get_production_hycom_Ndays.sh

  Extract 2D U field along a section
  Note variables are on staggered hycom grid in the output archv files
  U,V (i,j) is different from MOM6

  This is updated code to allow slanted / zigzagging (with steps) sections
  Extract 2D U field along a poly-section 
  comprised of several connected segments

  For transports: only 1/2 of grid cell at the ends of the segments

  GOFS3.1-93.0 restart files
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

f_plt  = True        # Plot original and interpolated sections for checking
f_cont = False       # True - load saved and start from last saved
f_save = False

expt   = '93.0'
YR1    = 2021
#sctnm = 'Yucatan2'  # slanted section
#sctnm = 'DavisStr2'
#sctnm = 'Fram79'
#sctnm = 'Fram80'
sctnm = 'FlorCabl'
#sctnm = 'DavisStr'
#sctnm = 'DavisS2'
#sctnm = 'Fram79s2'
#sctnm = 'BarentsS'
#sctnm = 'BeringS'
#sctnm = 'DenmarkS'
#sctnm = 'IclShtl'
#sctnm = 'ShtlScot'
#sctnm = 'LaManch'
#sctnm = 'NAtl39'

fld2d = 'Unrm'
dday  = 7       # time stepping for data processing/analysis

huge   = 1.e15
rg     = 9806.

dnmb1 = mtime.datenum([YR1,1,1])
dnmb2 = mtime.datenum([YR1,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)

pthrun  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/GOFS3.1/restart/'
pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/GOFS3.1/'
floutp  = f"gofs31-930_{fld2d}xsct_{dv1[0]}" + \
          f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
ffout = pthoutp + floutp

print(f"\nExtracting GOFS3.1 {fld2d} {sctnm} \n")

pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
IDM = HH.shape[1]
JDM = HH.shape[0]


STR = mom6vld.ocean_straits()
nlegs = STR[sctnm]["nlegs"]
I1    = STR[sctnm]["xl1"]
I2    = STR[sctnm]["xl2"]
J1    = STR[sctnm]["yl1"]
J2    = STR[sctnm]["yl2"]
Ni    = np.zeros((nlegs))
Nj    = np.zeros((nlegs))
IJ    = np.zeros((nlegs+1,2))
for kk in range(nlegs):
  Ni[kk]     = I2[kk]+1
  Nj[kk]     = J2[kk]+1
  IJ[kk,0]   = I1[kk]
  IJ[kk,1]   = J1[kk]
  IJ[kk+1,0] = I2[kk]
  IJ[kk+1,1] = J2[kk]

# Interpolate onto fixed z-levels:
ZZi = np.concatenate((np.arange(0.,     -21,    -1),
                      np.arange(-22.,   -202.,  -2),
                      np.arange(-205.,  -250.,  -5),
                      np.arange(-250.,  -510.,  -10),
                      np.arange(-525.,  -1000., -25),
                      np.arange(-1000., -2000., -50),
                      np.arange(-2000., -3000., -100),
                      np.arange(-3000., -4000., -250),
                      np.arange(-4000., -5000., -500),
                      np.arange(-5000., -9001., -1000)))

import mod_misc1 as mmisc
import mod_mom6 as mom6util
importlib.reload(mom6util)

# Time array
irc  = -1
nrec = 48
TPLT = np.zeros(nrec, dtype=[('dnmb', float),
                              ('date', int, (4,)),
                              ('yrday', float)])
for imo in range(dv1[1], dv2[1]+1):
  for mday in range(dv1[2], 29, dday):
    irc += 1
    dnmb = mtime.datenum([YR1,imo,mday])
    DV   = mtime.datevec(dnmb)
    _, jday = mtime.dnmb2jday(dnmb)
    TPLT['dnmb'][irc]   = dnmb
    TPLT['yrday'][irc]  = jday
    TPLT['date'][irc,:] = DV[0:4]


nrec = len(TPLT)

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
  pthbin  = pthrun
  flhycom = f"restart_r{YR}{MM:02d}{DD:02d}00_930"
  fina    = pthbin + flhycom + '.a'
  finb    = pthbin + flhycom + '.b'

  if irec == 0:
  # Calc segment lengths:
    DX, DY = mhycom.dx_dy(LON,LAT)
# Define segment lengths, norms, etc:
# positive: norm vector is to the left as follow the section
# negative: to the right
    SGMT     = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive')
#    SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='negative')
    II       = SGMT.I_indx
    JJ       = SGMT.J_indx
    hLsgm1   = SGMT.half_Lsgm1
    hLsgm2   = SGMT.half_Lsgm2
    Vnrm1    = SGMT.half_norm1
    Vnrm2    = SGMT.half_norm2
    nLeg     = SGMT.Leg_number
    LegNorm  = SGMT.LegNorm
    XX       = LON[JJ,II]
    YY       = LAT[JJ,II]
    Hb       = HH[JJ,II]
    LSgm     = np.zeros((len(II)))  # total segment length = half1 + half2
    for ik in range(len(II)):
       LSgm[ik] = hLsgm1[ik] + hLsgm2[ik]

    II_hf, JJ_hf, XX_hf, \
    YY_hf, Hb_hf, LSgm_hf = mom6vld.segm_half_coord(II, JJ, \
                            hLsgm1, hLsgm2, XX, YY, Hb)

# Define weights for U and V for each segment:
    nsgm  = len(II)
    whtU1 = np.zeros((nsgm))
    whtV1 = np.zeros((nsgm))
    whtU2 = np.zeros((nsgm))
    whtV2 = np.zeros((nsgm))
    for isgm in range(nsgm):
      nrm1   = Vnrm1[isgm]  # normal for 1st half of the segment
      nrm2   = Vnrm2[isgm]  # --"-- for 2nd half
      ii0    = II[isgm]
      jj0    = JJ[isgm]
      hlsgm1 = hLsgm1[isgm]
      hlsgm2 = hLsgm2[isgm]
      lnrm1  = np.sqrt(nrm1[0]**2 + nrm1[1]**2)
      lnrm2  = np.sqrt(nrm2[0]**2 + nrm2[1]**2)
# No zero-length norms:
      if lnrm1 < 1.e-30 and lnrm2 < 1.e-30:
        raise Exception(f"segm {isgm} has 0 normal vectors")

      whtU1[isgm] = nrm1[0]
      whtU2[isgm] = nrm2[0]
      whtV1[isgm] = nrm1[1]
      whtV2[isgm] = nrm2[1]

  if dnmb <= dnmbL:
    print(f"Skipping rec {irec+1} {YR}/{MM}/{DD} --->")
    continue

  if not os.path.isfile(fina):
    print(f"!!! missing {fina}")
    print(' ===== skipping =====>')
    continue

  # Read layer thicknesses:
  print(f"Reading thknss {fina}")
  dH, kdm = mhycom.read_hycom_restart(fina, finb, 'dp', IDM, JDM, finfo=True)
  dH = dH/rg
  dH = np.where(dH>huge, np.nan, dH)
  dH = np.where(dH<0.001, 0., dH)
  ZZ, ZM = mhycom.zz_zm_fromDP(dH, f_btm=False, finfo=False)

# Read barotropic U
# Only for archv output !
  F,  _ = mhycom.read_hycom_restart(fina,finb,'ubavg', IDM, JDM)
  F = np.where(F>=huge, 0., F)
  Ubtrop = F.copy()
# Read barotropic V
  F, _ = mhycom.read_hycom_restart(fina,finb,'vbavg', IDM, JDM)
  F = np.where(F>=huge, 0., F)
  Vbtrop = F.copy()

# Read U/V:
  U3d, _ = mhycom.read_hycom_restart(fina, finb, 'u', IDM, JDM, finfo=False)
  U3d = np.where(U3d >= huge, np.nan, U3d)
  V3d, _ = mhycom.read_hycom_restart(fina, finb, 'v', IDM, JDM, finfo=False)
  V3d = np.where(V3d >= huge, np.nan, V3d)
# Add barotropic U:
  for kk in range(kdm):
    U3d[kk,:,:] = U3d[kk,:,:] + Ubtrop
    V3d[kk,:,:] = V3d[kk,:,:] + Vbtrop

# Extract fields for the section
# Note that all legs are combined in one 2D field along the section
# Collocate U and V with h-points
  Ui   = U3d[:,JJ,II]
  Uip1 = U3d[:,JJ,II+1]
  Ui   = np.where(np.isnan(Ui),0.,Ui)
  Uip1 = np.where(np.isnan(Uip1),0.,Uip1)
  U2d  = 0.5*(Ui+Uip1)

  Vj   = V3d[:,JJ,II]
  Vjm1 = V3d[:,JJ-1,II]
  Vj   = np.where(np.isnan(Vj),0.,Vj)
  Vjm1 = np.where(np.isnan(Vjm1),0.,Vjm1)
  V2d  = 0.5*(Vj+Vjm1)

  dH2d = dH[:,JJ,II].squeeze()
  ZZ2d = ZZ[:,JJ,II].squeeze()
  ZM2d = ZM[:,JJ,II].squeeze()
  dZ   = abs(np.diff(ZZ2d, axis=0))

# Arrays for half segments
# For plotting 2D section keep original sign of U and V
  UV1 = U2d.copy()*0.
  UV2 = U2d.copy()*0.
  for isgm in range(nsgm):
    UV1[:,isgm] = U2d[:,isgm]*abs(whtU1[isgm]) + \
                  V2d[:,isgm]*abs(whtV1[isgm])
    UV2[:,isgm] = U2d[:,isgm]*abs(whtU2[isgm]) + \
                  V2d[:,isgm]*abs(whtV2[isgm])

# For transport - use normal vector direction
  Unrm1 = U2d.copy()*0.
  Unrm2 = U2d.copy()*0.
  for isgm in range(nsgm):
    Unrm1[:,isgm] = U2d[:,isgm]*whtU1[isgm] + \
                    V2d[:,isgm]*whtV1[isgm]
    Unrm2[:,isgm] = U2d[:,isgm]*whtU2[isgm] + \
                    V2d[:,isgm]*whtV2[isgm]

#
# Vol transport:
  VFlx1 = mom6util.vol_transp_2Dsection(hLsgm1, ZZ2d, Unrm1)
  VFlx2 = mom6util.vol_transp_2Dsection(hLsgm2, ZZ2d, Unrm2)
# Net transport: combine half grid cells:
  VFlx = VFlx1 + VFlx2

#  mmisc.print_2col(VFlx1,VFlx2, prc=3)

# Data for plotting 2D sections
# Project U on the normal vector for the main section line
  UV2d = np.zeros((kdm,nsgm))
  for isgm in range(nsgm):
# indices for 1st half segment 
    uu = U2d[:,isgm]
    vv = V2d[:,isgm]
    Snrm = LegNorm[isgm,:]
    UV2d[:,isgm] = uu*Snrm[0] + vv*Snrm[1]

# Interpolate on Z levels:
#  UV2di = mom6vld.interp_2Dsect_segmhalf(UV2d, ZZi, ZZ2d, ZM2d, Hb)
  UV2di = mom6vld.interp_2Dsect(UV2d, ZZi, ZZ2d, ZM2d, Hb)


  # 
  # Plot section
  btx = 'extrUVxsect_gofs.py'
  if f_plt:
    plt.ion()

    importlib.reload(mom6vld)
    import mod_utils as mutil
    import plot_sect as psct

    cmpr = mutil.colormap_ssh(nclrs=100)
    rmin = STR[sctnm]["umin"]
    rmax = STR[sctnm]["umax"]

    VFtot = np.nansum(VFlx)*1e-6
    XI  = np.arange(0, nsgm, 1, dtype=int)

    IJs = np.stack((II,JJ)).T

    stl = f"GOFS3.1-{expt} Unrm {sctnm}, {YR}/{MM:02d}/{DD:02d}, VF={VFtot:4.1f} Sv"

    mom6vld.plot_xsect(XI, Hb, ZZ2d, UV2d, HH, stl=stl,\
                       rmin = rmin, rmax = rmax, clrmp=cmpr,\
                       zstart=-200., lstart=10, IJs=IJs, btx=btx)

    stl = f"GOFS3.1-{expt} Unrm {sctnm} Interpelodated, {YR}/{MM:02d}/{DD:02d}"

    mom6vld.plot_xsect(XI, Hb, ZZi, UV2di, HH, stl=stl, fgnmb=2,\
                       rmin = rmin, rmax = rmax, clrmp=cmpr,\
                       IJs=IJs, btx=btx)


    AB = STOP

  if irec == 0:
    F2D  = mom6vld.UTS2D(dnmb, II, JJ, XX, YY, \
                         LSgm, ZZi, Hb, UV2di)
    UFLX = mom6vld.TRANSP(dnmb, VFlx, XX, YY, Hb)
  else:
    F2D.add_array(dnmb, UV2di)
    UFLX.add_array(dnmb, VFlx)


  rtime0 = time.time()
  dtime  = (rtime0-rtimeS)/60.
  print('{0:4d}: {1}/{2:02d}/{3:02d} '.format(irec, YR, MM, DD) + \
        'min/max {4} = {0:5.2f}/{1:5.2f}, {4}_intrp= {2:5.2f}/{3:5.2f}'.\
        format(np.nanmin(UV2d), np.nanmax(UV2d),\
               np.nanmin(UV2di), np.nanmax(UV2di), fld2d))
  print('  1 day processed, {0:8.3f} min ...'.format(dtime))

# Save:
  if irec%10 == 0 or irec == nrec-1:
    if f_save:
      print('Saving to '+ffout)
      with open(ffout,'wb') as fid:
        pickle.dump([F2D, UFLX],fid)
    else:
      print('Not saved ...')

print(' All done ' )

f_chck=False
if f_chck:
  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
  ax1.plot(XX,VFlx*1e-6)


