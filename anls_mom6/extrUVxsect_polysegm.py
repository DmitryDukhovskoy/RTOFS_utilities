"""
  Extract 2D U field along a poly-section 
  comprised of several connected segments
  Slanted sections are allowed
  
  For transports: only 1/2 of grid cell at the ends of the segments

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

# Check orientation of the line/norms 
# curve_ornt
f_cont = False     # load saved and start from last saved
#sctnm = 'Fram79'
#sctnm = 'DavisStr'
#sctnm = 'Yucatan'  # 2-leg right angle section
#sctnm = 'Yucatan2'  # slented section
#sctnm = 'DavisS2'
#sctnm = 'Fram79s2'
#sctnm = 'BarentsS'
#sctnm = 'BeringS'
#sctnm = 'DenmarkS'
#sctnm = 'IclShtl'
#sctnm = 'ShtlScot'
#sctnm = 'LaManch'
sctnm = 'NAtl39'
fld2d = 'Unrm'
dday  = 5       # time stepping for data processing/analysis
dnmb1 = mtime.datenum([2021,1,1])
dnmb2 = mtime.datenum([2021,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)

expt  = '003'
hg    = 1.e15

pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'
pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/' + \
          'MOM6_CICE6/expt{0}/'.format(expt)
#floutp = 'mom6-{4}_u2dsect_{0}{1:02d}-{2}{3:02d}_{5}.pkl'.\
#         format(dv1[0], dv1[1], dv2[0], dv2[1], expt, sctnm)
floutp = f"mom6-{expt}_{fld2d}VFlx_{dv1[0]}" + \
         f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
ffout = pthoutp + floutp

print(f"\nExtracting {fld2d} {sctnm} \n")

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
                      np.arange(-2000., -9001., -1000)))

import mod_misc1 as mmisc
import mod_mom6 as mom6util
importlib.reload(mom6util)

TPLT = mom6util.create_time_array(dnmb1, dnmb2, dday, date_mat=True)
nrec = len(TPLT)

dnmbL = 0
if f_cont:
  print('Start from last saved record')
  print('Loading ' + ffout)
  try:
    with open(ffout,'rb') as fid:
      F2D, UFLX = pickle.load(fid)
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

# Define segment lengths, norms, etc:
# positive: norm vector is to the left as follow the section
# negative: to the right
    DX, DY   = mom6util.dx_dy(LON,LAT)
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

  # Read layer thicknesses:
  dH     = mom6util.read_mom6(flin, 'h', finfo=False)
  ssh    = mom6util.read_mom6(flin, 'SSH', finfo=False)
  ZZ, ZM = mom6util.zz_zm_fromDP(dH, ssh, f_intrp=True, finfo=False)
  #ZZ, ZM = mom6util.get_zz_zm(flin)

# Read U/V only normal velocity needed:
  U3d = mom6util.read_mom6(flin, 'u', finfo=False)
  V3d = mom6util.read_mom6(flin, 'v', finfo=False)

# Extract fields for the section
# Note that all legs are combined in one 2D field along the section
  Uj   = U3d[:,JJ,II]
  Ujm1 = U3d[:,JJ,II-1]
  Uj   = np.where(np.isnan(Uj),0.,Uj)
  Ujm1 = np.where(np.isnan(Ujm1),0.,Ujm1)
  U2d  = 0.5*(Uj+Ujm1)

  Vi   = V3d[:,JJ,II]
  Vim1 = V3d[:,JJ-1,II]
  Vi   = np.where(np.isnan(Vi),0.,Vi)
  Vim1 = np.where(np.isnan(Vim1),0.,Vim1)
  V2d  = 0.5*(Vi+Vim1) 
     
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

  f_prnt=False
  if f_prnt: 
    ij= 6
    i = II[ij]
    j = JJ[ij]
#    mmisc.print_1col(ZZ[:,j,i])
    mmisc.print_2col(ZZ[:,j,i],dH[:,j,i], prc=3)
    mmisc.print_2col(ZZ2d[:,ij],dZ[:,ij], prc=3)
    mmisc.print_2col(U2d[:,ij],V2d[:,ij], prc=3)
    mmisc.print_2col(Unrm1[:,ij],Unrm2[:,ij], prc=3)
    mmisc.print_2col(VFlx1,VFlx2, prc=3)

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
  # Check: Plot section
  btx = 'extrUVxsect_polysegm.py' 
  f_chck = False
  if f_chck:
  # Draw section with half-segments and norm vectors
    ax1 = mom6vld.plot_section_map(II, JJ, IJ, Vnrm1, Vnrm2, II_hf, JJ_hf, \
                                   HH, fgnmb=1, btx=btx)
    F = STOP

  f_plt = False
  if f_plt:
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
      rmin = -0.2
      rmax = 0.2

# Get interf depths for half-segments
#    ZZ_hf = mom6vld.segm_half_zintrf(II, JJ, ZZ2d)

#    nHf = len(Hb_hf)
    XI  = np.arange(0, nsgm, 1, dtype=int)

    IJs = np.array([II,JJ]).transpose()
    sttl = 'UV (x+ y+) '
#    stl = ('0.08 MOM6-CICE6-{0} {1} {2}, {3}/{4}/{5}'.\
#            format(expt, sttl, sctnm, YR, MM, DD))
    stl = f"0.08 MOM6-CICE6-{expt} {sttl} {sctnm}, {YR}/{MM:02d}/{DD:02d}"

    mom6vld.plot_xsect(XI, Hb, ZZ2d, UV2d, HH, stl=stl,\
                       rmin = rmin, rmax = rmax, clrmp=cmpr,\
                       IJs=IJs, btx=btx)

    stl = ('0.08 MOM6-CICE6-{0} Interpolated  {1} {2}, {3}/{4}/{5}'.\
            format(expt, sttl, sctnm, YR, MM, DD))

    mom6vld.plot_xsect(XI, Hb, ZZi, UV2di, HH, stl=stl, fgnmb=2,\
                       rmin = rmin, rmax = rmax, clrmp=cmpr,\
                       IJs=IJs, btx=btx)
    C = G

  # 2D field
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
    print('Saving to '+ffout)
    with open(ffout,'wb') as fid:
      pickle.dump([F2D, UFLX],fid)

print(' All done ' )


f_chck2=False
if f_chck2:
  VFlx = UFLX.trnsp
  VFlx_mn = np.mean(np.nansum(VFlx, axis=1)) # mean Transp m3/s





