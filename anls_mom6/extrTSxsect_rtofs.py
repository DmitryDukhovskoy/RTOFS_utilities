"""
  Extract 2D T/S field along a poly-section 
  comprised of several connected segments
  Slanted sections are allowed

  RTOFS
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

import mod_read_hycom as mhycom
import mod_time as mtime
from mod_utils_fig import bottom_text

import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)

# Check orientation of the line/norms 
# curve_ornt
#sctnm = 'DavisS2'
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
#sctnm = 'NAtl39'
#======= Ocean Sections =====
#sctnm = 'BaffNAFram'
#sctnm = 'AlaskaIcld' 
sctnm = 'GoMCarib'

f_save    = True
f_cont    = False    # load saved and start from last saved
f_chcksgm = False  # plot section with all segments and norms
f_plt     = True # Plot sections on MOM grid and interpolated to Z

fld2d = 'salt'
#fld2d = 'potT'
dday  = 7       # time stepping for data processing/analysis

expt  = 'product'
sfx   = 'n-24'
YR1   = 2023
hg    = 1.e15
huge  = 1.e15
rg    = 9806.

dnmb1 = mtime.datenum([YR1,1,2])
dnmb2 = mtime.datenum([YR1,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)
ds1   = mtime.datestr(dnmb1)
ds2   = mtime.datestr(dnmb2)

# HYCOM notations:
fld = fld2d
if fld2d == 'salt':
  fld = 'salin'
elif fld2d == 'potT':
  fld = 'temp'

print(f"\nExtracting {sctnm} {fld} RTOFS-{expt} {ds1}-{ds2}\n") 

pthrun  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/wcoss2.prod/'
pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/RTOFS_production/'
floutp  = f"rtofs-{expt}_{fld}xsct_{dv1[0]}" + \
          f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
ffout   = pthoutp + floutp

pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
idm = HH.shape[1]
jdm = HH.shape[0]

# Either strait or ocean section:
STR = mom6vld.ocean_straits()
if sctnm in STR:
  oc_strait = True

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

else:
  STR = mom6vld.ocean_sections()
  if sctnm in STR:
    oc_strait = False
    NPsct   = STR[sctnm]["NP"]
    Is      = STR[sctnm]["II"]
    Js      = STR[sctnm]["JJ"]
    Is      = np.array(Is)
    Js      = np.array(Js)
# Replace with N. Pole values if needed:
# MOM6 and RTOFS have differen # of jdm
    if NPsct:
      dJ = abs(Js - jdm)
      Js = np.where(dJ < 3, jdm-1, Js)

    nlegs   = len(Is) - 1
    IJ      = np.zeros((nlegs+1,2))
    IJ[:,0] = Is
    IJ[:,1] = Js
    
  else:
    raise Exception(f"Name {sctnm} is not defined as a strait or section")

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
dday = 7
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
  pthbin  = pthrun + f"rtofs.{YR}{MM:02d}{DD:02d}/"
  flhycom = f"rtofs_glo.t00z.{sfx}.archv"
  fina    = pthbin + flhycom + '.a'
  finb    = pthbin + flhycom + '.b'

  if irec == 0:

# Define segment lengths, norms, etc:
# positive: norm vector is to the left as follow the section
# negative: to the right
    DX, DY = mhycom.dx_dy(LON,LAT)
    SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive')
#    SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='negative')
    II     = SGMT.I_indx
    JJ     = SGMT.J_indx
    hLsgm1 = SGMT.half_Lsgm1
    hLsgm2 = SGMT.half_Lsgm2
    Vnrm1  = SGMT.half_norm1
    Vnrm2  = SGMT.half_norm2
    nLeg   = SGMT.Leg_number
    XX     = LON[JJ,II]
    YY     = LAT[JJ,II]
    Hb     = HH[JJ,II]
    LSgm   = np.zeros((len(II)))  # total segment length = half1 + half2
    for ik in range(len(II)):
       LSgm[ik] = hLsgm1[ik] + hLsgm2[ik]

  if dnmb <= dnmbL:
    print(f"Skipping rec {irec+1} {YR}/{MM}/{DD} --->")
    continue

  # Read RTOFS
  # Read layer thicknesses:
  print(f"Reading thknss {fina}")
  dH, idm, jdm, kdm = mhycom.read_hycom(fina, finb, 'thknss', finfo=False)
  dH = dH/rg
  dH = np.where(dH>huge, np.nan, dH)
  dH = np.where(dH<0.001, 0., dH)
  ZZ, ZM = mhycom.zz_zm_fromDP(dH, f_btm=False, finfo=False)

  print(f"Reading {fld} {fina}")
  A3d, _, _, _ = mhycom.read_hycom(fina, finb, fld, finfo=False)
  A3d = np.where(A3d >= huge, np.nan, A3d)

# Subset to section:
  dH2d  = dH[:,JJ,II].squeeze()
  Hb    = HH[JJ,II].squeeze()
  ZZ2d  = ZZ[:,JJ,II].squeeze()
  ZM2d  = ZM[:,JJ,II].squeeze()
  dZ    = abs(np.diff(ZZ2d, axis=0))

  A2d  = A3d[:,JJ,II]
  A2d  = mom6util.fill_bottom(A2d, ZZ2d, Hb, fill_land=True)
  A2di = mom6vld.interp_2Dsect(A2d, ZZi, ZZ2d, ZM2d, Hb) # interpolate on z leveles

# Check segments etc:
# Draw section line with all segments, norm vectors and
  if f_chcksgm:
    IIhf = []  # half-segments
    JJhf = []
    nI = len(Hb)
    XI = np.arange(0, nI, 1, dtype=int)

    ax2 = mom6vld.plot_section_map(II, JJ, IJ, Vnrm1, Vnrm2, IIhf, JJhf, \
                       HH, fgnmb=3, XI=XI, dX=100, btx='extrTSxsect_polysegm.py')
    A = STOP

  # Plot section
  btx = 'extrTSxsect_polysegm.py'
  if f_plt:
    plt.ion()

    importlib.reload(mom6vld)
    import mod_utils as mutil
    import plot_sect as psct
    import mod_colormaps as mcmp
    importlib.reload(mcmp)

#    cmpr = mutil.colormap_ssh(nclrs=100)
    if fld2d == 'salt':
#      cmpr = mutil.colormap_salin(clr_ramp=[0.94,0.95,1])
      cmpr = mcmp.colormap_haline()
      rmin = STR[sctnm]["smin"]
      rmax = STR[sctnm]["smax"]
    elif fld2d == 'potT':
#      cmpr = mutil.colormap_temp2()
      cmpr = mcmp.colormap_temp2()
      rmin = STR[sctnm]["tmin"]
      rmax = STR[sctnm]["tmax"]
    elif fld2d == 'Unrm':
      cmpr = mutil.colormap_ssh(nclrs=100)
      rmin = STR[sctnm]["umin"]
      rmax = STR[sctnm]["umax"]


    nI = len(Hb)
    XI = np.arange(0, nI, 1, dtype=int)
#    XI = XX.copy()

    IJs = np.array([II,JJ]).transpose()
    sttl = fld
    stl = f"0.08 RTOFS-{expt} {sttl} {sctnm}, {YR}/{MM:02d}/{DD:02d}"

    mom6vld.plot_xsect(XI, Hb, ZZ2d, A2d, HH, stl=stl,\
                       rmin = rmin, rmax = rmax, clrmp=cmpr,\
                       IJs=IJs, btx=btx)

    stl = ('0.08 RTOFS-{0} Interpolated  {1} {2}, {3}/{4}/{5}'.\
            format(expt, sttl, sctnm, YR, MM, DD))

    mom6vld.plot_xsect(XI, Hb, ZZi, A2di, HH, stl=stl, fgnmb=2,\
                       rmin = rmin, rmax = rmax, clrmp=cmpr,\
                       IJs=IJs, btx=btx)
    F = STOP

  # 2D field
  if irec == 0:
    F2D  = mom6vld.UTS2D(dnmb, II, JJ, XX, YY, \
                         LSgm, ZZi, Hb, A2di)
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
    if f_save: 
      print('Saving to '+ffout)
      with open(ffout,'wb') as fid:
        pickle.dump(F2D,fid)
    else:
      print('Saving turned off \n')

print(' All done ' )

f_chck=False
if f_chck:
#  II, JJ = mmisc.xsect_indx(IIv,JJv)
  xl1 = np.min(IJ[:,0])
  xl2 = np.max(IJ[:,0])
  yl1 = np.min(IJ[:,1])
  yl2 = np.max(IJ[:,1])
  dxy = 50

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
  ax1.contour(HH,[0], colors=[(0,0.,0)])
  ax1.contour(HH,[-1000,-500], colors=[(0.8,0.8,0.8)], linestyles='solid')
  ax1.axis('scaled')
  ax1.set_xlim([xl1-dxy, xl2+dxy])
  ax1.set_ylim([yl1-dxy, yl2+dxy])
  ax1.plot(IJ[:,0],IJ[:,1],'b-')
  ax1.plot(IJ[:,0],IJ[:,1],'.', ms=15, color=(0.,0.5,0.8))

  nsgm = len(II)
  for isgm in range(nsgm):
    ii0 = II[isgm]
    jj0 = JJ[isgm]
    if isgm < nsgm-1:
      ip1 = II[isgm+1]
      jp1 = JJ[isgm+1]
    else:
      ip1 = ii0
      jp1 = jj0
    if isgm > 0:
      im1 = II[isgm-1]
      jm1 = JJ[isgm-1]
    else:
      im1 = ii0
      jm1 = jj0
    ih1 = 0.5*(im1+ii0)
    jh1 = 0.5*(jm1+jj0)
    ih2 = 0.5*(ip1+ii0)
    jh2 = 0.5*(jp1+jj0)

    ax1.plot(ii0, jj0, '.', ms=10, color=(0., 0.2, 0.5))
    if isgm > 0:
      ax1.plot([ih1,ii0],[jh1,jj0],'r-') # 1st half of the segment
    if isgm < nsgm-1:
      ax1.plot([ii0,ih2],[jj0,jh2],'g-') # 2nd half
# Plot norm
    scl = 0.2
    v1  = scl*Vnrm1[isgm,:]
    v2  = scl*Vnrm2[isgm,:]
    in1 = 0.5*(ii0-ih1)+ih1
    jn1 = 0.5*(jj0-jh1)+jh1
    in2 = 0.5*(ih2-ii0)+ii0
    jn2 = 0.5*(jh2-jj0)+jj0
    if isgm > 0:
      ax1.plot([in1, in1+v1[0]],[jn1, jn1+v1[1]],'-', color=(0.7,0.2,0))
    if isgm < nsgm-1:
      ax1.plot([in2, in2+v2[0]],[jn2, jn2+v2[1]],'-', color=(0.,1.,0.4))

    bottom_text(btx,pos=[0.08, 0.08])    

f_chck2=False
if f_chck2:
  VFlx = UFLX.trnsp
  VFlx_mn = np.mean(np.nansum(VFlx, axis=1)) # mean Transp m3/s





