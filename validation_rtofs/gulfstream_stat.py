"""
  Computed MHD scores for Gulf Stream
  Find Gulfstream north wall using definition of
  Halkin & Rossby, 1985 The Structure and Transport of the Gulf Stream at 73°W
  intersection of the 12C and 400m isobath

  NAVO Guld Stream north wall data:
  From Gwen: 
   On Hera, you can find the NAVO Gulf Stream north wall data files at /scratch2/NCEPDEV/ovp/Lichuan.Chen/DCOM/$YYYYMMDD/wtxtbul.  We use the gs*.sub to generate the plots at https://polar.ncep.noaa.gov/global/fronts/.  

   You can also find the data files on WCOSS2.  They are located at /lfs/h1/ops/prod/dcom/$YYYYMMDD/wtxtbul (e.g., /lfs/h1/ops/prod/dcom/20230407/wtxtbul/gs091nw.sub). 

# hour and output to plot wrt to rdate - f/cast date:
# -24, ..., 0, ..., 12
# for hr <=0, output from nowcast/hindcast "n-24" - initial after incup,
#    n-12, ..., n00 - initial for f/cast
# for hr >0: f01, ..., f12, ... f196 (if run that long) - forecasts


  Dmitry Dukhovskoy, NOAA/NCEP/EMC July 2023

"""
import os
import numpy as np
import sys
import importlib
import matplotlib
import matplotlib.pyplot as plt
import datetime
import pickle

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hausdorff')

import mod_read_hycom as mhycom
import mod_misc1 as mmisc
import mod_hausdorff_distance as mmhd
import mod_utils_fig as mufig
import mod_utils as mutil
import mod_rtofs as mrtofs
importlib.reload(mmisc)
import mod_gulfstream as mgulf
importlib.reload(mgulf)
import mod_time as mtime
importlib.reload(mtime)
from mod_utils_fig import bottom_text

# Interpolation depth
z0 = -400.

#expt   = 'paraD5'  # paraXX - parallel runs or product - production
expt    = 'product'
rdateS  = '20230303'
rdateE  = '20230531'

#rdate0  = '20230419'
sfx     = 'n-24'

# Figure output directory:
f_save = 2   # =0 - not save, =1 rewrite/save, =2 - cotinue from saved
f_intract = True  # False - no figures shown
#pthfig  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/' + expt + '/fig/'
pthscr  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/'
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
pthsave = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/data_outp/validation_rtofs/'

if not f_intract:
  print('Interactive mode is OFF')
  matplotlib.use('Agg')
  plt.close('all')
  plt.ioff()
else:
  plt.ion()

# Function to print mouse click event coordinates
def onclick(event):
   print([event.xdata, event.ydata])


# Define region of interest:
# For the HYCOM  0.08 Global grid !!! 
II = [2565, 2565, 2950, 2950]
JJ = [1809, 2190, 2190, 1809]

# Region to exclude shallow regions
ISH = [2562, 2569, 2574, 2579, 2581, 2609, 2721, 3012, 3012, 2562]
JSH = [1808, 1812, 1816, 1823, 1866, 1911, 1991, 2012, 1772, 1770]



dnmbS = int(mtime.rdate2datenum(rdateS))
dnmbE = int(mtime.rdate2datenum(rdateE))

flnmout  = 'MHD' + expt + '.pkl' 
fmhd_out = pthsave + flnmout 

icnt  = 0
MHDGS = []
TM    = []
dlast = -1

if f_save > 1:
  print('Loading saved ' + fmhd_out)
  with open(fmhd_out,'rb') as fid:
    [MHDGS,TM] = pickle.load(fid)
    dlast = int(TM[-1])
    MHDGS = list(MHDGS)
    TM    = list(TM)
    rdate0 = mtime.dnumb2rdate(dlast, ihours=False)
    print('Last saved record = ' + rdate0)

for dnmb0 in range(dnmbS,dnmbE+1):
  rdate0 = mtime.dnumb2rdate(dnmb0, ihours=False)
  print('Processing ' + rdate0)

  if dnmb0 <= dlast:
    icnt += 1
    print(rdate0 + ' processed, skipping ...')
    continue

  if expt == 'product':
    pthbin = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/wcoss2.prod' +\
             '/rtofs.' + rdate0 + '/'
  else:
    pthbin = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/' + expt +\
         '/rtofs.' + rdate0 + '/'

  flarchv= 'rtofs_glo.t00z.' + sfx + '.archv'
  fina   = pthbin + flarchv + '.a'
  finb   = pthbin + flarchv + '.b'

  if not os.path.isfile(fina):
    print('File not found ' + fina)
    print('skipping ...')
    continue

  if len(rdate0) > 0:
    YR     = int(rdate0[0:4])
    MM     = int(rdate0[4:6])
    DD     = int(rdate0[6:8])
    yrday  = mtime.rdate2jday(rdate0)
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


  huge = 1.e20
  rg   = 9806.
  IDM, JDM, KDM = mhycom.hycom_dim(fina,finb)

  ftopo = 'regional.depth'
  fgrid = 'regional.grid'
  LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

  if not 'X' in vars() or not 'MS' in vars() or not 'Msh' in vars():
    X, Y = np.meshgrid(np.arange(IDM), np.arange(JDM))
    MS, IM, JM    = mmisc.inpolygon_v2(X,Y,II,JJ)  # Gulfstream region
    Msh, Ish, Jsh = mmisc.inpolygon_v2(X,Y,ISH,JSH)  # shelves to exclude

  # Read layer pressures:
  dH = np.zeros((KDM,JDM,IDM))
  for kk in range(1,KDM+1):
    F,nn,mm,lmm = mhycom.read_hycom(fina,finb,'thknss',rLayer=kk)

    F = F/rg
    F[np.where(F>huge)] = np.nan
    F[np.where(F<0.001)] = 0.
    dH[kk-1,:,:] = F

  ZZ, ZM = mhycom.zz_zm_fromDP(dH)
  kdm = ZM.shape[0]
  jdm = ZM.shape[1]
  idm = ZM.shape[2]

  # Read T:
  A3d  = np.array([])
  fld  = 'temp'
  for lvl in range (1,kdm+1):
    F,n1,m1,l1 = mhycom.read_hycom(fina,finb,fld,rLayer=lvl)
    F[np.where(F>huge)] = np.nan
    if A3d.size == 0:
      A3d = F.copy()
      A3d = np.expand_dims(A3d, axis=0)
    else:
      F = np.expand_dims(F, axis=0)
      A3d = np.append(A3d, F, axis=0)

  # 2D Arrays with S/T above and below z0 for interpolation
  print('Searching top/btm values for {0}'.format(z0))
  Atop = np.zeros((jdm,idm))*np.nan
  Abtm = np.zeros((jdm,idm))*np.nan
  Ztop = np.zeros((jdm,idm))*np.nan
  Zbtm = np.zeros((jdm,idm))*np.nan
  for kk in range(kdm-1):
    zm1 = np.squeeze(ZM[kk,:,:])
    zm2 = np.squeeze(ZM[kk+1,:,:])
    if np.nanmin(zm2) > z0:
      continue
    if np.nanmax(zm1) < z0:
      continue

    print('lvl={0}  max(z1)/min(z2) = {1:5.1f}/{2:5.1f}'.\
           format(kk+1,np.nanmax(zm1),np.nanmin(zm2)))

    [J,I] = np.where((zm1 >= z0) & (zm2 < z0))
    Atop[J,I] = A3d[kk,J,I]
    Abtm[J,I] = A3d[kk+1,J,I]
    Ztop[J,I] = ZM[kk,J,I]
    Zbtm[J,I] = ZM[kk+1,J,I]

  # linear interp
  # Lagrange polynomial:
  Aint = Atop*((z0-Zbtm)/(Ztop-Zbtm)) + \
         Abtm*((z0-Ztop)/(Zbtm-Ztop))


  T12 = Aint.copy()
  # For continuous contour 
  # along the shelf fill shelf with cold water:
  T12[np.where( (HH<0) & (HH>=z0) & (LAT>24.2) )] = -2.
  T12[np.where(MS==0)] = np.nan
  T12[np.where((Msh==1) & (HH>z0))] = np.nan
  # Patches near Bahamas
  T12[1808:1815, 2571:2586] = np.nan
  TCNT = mgulf.derive_contour(T12, tz0=12.)


  # Read NAVO Gulf Stream northwall coordinates:
  # adjust # of days around rdate to look for NAVO
  # fronts using missing=... 
  importlib.reload(mrtofs)
  importlib.reload(mgulf)
  pthnavo = '/scratch2/NCEPDEV/ovp/Lichuan.Chen/DCOM/'
  XNW_navo, YNW_navo, navonmb = mgulf.read_navogs(rdate0, pthnavo, missing=1)
  if navonmb < 0.:
    print('NAVO not found skipping {0}/{1}/{2}\n'.format(YR,MM,DD))
    continue

  print('Mapping NAVO lon/lat --> RTOFS index space ...')
  INW_navo, JNW_navo = [], []
  knavo = len(XNW_navo)
  for ik in range(knavo):
    if ik%50 == 0:
      print('  {0:5.2f}% done ...'.format(float(ik)/float(knavo)*100.))

    xn0 = XNW_navo[ik]
    yn0 = YNW_navo[ik]

  #  inavo, jnavo = mutil.interp_indx_lonlat(xn0, yn0, LON, LAT)
    inavo, jnavo = mrtofs.interp_lonlat2indx(xn0, yn0, LON, LAT)
    INW_navo.append(inavo)
    JNW_navo.append(jnavo)

  INW_navo = np.array(INW_navo)
  JNW_navo = np.array(JNW_navo)

#
# Compute Modified Hausdorff Distance
# Truncate RTOFS contour to NAVO
# Make sure RTOFS contour is longer than NAVO
# Next, map RTOFS Gulf Stream indices --> lon/lat space
# Compute MHD
  INW_rtofs, JNW_rtofs = mgulf.adjust_rtofs_gsnw(TCNT[:,0], TCNT[:,1],\
                         INW_navo, JNW_navo, nskip=-1)

  print('Mapping RTOFS GS northwall index ---> Lon/Lat space ...')
  XNW_rtofs, YNW_rtofs = [], []
  krt = len(INW_rtofs)
  for ik in range(krt):
    if ik%50 == 0:
      print('  {0:5.2f}% done ...'.format(float(ik)/float(krt)*100.))

    xrt, yrt = mrtofs.interp_indx2lonlat(INW_rtofs[ik], JNW_rtofs[ik], LON, LAT)
    XNW_rtofs.append(xrt)
    YNW_rtofs.append(yrt)
                
  XNW_rtofs = np.array(XNW_rtofs)
  YNW_rtofs = np.array(YNW_rtofs)

  P = np.zeros((krt,2))
  Q = np.zeros((knavo,2))
  P[:,0] = XNW_rtofs
  P[:,1] = YNW_rtofs
  Q[:,0] = XNW_navo
  Q[:,1] = YNW_navo
  mhdGS = mmhd.modifHD(P, Q, geo2cart=True)
  print('{0} vs NEMO Gulf Stream MHD = {1:7.2f}km'.format(expt, mhdGS))

  MHDGS.append(mhdGS)
  TM.append(dnmb0)

  btx = 'gulfstream_stat.py'
  f_chck = False
  if f_chck:
    MHDinfo = 'Modif Hausdorff Dist RTOFS/NAVO = {0:5.1f} km'.format(mhdGS)
    ctl = '{7} SST GS front 12C/400m, {5} {1}/{2:02d}/{3:02d}:{4:02d}\n {6}'.\
           format(z0, YR, MM, DD, hr, fhcm, MHDinfo, expt)
    DV = mtime.datevec(navonmb)
    ssinf = 'RTOFS: {0}/{1:02d}/{2:02d}:{3:02d}\n'.format(YR,MM,DD,hr)
    ssinf = ssinf + 'NAVO: {0}/{1:02d}/{2:02d}:{3:02d}\n'.\
             format(DV[0],DV[1],DV[2],DV[3])

    mgulf.plot_gulfstr(TCNT, A3d, INW_navo, JNW_navo, LON, LAT, \
                       INW_rtofs, JNW_rtofs, sinfo=ctl, btx=btx, ssinf=ssinf)

  icnt += 1
  if icnt%10 == 0 and f_save > 0:
    print('Saving --> ' + fmhd_out)
    with open(fmhd_out,'wb') as fid:
      pickle.dump([np.array(MHDGS), np.array(TM)], fid)
      
  
print('Calculation of MHD completed')
MHDGS = np.array(MHDGS)
TM    = np.array(TM)

if f_save > 0:
  print('Saving --> ' + fmhd_out)
  with open(fmhd_out,'wb') as fid:
    pickle.dump([MHDGS,TM], fid)

print('    All Done')



