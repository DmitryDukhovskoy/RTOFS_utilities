"""
  Utility functions for reading/plotting GDEM climatology
  How to read gdem .short
  see code:
  /scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/ncoda/sorc/ncoda_qc/libsrc/prfobs
  gdem_mod.f

"""
import numpy as np

def gdem_grid():
  """
    GDEM parameters
  """
  NJ   = 689
  NI   = 1441
  NK   = 78
  lat1 =  -82.
  lon1 = 0.
  dx   = 0.25

  LON = np.linspace(lon1, 360., num=NI, endpoint=True)
  LAT = np.linspace(lat1, 90., num=NJ, endpoint=True)

  ZM = [0.,    2.,    4.,    6.,    8.,   10.,   15.,\
       20.,   25.,   30.,   35.,   40.,   45.,   50.,\
       55.,   60.,   65.,   70.,   75.,   80.,   85.,\
       90.,   95.,  100.,  110.,  120.,  130.,  140.,\
      150.,  160.,  170.,  180.,  190.,  200.,  220.,\
      240.,  260.,  280.,  300.,  350.,  400.,  500.,\
      600.,  700.,  800.,  900., 1000., 1100., 1200.,\
     1300., 1400., 1500., 1600., 1800., 2000., 2200.,\
     2400., 2600., 2800., 3000., 3200., 3400., 3600.,\
     3800., 4000., 4200., 4400., 4600., 4800., 5000.,\
     5200., 5400., 5600., 5800., 6000., 6200., 6400.,\
     6600. ]

  ZM = -np.array(ZM)

  return LON, LAT, ZM

def read_gdem_Lmask():
  """
    Land mask
  """
  import pickle
  pthout   =  '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/GDEM/'
  gdem_out = pthout + 'GDEM_sealand_mask.pkl'
  print('Loading GDEM land mask --> ' + gdem_out)

  with open(gdem_out,'rb') as fid:
    Lmsk = pickle.load(fid)

  return Lmsk

def read_gdem3d(pthgdem, month, fldrd, NI, NJ, NK):
  """
    Read GDEM binary 3D fields
    fldrd - field to read: salt temp
    month - month to read
    NI, NJ, NK - dimensions of the 3D array
    in this GDEM data are short integers with offset and scale
  """
  import struct

  NIJK = NI*NJ*NK
  fina = '{0}gdem_{1}{2:02d}.short'.format(pthgdem,fldrd,month)
  print('Reading GDEM from ' + fina)
 
  try:
    fga = open(fina,'rb')
  except:
    print('Could not open '+finb)

  fga.seek(0)
#  print(struct.unpack('i', fga.read(4)))
  fdump  = np.fromfile(fga, dtype='>i4', count=1)[0]
  offst  = np.fromfile(fga, dtype='>f', count=1)[0]
  scale  = np.fromfile(fga, dtype='>f', count=1)[0]
  mssval = np.fromfile(fga, dtype='>i4', count=1)[0]
  fdump2  = np.fromfile(fga, dtype='>i4', count=1)[0]
#  print('Record {0} Bt'.format(fdump))  
#  print('offst={0}'.format(offst))  
#  print('scale={0}'.format(scale))  
#  print('mssval={0}'.format(mssval))  
#  print('Written {0} Bt'.format(fdump2))  

  fdump  = np.fromfile(fga, dtype='>i4', count=1)[0]
  DATA = np.fromfile(fga, dtype='>i2', count=NIJK)
 
  fga.close()   

  DATA.astype(float)
  DATA = DATA*scale + offst
#  DATA[np.where(DATA<0.)] = np.nan

# AA = np.zeros((NK,NJ,NI))*np.nan
# icc = -1
# for kk in range(NK):
#   for jj in range(NJ):
#     for ii in range(NI):
#       icc = icc+1
#       AA[kk,jj,ii] = DATA[icc]

  AA = np.reshape(DATA,(NK,NJ,NI), order='C')

#  aa=np.squeeze(AA[0,:,:])
#  plt.clf()
#  plt.pcolormesh(aa)

  return AA

def gdem_profile(pthgdem, month, fldrd, lon0, lat0):
  """
  Extract GDEM profile at x=lon0, y=lat0
  """
  LONG, LATG, ZMG = gdem_grid()
  NI = LONG.shape[0]
  NJ = LATG.shape[0]
  NK = ZMG.shape[0]

  A3d = read_gdem3d(pthgdem, month, fldrd, NI, NJ, NK)
  if lon0 < 0.:
    II = np.where(LONG>180.)[0]
    LONG[II] = LONG[II]-360.

  dst = np.square((LONG-lon0)**2)
  I0  = np.where(dst == np.min(dst))[0][0]
  lng = LONG[I0]
  dst = np.square((LATG-lat0)**2)
  J0  = np.where(dst == np.min(dst))[0][0]
  ltg = LATG[J0]
  dst = np.square((lng-lon0)**2+(ltg-lat0)**2)
  if dst > 0.25:
    print(' ERROR finding closest GDEM pnt for {0:6.3f}N {0:6.3f}E'.\
           format(lat0, lon0))
    raise Exception('Check closest point in GDEM: {0:6.3f}N {0:6.3f}E'.\
           format(ltg, lng))

  Aprf = np.squeeze(A3d[:,J0,I0])

  return ZMG, Aprf

def find_gdem_indx(XX, YY, LON, LAT, Hb0 = []):
  """
    Find GDEM indices corresponding to X, Y coordinates
    LON, LAT - GDEM coordinates 1D arrays

  """
  import mod_misc1 as mmisc

  nI   = len(XX)
  II   = []
  JJ   = []
  XXG  = []
  YYG  = []
  Hbtm = []
  nLN  = len(LON)
  nLT  = len(LAT)
  LAT  = np.where(LAT>=90., 89.99, LAT)

  print(f"Finding section {nI} indices on GDEM grid ...")
  for ix in range(nI):
    if ix%200 == 0:
      npr = ix/nI*100
      print(f"   {npr:.2f}% done ...")

    x0 = XX[ix]
    y0 = YY[ix]

    DMIN = np.zeros((nLT))-999.
    IMIN = np.zeros((nLT))-999.
    for jj in range(nLT):
      ltt = np.zeros((nLN)) + LAT[jj]
      dst = mmisc.dist_sphcrd(y0, x0, ltt, LON)

      imn      = np.argmin(dst)
      DMIN[jj] = dst[imn]
      IMIN[jj] = imn

    jj0 = int(np.argmin(DMIN))
    ii0 = int(IMIN[jj0])
    if len(II) == 0:
      II.append(ii0)
      JJ.append(jj0)
      XXG.append(LON[ii0])
      YYG.append(LAT[jj0])
      if len(Hb0) > 0:
        Hbtm.append(Hb0[ix])
    else:
# check duplicates:
      if not (II[-1] == ii0 and JJ[-1] == jj0):
        II.append(ii0)
        JJ.append(jj0)
        XXG.append(LON[ii0])
        YYG.append(LAT[jj0])
# Extract bottom:
        if len(Hb0) > 0:
          Hbtm.append(Hb0[ix])
        
  II   = np.array(II, dtype=int)
  JJ   = np.array(JJ, dtype=int) 
  XXG  = np.array(XXG)
  YYG  = np.array(YYG) 
  Hbtm = np.array(Hbtm)

  nIX = len(II)
  Lsgm = np.zeros((nIX))
  for ii in range(nIX-1):
    dst = mmisc.dist_sphcrd(YYG[ii],XXG[ii],YYG[ii+1],XXG[ii+1])
    Lsgm[ii] = dst

  Lsgm[-1] = dst

  return II, JJ, XXG, YYG, Lsgm, Hbtm







