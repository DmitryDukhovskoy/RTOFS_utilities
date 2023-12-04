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


