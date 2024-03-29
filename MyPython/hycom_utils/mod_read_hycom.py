"""
  Read HYCOM *.[ab] files
  grid, topo, 
  
  Dmitry Dukhovskoy, NOAA/NWS/NCEP/EMC
   
"""
import os
from netCDF4 import Dataset as ncFile
import numpy as np
import sys

sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython')
import mod_misc1 as mmsc1


def read_grid_topo(pthtopo,ftopo,fgrid,f_lon=180,pthgrid=None, dpth_neg=True):
# Get topo and lon/lat (plon/plat) for plotting 
# To get other lon/lat (ulon, ulat, vlon, vlat) - need to 
# read full grid file
#
# f_lon = 180 - convert lon to [-180, 180]
#       = 360 - convert lon to [0, 360]
#       = else - leave as in the grid file
#
# read HYCOM grid and topo files *.[ab]
  if pthgrid is None:
    pthgrid = pthtopo
  fltopoa = pthtopo+ftopo+'.a'
  fltopob = pthtopo+ftopo+'.b'
  flgrida = pthgrid+fgrid+'.a'
  flgridb = pthgrid+fgrid+'.b'

# Read dim, lon, lat
  fgb = open(flgridb,'r')
  fgb.seek(0)
  data = fgb.readline().split()
  IDM = int(data[0])
  data = fgb.readline().split()
  JDM = int(data[0])
  IJDM = IDM*JDM
  fgb.close()

  npad =4096-IJDM%4096

  print('Reading HYCOM grid/topo:{0} {1} '.format(fgrid,ftopo))
  print('Grid: IDM={0}, JDM={1}'.format(IDM,JDM))

# Read direct access binary file
  fga = open(flgrida,'rb')
  fga.seek(0)
  plon = np.fromfile(fga,dtype='>f',count=IJDM)

  if f_lon == 180:
    while any(plon > 180.):
      plon = np.where(plon<=180.,plon,plon-360.)

    while any(plon < -180.):
      plon = np.where(plon >= -180., plon, plon+360.)

  if f_lon == 360:
    while any(plon < 0.):
      plon = np.where(plon >= 0., plon, plon+360.)
    while any(plon > 360.):
      plon = np.where(plon <= 360., plon, plon-360.)

  plon = plon.reshape((JDM,IDM))


  fga.seek(4*(npad+IJDM),0)
  plat = np.fromfile(fga, dtype='>f',count=IJDM)
  plat = plat.reshape((JDM,IDM))

  fga.close()

  print('Min/max lon = {0}, {1}, Min/max lat = {2}, {3}'.format(np.min(plon),\
        np.max(plon),np.min(plat), np.max(plat)))

# Read bottom topography:
# Big endian float 32
  fbt = open(fltopoa,'rb')
  fbt.seek(0)
  HH = np.fromfile(fbt, dtype='>f', count=IJDM)
  HH = HH.reshape((JDM,IDM))
  fbt.close()

  if dpth_neg:
    HH[HH<1.e10] = -1.*HH[HH<1.e10]
    HH[HH>1.e10] = 100.

  print('Min/max Depth = {0}, {1}'.format(np.min(HH),np.max(HH)))

  return plon,plat,HH


def read_topo(pthtopo,ftopo,IDM,JDM,dpth_neg=True):
# Get topo only need to know I and J dimensions
#
# f_lon = 180 - convert lon to [-180, 180]
#       = 360 - convert lon to [0, 360]
#       = else - leave as in the grid file
#
# read HYCOM grid and topo files *.[ab]
  fltopoa = pthtopo+ftopo+'.a'
  fltopob = pthtopo+ftopo+'.b'

  IJDM = IDM*JDM
  npad =4096-IJDM%4096

  print('Reading HYCOM topo {0} '.format(ftopo))
  print('Grid: IDM={0}, JDM={1}'.format(IDM,JDM))

# Read bottom topography:
# Big endian float 32
  fbt = open(fltopoa,'rb')
  fbt.seek(0)
  HH = np.fromfile(fbt, dtype='>f', count=IJDM)
  HH = HH.reshape((JDM,IDM))
  fbt.close()

  if dpth_neg:
    HH[HH<1.e10] = -1.*HH[HH<1.e10]
    HH[HH>1.e10] = 100.

  print('Min/max Depth = {0}, {1}'.format(np.min(HH),np.max(HH)))

  return HH

def readnc_grid_topo(pthgrid,flgrid,f_lon=180,dpth_neg=True):
  """
  Read grid topo from NetCDF 
  If dpth_neg: make depths <0 and land 100 (default)
  If not dpth_neg: depths>0, land - huge 
  """
  from netCDF4 import Dataset as NetCDFFile

  hg = 1.e30
  nc = NetCDFFile(pthgrid+flgrid)
  lonc = nc.variables['Longitude'][:]
  latc = nc.variables['Latitude'][:]
  HH   = nc.variables['Bathymetry'][:]
  print('Reading ',pthgrid+flgrid)

  if f_lon == 180:
    while any(lonc.flatten() > 180.):
      lonc = np.where(lonc<=180.,lonc,lonc-360.)

    while any(lonc.flatten() < -180.):
      lonc = np.where(lonc >= -180., lonc, lonc+360.)

  if f_lon == 360:
    while any(lonc.flatten() < 0.):
      lonc = np.where(lonc >= 0., lonc, lonc+360.)
    while any(lonc.flatten() > 360.):
      lonc = np.where(lonc <= 360., lonc, lonc-360.)

  if dpth_neg:
    if np.min(HH) > 0.:
      HH[HH<1.e10] = -1.*HH[HH<1.e10]
      HH[HH>1.e10] = 100.
  else:
    if np.min(HH) < 0.:
      HH[HH>=0] = hg
      HH[HH<0]  = -1.*HH[HH<1.e10]

  print('Min/max lon = {0}, {1}, Min/max lat = {2}, {3}'.format(np.min(lonc),\
        np.max(lonc),np.min(latc), np.max(latc)))
  print('Min/max Depth = {0}, {1}'.format(np.min(HH),np.max(HH)))

  return lonc,latc,HH

####
def land3d_pmask(fina,finb, fland=0, fsea=1, fld='temp'):
  """
  Create 3D land mask for all model layers 
  for p-point variables
  default: land = 0, sea = 1
  """

  print('land3d_pmaks:  Creating 3D Land mask p-points')
  F, idm, jdm, kdm = read_hycom(fina,finb,'temp')
  F[F>1.e20] = np.nan

  Lmsk = np.zeros((jdm,idm,kdm))
  for kk in range(kdm):
    aa = np.squeeze(F[:,:,kk])
    aa = np.where(np.isnan(aa),fland,fsea)
    Lmsk[:,:,kk] = aa

  return Lmsk


def hycom_dim(fina,finb,fld='temp'):
  """
  Obtain hycom dimensions from *b file
  """
  try: 
    fgb = open(finb,'r')
  except:
    print('Could not open '+finb)

  fgb.seek(0)
  nl0 = 0
  while nl0 < 100:
    nl0 += 1
    data = fgb.readline().split()
    if len(data) < 2:
      continue
    adim = data[1]
    ii = adim.find('idm')
    if ii>0:
      break

  if ii<0:
    fgb.close()
    raise Exception('No idm found: Reading ',finb)

  IDM = int(data[0])
  data = fgb.readline().split()
  JDM = int(data[0])
  IJDM = IDM*JDM
#  fgb.seek(0)

  npad =4096-IJDM%4096

  aa = fgb.readline().split()

# Find location of the requested field:
  cntr= 0
  nf = len(fld)
  FLOC=[]
  while cntr<1e6:
    aa = fgb.readline().split()
    if len(aa) == 0:  # end of file
      break
    cntr += 1
    aname = aa[0]
    ii = aname.find(fld)
    if ii >= 0 and len(aname)==nf:
      FLOC.append(cntr)

  fgb.close()

  nrec = len(FLOC)
  if nrec == 0:
    raise Exception('read_hycom: Field {0} not found in {1}'.format(fld,finb))

# N. of v. layers
  KDM = len(FLOC)
  print('Grid: IDM={0}, JDM={1}, KDM={2}'.format(IDM,JDM,KDM))

  return IDM, JDM, KDM


def read_hycom(fina,finb,fld,Rtrc=None,rLayer=None,finfo=True):
  """
  reads hycom binary archive files (model output), 
  returns specified field 'fld'
  and dimensions of the grid
  Rtrc - passive tracer # to read, if there are tracers 
  rLayer - layer number to read, otherwise all layers are read - more time/memory
          numbering is lr=1,...,nlayers in the model
          rLayer is a layer number not an index, i.e.
          rLayer = 1 is 1st layer corresponds to index 0 in python

%  2017: added options for nlayer, n tracers
% if several tracers, options are:
% read tracer N: 'r_tracer',1
% read all tracers by default
% any variable can be read in 1 layer Nl:
%       'r_layer',1
%
% If all tracers are read, recommended to specify 
% only 1 layer to read 
%
  """
  try: 
    fgb = open(finb,'r')
  except:
    print('Could not open '+finb)

  fgb.seek(0)
  nl0 = 0
  while nl0 < 100:
    nl0 += 1
    data = fgb.readline().split()
    if len(data) < 2:
      continue
    adim = data[1]
    ii = adim.find('idm')
    if ii>0:
      break

  if ii<0:
    fgb.close()
    sys.exit('No idm found: Reading ',finb)

  IDM = int(data[0])
  data = fgb.readline().split()
  JDM = int(data[0])
  IJDM = IDM*JDM
#  fgb.seek(0)

  npad =4096-IJDM%4096

  if finfo:
    print('Reading HYCOM :{0} '.format(finb))

  aa = fgb.readline().split()

# Find location of the requested field:
  cntr= 0
  nf = len(fld)
  FLOC=[]
  while cntr<1e6:
    aa = fgb.readline().split()
    if len(aa) == 0:  # end of file
      break
    cntr += 1
    aname = aa[0]
    ii = aname.find(fld)
    if ii >= 0 and len(aname)==nf:
      FLOC.append(cntr)

  fgb.close()

  nrec = len(FLOC)
  if nrec == 0:
    raise Exception('read_hycom: Field {0} not found in {1}'.format(fld,finb))

# N. of v. layers
  """
   If fld = tracer and # of tracers >1
   need to distinguish # of layers 
   vs # of tracers*layers
   if strmatch(fld,'tracer')
  """
  ll = len(FLOC)
  if finfo:
    print('Grid: IDM={0}, JDM={1}, KDM={2}'.format(IDM,JDM,ll))

  FLOC = np.array(FLOC)
  if ll == 1:
    nVlev = 1
    nTR = 1
  else:
    dI = np.diff(FLOC)
    dmm = np.where(dI>1)
    dindx = dmm[0]
    nTR = dindx[0]+1         # N. of tracers in 1 layer
    nVlev = ll/nTR        # # of v. layers

# breakpoint()
  if nTR != 1 and finfo:
    print('read_hycom: Found {0} variables {1}  per layer'.format(nTR,fld))

# Find layers to read, if specified
# and tracers, if specified
  lr1=-1
  lr2=-1
  if Rtrc is not None:
    if nTR < Rtrc:
      raise Exception('Number of saved tracers {0} < requested {1}'.\
                      format(nTR,Rtrc))
    dmm = np.copy(FLOC)
    FLOC = dmm[Rtrc-1::nTR]

    if lr1 < 0 or lr2 < 0 :
      lr2 = FLOC.shape[0]
      ll = lr2

  if rLayer is not None:
    lr1 = rLayer
    lr2 = lr1

# If a layer Number not specified - read all
  if lr1 < 0 or lr2 < 0:
    lr1 = 1
    lr2 = ll

#  print('Reading {0}, Layers: {1}-{2}'.format(fld,lr1,lr2))
  fga = open(fina,'rb')
  F   = np.array([])
  ccL = -1
  huge = 0.0001*2.**99
  for ii in range(lr1,lr2+1):
    fga.seek(0)
    k0 = FLOC[ii-1]-1
    fga.seek(k0*(npad+IJDM)*4,0)
    dmm = np.fromfile(fga, dtype='>f',count=IJDM) # read 1 layer
    amin = np.min(dmm[np.where(dmm<huge)])
    amax = np.max(dmm[np.where(dmm<huge)])
    if finfo:
      print('Reading {0} k={1} min={2:12.4f} max={3:12.4f}'.\
            format(fld,ii,amin,amax))
    dmm = dmm.reshape((JDM,IDM))
    ccL += 1
#   print('lr={0}, ccL={1}'.format(ii,ccL))
#   breakpoint()
    if ccL == 0:
      F = np.copy(dmm)
      F = np.expand_dims(F, axis=0)
    else:
      dmm = np.expand_dims(dmm, axis=0)
#      print(dmm.shape)
#      print(F.shape)
      F = np.append(F, dmm, axis=0)

  if ll == 0:
    print('!!! read_hycom: {0} not found in {1} ERR'.format(fld,fina))
    print('!!! read hycom: check fields in ',finb)

  fga.close()

  return np.squeeze(F), IDM, JDM, ll


def read_hycom_restart(fina, finb, fld, IDM, JDM, \
                       Rtrc=None, rLayer=None, t_level=1, \
                       sprt='=', finfo=True):
  """
  reads hycom binary restart files (model output), 
  returns specified field 'fld' 
  if fields are saved at several (2) time levels
  only 1 (default = 1st) is read in
  IDM,JDM - J and I dimensions 
  sprt - separator symbol in *.b restart file 
         where numbers begin in the string, e.g.:
   temp    : layer,tlevel,range =  39 2  -2.3688760E+00  1.0018643E+02
 
  if several tracers, options are:
  read tracer N: 'r_tracer',1
  read all tracers by default
  any variable can be read in 1 layer Nl:
        'r_layer',1
 
  If all tracers are read, recommended to specify 
  only 1 layer to read 
 
  """
  
  IJDM = IDM*JDM
  npad = 4096-IJDM%4096
  toto = np.ones((npad,1),dtype=np.float32)

  try: 
    fgb = open(finb,'r')
  except:
    print('Could not open '+finb)

# Find all locations of required fields in *.b
  fgb.seek(0)
  nl0 = 0
  cntr=0      # poisition record # in restart
  nrec=-1     # conseq. # of the output field
  FLOC = np.array([], dtype=int)
  while nl0 >= 0:
    nl0 += 1
    data = fgb.readline().split()
    if len(data) == 0:
      break

    fldrd = data[0]
    if fldrd[0:7] == 'RESTART':    # skip 2-line header
      continue
    
    try:
      eqloc  = data.index(sprt)
      tm_loc = eqloc+2
      lr_loc = eqloc+1 
    except:
      raise Exception('ERR: separator in line {0} not found, check restart*.b'.\
                       format(sprt))

    cntr += 1
    tmlev  = int(data[tm_loc])
    lrnmb  = int(data[lr_loc])
    if data[0]==fld and tmlev==t_level:
      nrec += 1
      FLOC = np.append(FLOC,[cntr])
        
  if finfo:
    print('RESTART: {0}, time_level={1}, records found={2}'.\
         format(fld,t_level,nrec+1))

  fgb.close()

  if FLOC.shape[0] == 0:
    raise Exception('RESTART read failed: could not find var {0} in {1}'.\
                     format(fld,finb))

# N. of v. layers
  """
  If fld = tracer and # of tracers >1
  need to distinguish # of layers 
  vs # of tracers*layers
  if strmatch(fld,'tracer')
  """
  ll = FLOC.shape[0]
  if ll == 1:
    nVlev = 1
    nTR = 1
  else:
    dI = np.diff(FLOC)
    dmm = np.where(dI>0)
    dindx = dmm[0]
    nTR = dindx[0]+1         # N. of tracers in 1 layer
    nVlev = ll/nTR        # # of v. layers

# breakpoint()
  if nTR != 1:
    print('read_hycom_restart: Found {0} variables {1}  per layer'.\
          format(nTR,fld))
 
# Find layers to read, if specified
# and tracers, if specified
  lr1=-1
  lr2=-1
  if Rtrc is not None:
    if nTR < Rtrc:
      sys.exit('Number of saved tracers {0} < requested {1}'.format(nTR,Rtrc))
    dmm = np.copy(FLOC)
    FLOC = dmm[Rtrc-1::nTR]

    if lr1 < 0 or lr2 < 0 :
      lr2 = FLOC.shape[0]
      ll = lr2

  if rLayer is not None:
    lr1 = rLayer
    lr2 = lr1

# If a layer Number not specified - read all
  if lr1 < 0 or lr2 < 0:
    lr1 = 1
    lr2 = ll

  if finfo:
    print('Reading {0}, Layers: {1}-{2}'.format(fld,lr1,lr2))

  fga = open(fina,'rb')
  F = []
  ccL = -1
  for ii in range(lr1,lr2+1):
    fga.seek(0)
    k0 = FLOC[ii-1]-1
    fga.seek(k0*(npad+IJDM)*4,0)
    dmm = np.fromfile(fga, dtype='>f',count=IJDM) # read 1 layer
    dmm = dmm.reshape((JDM,IDM))
    ccL += 1
#   print('lr={0}, ccL={1}'.format(ii,ccL))
#   breakpoint()
    if ccL == 0:
      F = np.copy(dmm)
      F = np.expand_dims(F, axis=0)
    else:
#      F = np.dstack((F,dmm))
      dmm = np.expand_dims(dmm, axis=0)
      F = np.append(F, dmm, axis=0)

  if ll == 0:
    print('!!! read_hycom: {0} not found in {1} ERR'.format(fld,fina))
    print('!!! read hycom: check fields in ',finb)

  fga.close()

  return np.squeeze(F), ll 


def get_zz_zm(fina,finb,f_btm=True):
  """
    Derive layer depths: ZZ - interface depths,
    ZM - mid cell depths
    f_btm = true - ZZ=nan below bottom
  """
  rg = 9806.
  print(' Deriving ZZ, ZM from '+fina)

  fld = 'thknss'
  F, nn, mm, ll = read_hycom(fina,finb,fld,rLayer=1)
  F = F/rg
  F[np.where(F>1.e20)] = np.nan
  F[np.where(F<0.001)] = 0.
  ZZ = np.zeros((ll+1,mm,nn))
  ZZ[1,:,:] = -F
  ZM = np.zeros((ll,mm,nn))

  for kk in range(1,ll):
    F, nn, mm, ll = read_hycom(fina,finb,fld,rLayer=kk+1)
    F = F/rg
    F[np.where(F>1.e20)] = np.nan
    F[np.where(F<0.001)] = 0.
    ZZ[kk+1,:,:] = ZZ[kk,:,:]-F

# Land mask if requested:
    if f_btm:
      jb, ib = np.where(F<0.001)
      ZZ[kk+1,jb,ib] = np.nan
  
# Depths of the middle of the grid cells:
  for kk in range(ll):
    ZM[kk,:,:] = 0.5*(ZZ[kk+1,:,:]+ZZ[kk,:,:])

  return ZZ, ZM

def read_2Dbinary(fina,IDM,JDM,rLayer):
  """
    Read hycom binary 2D fields
    written with zaiowr hycom subroutine
    for 3D fields - read by layers
    only binary (*.a) file is needed
  """
  print('Reading {0}, Layer: {1}'.format(fina,rLayer))

  IJDM = IDM*JDM
  npad =4096-IJDM%4096
  fga = open(fina,'rb')
  fga.seek(0)
  k0 = rLayer-1
  fga.seek(k0*(npad+IJDM)*4,0)
  dmm = np.fromfile(fga, dtype='>f',count=IJDM) # read 1 layer
  dmm = dmm.reshape((JDM,IDM))
  F = np.copy(dmm)
#  F[np.where(F>1e25)]=np.nan

  fga.close()

  return F

def zz_zm_fromDP(dH, f_btm=True, finfo=True):
  """
    Calculate ZZ, ZM from layer thkcness (dH) 
    ZZ - interface depths,
    ZM - mid cell depths
    f_btm = true - nan below bottom & land
            false - no nans, all ZZ=zbtm or 0
  """
  print('Deriving ZZ, ZM from dH ...')
  ll = dH.shape[0]
  mm = dH.shape[1]
  nn = dH.shape[2]

  ZZ = np.zeros((ll+1,mm,nn))
  ZM = np.zeros((ll,mm,nn))
  for kk in range(ll):
    ZZ[kk+1,:,:] = ZZ[kk,:,:]-dH[kk,:,:]

    if f_btm:
      [JJ,II] = np.where(dH[kk,:,:] <= 1.e-30)
      ZZ[kk+1,JJ,II] = np.nan
      if kk==0: 
        ZZ[kk,JJ,II] = np.nan
    else:
      [JJ,II] = np.where(np.isnan(dH[kk,:,:]))
      ZZ[kk+1,JJ,II] = ZZ[kk,JJ,II]

# Depths of the middle of the grid cells:
    ZM[kk,:,:] = 0.5*(ZZ[kk+1,:,:]+ZZ[kk,:,:])
    if finfo:
      print(' kk={0} min/max ZZ = {1:5.1f}/{2:5.1f}'.\
           format(kk+1,np.nanmin(ZZ[kk,:,:]),np.nanmax(ZZ[kk,:,:])))

  return ZZ, ZM

def zz_zm_fromDP1D(dH, f_btm=True, finfo=True):
  """
    Calculate ZZ, ZM from layer thickness (dH) for 1D profile
    ZZ - interface depths,
    ZM - mid cell depths
    f_btm = true - nan below bottom & land
            false - no nans, all ZZ=zbtm or 0
  """
  if finfo:
    print('Deriving ZZ, ZM from dH ...')
  ll = dH.shape[0]

  ZZ = np.zeros((ll+1))
  ZM = np.zeros((ll))
  for kk in range(ll):
    ZZ[kk+1] = ZZ[kk]-dH[kk]

    if f_btm:
      if dH[kk] <= 1.e-30:
        ZZ[kk+1] = np.nan
        if kk==0: 
          ZZ[kk] = np.nan
    else:
      if dH[kk] <= 1.e-30:
        ZZ[kk+1] = ZZ[kk]

# Depths of the middle of the grid cells:
    ZM[kk] = 0.5*(ZZ[kk+1]+ZZ[kk])
    if finfo:
      print(' kk={0} ZZ(up/btm) = {1:5.1f}/{2:5.1f}, ZM = {3:5.1f}'.\
           format(kk+1,ZZ[kk],ZZ[kk+1],ZM[kk]))

  return ZZ, ZM

def dx_dy(LON,LAT):
  """
    Find horizontal grid spacing from LON,LAT 2D arrays 
    hycom grid
  """
  IDM = LON.shape[1]
  JDM = LON.shape[0]
  print('Calculating DX, DY for HYCOM idm={0}, jdm={1}'.format(IDM,JDM))
  
  DX = np.zeros((JDM,IDM))
  DY = np.zeros((JDM,IDM))

  for ii in range(IDM-1):
    LT1 = LAT[:,ii]
    LT2 = LAT[:,ii+1]
    LN1 = LON[:,ii]
    LN2 = LON[:,ii+1]
    dx  = mmsc1.dist_sphcrd(LT1,LN1,LT2,LN2)
    DX[:,ii] = dx

  DX[:,IDM-1] = dx

  for jj in range(JDM-1):
    LT1 = LAT[jj,:]
    LT2 = LAT[jj+1,:]
    LN1 = LON[jj,:]
    LN2 = LON[jj+1,:]
    dy  = mmsc1.dist_sphcrd(LT1,LN1,LT2,LN2)
    DY[jj,:] = dy
    
  DY[JDM-1,:] = dy

  print('Min/max DX, m = {0:12.4f} / {1:12.4f}'.format(np.min(DX),np.max(DX)))
  print('Min/max DY, m = {0:12.4f} / {1:12.4f}'.format(np.min(DY),np.max(DY)))

  return DX, DY
 

def find_indx_lonlat(x0,y0,LON,LAT):
  """
  Find closest grid point to lon/lat coordinate
  x0, y0 - geographic coordinates
  """
  if x0 > 180.:
    x0 = x0-360.

  dmm = np.sqrt((LON-x0)**2+(LAT-y0)**2)
  jj0, ii0 = np.where(dmm == np.min(dmm))

  return ii0[0], jj0[0]
  

def read_targ_dens(finb, ss, sep, ifld, kki=-1, prntTD=True):
  """
   Read target density from HYCOM *b file, finb
   ss is a piece of a string where to start reading targ densities 
   sep - separator character after which numeric values begin in the string
   ifld  - field # of the density in the numeric values after separator sep
   kki= -1: layer index is 1 before density index otherwise
           specify kki = field # after sep
  
  e.g.: relax_sal.b:
   PHC 3.0 (NOAA WOA-98 with improved Arctic) Climatology    
  Expt 06.0 nhybrd=32 nsigma=14 ds00= 0.50 dp00= 3.00 dp00x= 120.0 dp00f=1.180  
  Layered averages w.r.t. Sigma-2,  levtop=2                                     
  Salinity                                                                       
  i/jdm = 1600 2520                                                              
   sal: month,layer,dens,range = 01 01 28.100 3.7587984 38.609257   % 
   sal: month,layer,dens,range = 01 02 28.900 3.7587984 38.609257    
  
   ...
  %
   fnm = 'relax_sal.b'
   TDENS=read_targ_dens(fnm,'sal','=',3); 3rd value is density

# Allowed density deviation from the target density:
# see blkdat.input: sigjump 
# 0.02   'sigjmp' = minimum density jump across interfaces  (kg/m**3)

  """

  try:
    fgb = open(finb,'r')
  except:
    print('Could not open '+finb)

  print('Reading target densities ' + finb)
  if kki == -1:
    kfld = ifld-1 # layer # before targ dens
  else:
    kfld = kki 

  fgb.seek(0)
  nl0 = 0
  kl = -1  # python layer index for t.dens
  TDENS = []
  while nl0 < 10000:
    nl0 += 1
    data = fgb.readline().split()
    if len(data) < 2:
      continue
    afld = data[0]
    if not afld == ss:
      continue

    isub = data.index(sep)
    sstr = data[isub:]
    klr  = int(sstr[kfld])
    tdns = float(sstr[ifld])
#
# Some *b may have several fields for same layer , restart - 2 time steps
# check that layer # is updated
    if klr > 1 and TDENS[kl] == tdns:
      continue

    kl = kl+1
    TDENS.append(tdns)

    if prntTD:
      print('>> Lr {0:03d} = {1:8.4f}'.format(klr, tdns))

  fgb.close()

  TDENS = np.array(TDENS)

  if kl == 0:
    print('!!!!!!! No target densities found in {0} !!!'.format(finb))

  print('Read {0} Target Densities'.format(klr))

  return TDENS
    



