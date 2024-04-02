"""
  MOM6 utilities
  reading grid
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
from netCDF4 import Dataset as ncFile

def read_mom6grid(fgrid, grdpnt='hgrid', grid='nonsymmetr'):
  """
    Read MOM6 grid
    coordinates on "supergrid" (half grid points of actual grid)
    default - coordinates of h-points (scalar) are derived
    options: 'qgrid' - coordinated in the upper-right corner
                      q-points, vorticity points

   https://mom6.readthedocs.io/en/main/api/generated/pages/Discrete_Grids.html
 
   q(i-1,j)   ------------  v(i,j) --------------  q(i,j)
     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  

   u(i-1,j)               * h(i,j)                u(i,j)

     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  
   q(i-1,j-1) ------------  v(i,j-1) ------------  q(i,j-1)


  grid nonsymmetr: All arrays are declared with the same shape (isd:ied,jsd:jed)
                 There are fewer staggered variables to the south-west of the 
                 computational domain. An operator applied at h-point locations 
                 involving u- or v- point data can not have as wide a stencil on 
                 the south-west side of the processor domain as it can on 
                 the north-east side.

  grid symmetr:   Arrays have different shapes depending on their staggering 
                  location on the Arakawa C grid. 

  """

  print(f"Reading MOM6 grid for: {grid} grdpnt={grdpnt}")
  print('Grid: ' + fgrid)
  nc  = ncFile(fgrid,'r')
# Read MOM supergrid:
  XX  = nc.variables['x'][:].data
  YY  = nc.variables['y'][:].data
  mm  = XX.shape[0]
  nn  = XX.shape[1]


  if grid == 'nonsymmetr':
    if grdpnt == 'qgrid':
      LON = XX[2::2, 2::2]
      LAT = YY[2::2, 2::2]
    elif grdpnt == 'hgrid':
      LON = XX[1::2, 1::2]
      LAT = YY[1::2, 1::2]
    elif grdpnt == 'ugrid':
      LON = XX[1::2, 2::2]
      LAT = YY[1::2, 2::2]
    elif grdpnt == 'vgrid':
      LON = XX[2::2, 1::2]
      LAT = YY[2::2, 1::2]
  elif grid == 'symmetr':
    if grdpnt == 'qgrid':
      LON = XX[0::2, 0::2]
      LAT = YY[0::2, 0::2]
    elif grdpnt == 'hgrid':
      LON = XX[1::2, 1::2]
      LAT = YY[1::2, 1::2]
    elif grdpnt == 'ugrid':
      LON = XX[1::2, 0::2]
      LAT = YY[1::2, 0::2]
    elif grdpnt == 'vgrid':
      LON = XX[0::2, 1::2]
      LAT = YY[0::2, 1::2]

  jdm = LAT.shape[0]
  idm = LAT.shape[1]

  print(f" MOM6 {grdpnt} {grid} grid dim (J,I): {jdm} x {idm}")
  LON = np.where(LON < -180., LON+360., LON)

  return LON, LAT

def read_mom6angle(fgrid, grdpnt='hgrid', grid='nonsymmetr'):
  """
    Read angle the grid makes with the true East direction (lat line)
    radians
  """
  print(f"Reading MOM6 grid angle (radians) for: {grid} grdpnt={grdpnt}")
  print('Grid: ' + fgrid)
  nc  = ncFile(fgrid,'r')
# Read MOM supergrid:
  angle = nc.variables['angle_dx'][:].data
  mm  = angle.shape[0]
  nn  = angle.shape[1]

  if grid == 'nonsymmetr':
    if grdpnt == 'qgrid':
      alpha = angle[2::2, 2::2]
    elif grdpnt == 'hgrid':
      alpha = angle[1::2, 1::2]
    elif grdpnt == 'ugrid':
      alpha = angle[1::2, 2::2]
    elif grdpnt == 'vgrid':
      alpha = angle[2::2, 1::2]
  elif grid == 'symmetr':
    if grdpnt == 'qgrid':
      alpha = angle[0::2, 0::2]
    elif grdpnt == 'hgrid':
      alpha = angle[1::2, 1::2]
    elif grdpnt == 'ugrid':
      alpha = angle[1::2, 0::2]
    elif grdpnt == 'vgrid':
      alpha = angle[0::2, 1::2]

  jdm = alpha.shape[0]
  idm = alpha.shape[1]

  print(f" MOM6 dx angle {grdpnt} {grid} grid dim (J,I): {jdm} x {idm}")

  return  alpha

def read_mom6depth(ftopo, f_negate=True):
  """
    Read MOM6 depths at h-pnts
  """
  print('Reading MOM6 depths ' + ftopo)
  nc  = ncFile(ftopo,'r')
  HH  = nc.variables['depth'][:].data

# Convert depth to negatives:
# Land > 0
  if f_negate:
    HH  = np.where(HH > 0, -HH, 100.)

  return HH

def read_mom6lmask(ftopo):
  """
    Read MOM6 land mask
    1 = ocean, 0 = land
  """
  print('Reading MOM6 land mask ' + ftopo)
  nc   = ncFile(ftopo,'r')
  Lmsk = nc.variables['wet'][:].data

  return Lmsk


def read_mom6(foutp, fld, rLayer=-1, fnan=True, finfo=True):
  """
    Read 2D or 3D field from MOM6 output
    Specify Lyaer #  to read for 3D fields (from 1, ..., kdm)
    otherwise all layers are read in 3D array
    Replace out of range data with nan's - default

  """
  print('Reading MOM6 {0}: {1}'.format(fld, foutp))

  huge = 1.e18
  nc   = ncFile(foutp, 'r')

# Check variable dimension
# Assumed: AA(Time, dim1, dim2, optional: dim3)
  ndim = nc.variables[fld].ndim - 1 # Time excluded
  
# Check that time dim=1:
  tdim = nc.variables['Time'][:].data.shape[0]
  if tdim != 1:
    raise Exception('in netcdf output  Time dimensions is not 1: {0}'.\
                    format(tdim))

  if ndim == 2:
    FF = nc.variables[fld][:].squeeze().data
  elif ndim == 3:
    if rLayer <= 0:
      FF  = nc.variables[fld][:].squeeze().data
    else:
      klr = rLayer-1
      FF  = nc.variables[fld][:,klr,:,:].squeeze().data
  else:
    print('Output field dim should be 2 or 3')
    raise Exception('Check {0} dimension = {1}'.format(fld,ndim))

  AA = np.where(FF >= huge, np.nan, FF)
  if finfo:
    print('{0} lr {3} min/max: {1} / {2}'.\
         format(fld,np.nanmin(AA),np.nanmax(AA),rLayer))

  if fnan:
    FF = AA

  return FF

def mom_dim(flin):
  """
    Get horiz/vert dimensions for mom6
  """  
  nc   = ncFile(flin, 'r')
  xh = nc.variables['xh'][:].data.shape[0]
  yh = nc.variables['yh'][:].data.shape[0]
  zl = nc.variables['zl'][:].data.shape[0]

  return xh, yh, zl

def zz_zm_fromDP(dH, ssh, f_intrp=False, f_btm=True, finfo=True, eps0=1.e-5):
  """
    Calculate ZZ, ZM from layer thkcness (dH) 
    ZZ - interface depths,
    ZM - mid cell depths
    Note that sum(dH) = water column height that includes SSH
    Therefore ZZ[0] = ssh

    f_btm = true - nan below bottom & land
            false - no nans, all ZZ=zbtm or 0
    f_intrp = True: for vertical interpolation, land/bottom make
              not nans and not equal depths, i.e. depth continue
              increasing at eps0 rate below bottom to keep
              monotonicity of zz, zm
    f_interp overrides f_btm flag: makes it False
  """
  if f_intrp:
    f_btm = False

  print('Deriving ZZ, ZM from dH ...')
  ll = dH.shape[0]
  mm = dH.shape[1]
  nn = dH.shape[2]

  ZZ = np.zeros((ll+1,mm,nn))
  ZM = np.zeros((ll,mm,nn))
  ssh = np.where(np.isnan(ssh),0.,ssh)
  ZZ[0,:,:] = ssh
  if f_btm:
    dH = np.where(np.isnan(dH),0., dH)
  else:
    dH = np.where(np.isnan(dH),eps0, dH)

  if f_intrp:
    dH = np.where(dH < eps0, eps0, dH)

  for kk in range(ll):
    ZZ[kk+1,:,:] = ZZ[kk,:,:]-dH[kk,:,:]

    if f_btm:
      [JJ,II] = np.where(dH[kk,:,:] <= eps0)
      ZZ[kk+1,JJ,II] = np.nan
      if kk==0:
        ZZ[kk,JJ,II] = np.nan
    elif not f_btm and not f_intrp:
      [JJ,II] = np.where(dH[kk,:,:] <= eps0)
      ZZ[kk+1,JJ,II] = ZZ[kk,JJ,II]

# Depths of the middle of the grid cells:
    ZM[kk,:,:] = 0.5*(ZZ[kk+1,:,:]+ZZ[kk,:,:])
    if finfo:
      print(' kk={0} min/max ZZ = {1:5.1f}/{2:5.1f}'.\
           format(kk+1,np.nanmin(ZZ[kk,:,:]),np.nanmax(ZZ[kk,:,:])))

  return ZZ, ZM

def get_zz_zm(fina, f_btm=True):
  """
    Derive layer depths: ZZ - interface depths,
    ZM - mid cell depths
    Note that sum(dH) = water column height that includes SSH
    Therefore ZZ[0] = ssh
    f_btm = true - ZZ=nan below bottom
  """
  print(' Deriving ZZ, ZM from '+fina)
  huge = 1.e18
  nc   = ncFile(fina, 'r')

# Check variable dimension
# Assumed: AA(Time, dim1, dim2, optional: dim3)
  ndim = nc.variables['h'].ndim - 1 # Time excluded

# Check that time dim=1:
  tdim = nc.variables['Time'][:].data.shape[0]
  if tdim != 1:
    raise Exception('in netcdf output  Time dimensions is not 1: {0}'.\
                    format(tdim))

# Read layer thicknesses:
  dH     = nc.variables['h'][:].squeeze().data
  dH     = np.where(dH >= huge, np.nan, dH)
  ssh    = nc.variables['SSH'][:].squeeze().data
  ssh    = np.where(ssh >= huge, np.nan, ssh)
  ZZ, ZM = zz_zm_fromDP(dH, ssh, f_btm=f_btm)  

  return ZZ, ZM


def create_time_array(date1, date2, dday, date_mat=False):
  """
    Create time array for plotting fields
    if date_mat: date1 and date2 are in mat format
    otherwise: dates = [year, month, day]

  """
  import mod_time as mtime
  importlib.reload(mtime)

  if not date_mat:
#    yr1, mo1, mday1 = date1 
#    yr2, mo2, mday2 = date2  
#
    dnmb1 = mtime.datenum(date1)
    dnmb2 = mtime.datenum(date2)
  else:
    dnmb1 = date1
    dnmb2 = date2

  ldays = np.arange(dnmb1,dnmb2+dday/100.,dday)   
  nrec  = ldays.shape[0]
  TPLT  = np.zeros(nrec, dtype=[('dnmb', float),
                                ('date', int, (4,)),
                                ('yrday', float)])

  for irc in range(nrec):
    dnmb = ldays[irc]
    DV   = mtime.datevec(dnmb)
    _, jday = mtime.dnmb2jday(dnmb)

    TPLT['dnmb'][irc]   = dnmb
    TPLT['yrday'][irc]  = jday
    TPLT['date'][irc,:] = DV[0:4]

  print('Creating TPLT Time array, start: ' + \
        '{0}/{1}/{2} {3}hr, '.format(TPLT['date'][0,0], \
        TPLT['date'][0,1], TPLT['date'][0,2], TPLT['date'][0,3]) + \
        'end: {0}/{1}/{2} {3}hr'.\
      format(TPLT['date'][-1,0], TPLT['date'][-1,1], TPLT['date'][-1,2],\
      TPLT['date'][-1,3]))

  return TPLT  


def dx_dy(LON,LAT):
  """
    Find horizontal grid spacing from LON,LAT 2D arrays 
    hycom grid
  """
  import mod_misc1 as mmsc1

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

 
def vol_transp_2Dsection(LSgm, ZZ, UV):
  """
    Compute Vol transport along a straight line (EW or SN orientation)
    Note this algorithm works for zigzaging section
    as long as ZZ, UV are collocated in the center of the segments
    of the contour line
   
    Use midpoint quadrature
    All inputs are 2D fields:
    LSgm - segment lengths (grid cell dx)
    ZZ   - interface depths
    UV   - normal U component 
  """
  klv  = UV.shape[0]
  nsgm = UV.shape[1]
  dZ   = abs(np.diff(ZZ, axis=0))

  Ctrp = np.zeros((klv,nsgm))
  for kk in range(klv):
    Ctrp[kk,:] = UV[kk,:]*dZ[kk,:]*LSgm

  volTr_1D = np.nansum(Ctrp, axis=0)

  return volTr_1D

def fill_bottom(A2d, ZZ, Hbtm, fill_land=True):
  """
    Fill bottom values for plotting
    ZZ - layer interface depths 1D or 2D
    fill_land - fill land or not with closest values
  """
  A2df = A2d.copy()
  
  if len(ZZ.shape) == 2:
    f2d = True
  else:
    f2d = False
    z0 = ZZ

  kdm = A2d.shape[0]
  idm = A2d.shape[1]
  iOc = np.where(Hbtm < 0.)[0]

  for ii in range(idm):
    hb0 = Hbtm[ii]
    if hb0 >= 0.:
      if fill_land:
        dii        = abs(iOc - ii)
        ix         = np.where(dii == np.min(dii))[0][0]
        iocn       = iOc[ix]
        A2df[:,ii] = A2d[:,iocn] 
        continue
      else:
        continue 

    if f2d:
      z0 = ZZ[:,ii]

    dZb = ZZ-hb0
    izb = min(np.where(dZb <= 1.e-3)[0])
    A2df[izb:,ii] = A2d[izb-1,ii]
    
  return A2df

def bottom2nan(A2d, ZZ, Hbtm):
  """
    Fill bottom/land values with nans for plotting
    ZZ - layer interface depths 1D or 2D
  """
  A2df = A2d.copy()
  
  if len(ZZ.shape) == 2:
    f2d = True
  else:
    f2d = False
    z0 = ZZ

  kdm = A2d.shape[0]
  idm = A2d.shape[1]

  for ii in range(idm):
    hb0 = Hbtm[ii]
    if hb0 >= 0.:
     A2df[:,ii] = np.nan

    if f2d:
      z0 = ZZ[:,ii]

    dZb = ZZ-hb0
    izb = min(np.where(dZb <= 1.e-3)[0])
    A2df[izb:,ii] = np.nan
    
  return A2df


