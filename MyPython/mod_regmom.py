"""
  Utility for regional MOM6-SIS2 
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pdb
import importlib
#import struct
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
from mod_utils_fig import bottom_text

def find_gridpnts_box(x0, y0, LON, LAT, dhstep=0.5):
  """
    Given pnt (x0,y0) find 4 grid points enclosing the pnt
    on a grid XX, YY - coordinates

     * -----------------  *
     |                    |
     |                    |
     |    *(x0,y0)        |
     |                    |
     |                    |
     |                    |
     * -----------------  *

    The algorithms works for both regular (straight orthogonal Mercator-type grids)
    and irregular (e.g. bipolar with slented grid lines) although in discontinuity regions
    where grid boxes are too narrow and slented may give errors 

    Specify offset for i,j indices if a subset XX, YY is ued
    instead of whole domain X and Y (to speed up distance calculation)
    to get the global indices 

    dhstep - approximate max grid horizontal stepping (in units of XX/YY)
             in order to accelerate searching algorithm 
             making it too small may result in errors as
             grid points outside this range will be discarded 
             
  """
  import mod_misc1 as mmisc1
  import mod_bilinear as mblnr

# Need 2D arrays for LON, LAT
# if 1D array - Mercator grid is assumed
  ndim = len(LON.shape)
  if ndim == 1:
    mm = len(LAT)
    nn = len(LON)
    LON = np.tile(LON,(mm,1))
    LAT = np.tile(LAT,(nn,1)).transpose()

  if x0 < 0.: x0 = x0 + 360.
  LON = np.where(LON < 0., LON+360., LON)
  assert np.max(LON) <=360.
  assert np.min(LON) >= 0.
  assert 0. < x0 < 360.

  mm, nn = LON.shape
# Subsample the region:
  dy = dhstep
  JJ,II = np.where((LAT > y0-dy) & (LAT < y0+dy))

  def find_closest_point(y0, x0, LON, LAT, JJ, II):
#      Find closest point within a subset of points given JJ, II indices
    XX = LON[JJ,II]
    YY = LAT[JJ,II]
    DD   = mmisc1.dist_sphcrd(y0,x0,YY,XX)
    kmin = np.argmin(DD)
#  jmin, imin = np.unravel_index(kmin, LON.shape)
#  jmin, imin = np.where(DD == np.min(DD))
    jmin = JJ[kmin]
    imin = II[kmin]
    xmin = LON[jmin, imin]
    ymin = LAT[jmin, imin]

    return jmin, imin, xmin, ymin

  jv1, iv1, xv1, yv1   = find_closest_point(y0, x0, LON, LAT, JJ, II) 
  xref = xv1-0.1
  yref = yv1-0.1
  xv1c, yv1c = mblnr.lonlat2xy_pnt(xv1, yv1, xref, yref)
  x0c, y0c   = mblnr.lonlat2xy_pnt(x0, y0, xref, yref)
# Find grid points around x0,y0:
# First guess for xv2: 
# start by finding 2nd point by moving +/- 1 grid point along I axis
  di  = np.sign(x0c-xv1c)
  din = di
  icc = 0 
  iv2 = iv1
  jv2 = jv1
  while din == di:
    iv2  = int(iv2 + di)
    xv2  = LON[jv2,iv2]
    yv2  = LAT[jv2,iv2]
    xv2c, yv2c = mblnr.lonlat2xy_pnt(xv2, yv2, xref, yref)
    din = np.sign(x0c-xv2c)
    icc += 1
#    if abs(din) > abs(di) and np.sign(din) == np.sign(di):
#      raise Exception(f"x0={x0} y0={y0} search direction wrong di={di} din={din}")
    if icc > 10:
      print(f"ERR: x0={x0} y0={y0} search direction wrong di={di} din={din}")
      raise Exception(f"cannot find enclosed box, icc={icc} tries")

# Define pnt orientation wrt to this line to find which 
# direction to search for another vertex
# along J axis
# if at the bndry - ignore, nothing can do
  Dxv2 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv2c,yv2c],[x0c,y0c]))
  F_bndry = False
  if jv1+1 >= mm:
    Djp1 = - Dxv2
    F_bndry = True
  else:
    y1p1 = LAT[jv1+1, iv1]
    x1p1 = LON[jv1+1, iv1]
    x1p1c, y1p1c = mblnr.lonlat2xy_pnt(x1p1, y1p1, xref, yref)
    Djp1 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv2c,yv2c],[x1p1c,y1p1c]))

# Singularity points:
# check if the grid points are same locations for varying j or i
# these points are over land - can be ignored
    D = mmisc1.dist_sphcrd(yv1, xv1, y1p1, x1p1)
    if D < 1.e-3:
      F_bndry = True
      Djp1 = - Dxv2
      
  if Djp1 != Dxv2:
# Try opposite direction from pnt 1 along J:
# Ignore the bndry
    if jv1-1 < 0:
      dirj = 1
      F_bndry = True
    else:
      y1m1 = LAT[jv1-1, iv1]
      x1m1 = LON[jv1-1, iv1]
      x1m1c, y1m1c = mblnr.lonlat2xy_pnt(x1m1, y1m1, xref, yref)
      Djm1 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv2c,yv2c],[x1m1c,y1m1c]))

      if Djm1 != Dxv2 and not F_bndry:
        print(f"ERR: Cannot find direction to include pnt {x0} {y0}")
        raise Exception(f"Unexpected: orientations same Djp1={Djp1} Djm1={Djm1}")

      dirj = -1
  else:
    dirj = 1

# Find a 3rd vertex along J-axis moving in direction = dirj 
# such that the pnt (x0,y0) projected onto J falls within the 
# segment (vrtx1 - vrtx3)
  iv3 = iv1
  jv3 = jv1 + dirj
  xv3 = LON[jv3,iv3]
  yv3 = LAT[jv3,iv3]
  xv3c, yv3c = mblnr.lonlat2xy_pnt(xv3, yv3, xref, yref)
  Bv  = mmisc1.construct_vector([xv1c, yv1c],[xv3c, yv3c]) 
  Av  = mmisc1.construct_vector([xv1c, yv1c],[x0c, y0c])
  prA, cosT, tht = mmisc1.vector_projection(Av, Bv)
  lBv = np.sqrt(np.dot(Bv.transpose(), Bv))[0][0]
  prA = abs(prA)

  icc = 0
  while prA > lBv:
# if projection is outside then extend the side to include the point in the box
# may run ouside the domain!
    jv3 = jv3 + dirj
    if jv3 >= mm or jv3 < 0: 
      F_bndry = True
      jv3 = jv3 - dirj
      break

    xv3 = LON[jv3,iv3]
    yv3 = LAT[jv3,iv3]
    xv3c, yv3c = mblnr.lonlat2xy_pnt(xv3, yv3, xref, yref)
    Bv  = mmisc1.construct_vector([xv1c, yv1c],[xv3c, yv3c]) 
    prA, cosT, tht = mmisc1.vector_projection(Av, Bv)
    lBv_new = np.sqrt(np.dot(Bv.transpose(), Bv))[0][0]
    if lBv_new < lBv:
      print(f"ERR: unexpected side xv3 not increasing, new: {lBv_new} old: {lBv}")
      raise Exception("unexpected edge length decreasing when choosing xv3")

    icc += 1
    if icc > 10:
      raise Exception(f"cannot find jv3 # steps={icc}")    

#ax1.plot(x0c,y0c,'r*')
#ax1.plot(xv1c,yv1c,'ro')
#ax1.plot(xv2c,yv2c,'go')
#ax1.plot(xv3c,yv3c,'bo')

# Check if first guess of xv2, yv2 was correct to include the point x0, y0:
  Dv2 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv3c,yv3c],[xv2c,yv2c]))
  Dp0 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv3c,yv3c],[x0c,y0c]))
  if Dp0 == 0: Dp0 = Dv2 # allow pnt to be on the line
  if Dv2 != Dp0:
# need to switch vertex x2 to opposite
# only one of those should be not 0:
    diri = np.sign(iv1-iv2)
    dirj = np.sign(jv1-jv2)

    icc = 0    
    while Dv2 != Dp0:
      iv2 = int(iv2 + diri)
      jv2 = int(jv2 + dirj)
      if iv2 < 0 or iv2 >= nn or jv2 < 0 or jv2 >= mm:
        F_bndry = True
        break
      xv2 = LON[jv2,iv2]
      yv2 = LAT[jv2,iv2]
      xv2c, yv2c = mblnr.lonlat2xy_pnt(xv2, yv2, xref, yref)
      Dv2 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv3c,yv3c],[xv2c,yv2c]))
  
      icc += 1 
      if icc > 10:
        raise Excpetion(f"adjusting vx2: too many steps {icc}")

# Vertx 4: best guess use vertex 2 and vertex 3:
  iv4 = iv2
  jv4 = jv3

# Construct a box enclosing the point x0, y0 and check: 
  IV = [iv1, iv2, iv4, iv3]
  JV = [jv1, jv2, jv4, jv3]
#  IV, JV = mblnr.sort_gridcell_indx(IV, JV)
#  IV0 = IV.astype(int)
#  JV0 = JV.astype(int)
  IV0 = IV.copy()
  JV0 = JV.copy()

# Convert box vertices into Cartesian coord wrt a reference pnt = centroid
  XX  = LON[JV,IV]
  YY  = LAT[JV,IV]
  XXc = 0.25*np.sum(XX)
  YYc = 0.25*np.sum(YY)
  XV, YV  = mblnr.lonlat2xy_wrtX0(XX, YY, XXc, YYc)
  x0c, y0c = mblnr.lonlat2xy_pnt(x0,y0, XXc, YYc)
  INp     = mmisc1.inpolygon_1pnt(x0c, y0c, XV, YV)
  INp2    = mmisc1.inpolygon_1pnt(x0, y0, XX, YY)

  if F_bndry:
    if not INp:
      print(f"WARNING: x0={x0:6.2f} y0={y0:6.2f} on " +\
            f"bndry of outerdomain, cannot enclose")
    ixx = np.array(IV).astype(int)
    jxx = np.array(JV).astype(int)
    return ixx, jxx

# Most of cases should be handled by the above algorithm
# for very non-orthogonal I/J grid lines near singularities in 
# tripolar grids the point may still be outside the box
# Try this:
# Adjust sides to include the pnt if it is not in the box
  SIDE2V = {
     "side1": [0, 1],
     "side2": [1, 2],
     "side3": [2, 3],
     "side4": [3, 0]
  }

  V2SIDE = {
     "vertex0": [4, 1],
     "vertex1": [1, 2],
     "vertex2": [2, 3],
     "vertex3": [3, 4]
  }
     
# Plot box:
# ax1.plot(XV,YV,'.-')    
# ax1.plot(0,0,'r*')    # pnt should be inside XV,YV
# Check the grid points, if the grid is not rectilinear a regular boxes may not work
# find centroid for a simple polygon, and check x0,y0 positioning compared to the centroid
# should be on the same side for all sides of the box
# At least one side should give different orientations of centroid and the point x0, y0
  icc = 0
  while not INp:
    xC = 0.25*(np.sum(XV))
    yC = 0.25*(np.sum(YV))
    for iside in range(1,5):
      k1, k2 = SIDE2V[f"side{iside}"]  
      x1 = XV[k1]
      x2 = XV[k2]
      y1 = YV[k1]
      y2 = YV[k2]
      iv1 = IV[k1]
      iv2 = IV[k2]
      jv1 = JV[k1]
      jv2 = JV[k2]
      xg1 = XX[k1]
      yg1 = YY[k1]
      xg2 = XX[k2]
      yg2 = YY[k2]

      Dcentr = np.sign(mmisc1.orientation([x1,y1],[x2,y2],[xC,yC]))
      Dpnt0  = np.sign(mmisc1.orientation([x1,y1],[x2,y2],[x0c,y0c]))
#      print(f"side {iside} Dcentr={Dcentr} Dpnt0={Dpnt0}")
      if Dpnt0 == 0: Dpnt0 = Dcentr  # Allow point to be on the edge
      if Dpnt0 == Dcentr: continue   # this side is ok
# Find which vertex to move by checking which one is opposite to x0,y0 and centroid
      if Dcentr != Dpnt0:
        dx1xc = np.sign(xC-x1)
        dx1x0 = np.sign(x0c-x1)
        dy1yc = np.sign(yC-y1)
        dy1y0 = np.sign(y0c-y1)
        dx2xc = np.sign(xC-x2)
        dx2x0 = np.sign(x0c-x2)
        dy2yc = np.sign(yC-y2)
        dy2y0 = np.sign(y0c-y2)
        div1 = djv1 = div2 = djv2 = 0
        if dx1xc != dx1x0: div1 = dx1x0
        if dx2xc != dx2x0: div2 = dx2x0
        if dy1yc != dy1y0: djv1 = dy1y0
        if dy2yc != dy2y0: djv2 = dy2y0

        icc = 0
        while Dcentr != Dpnt0:
          iv1 = int(iv1+div1)
          iv2 = int(iv2+div2)
          jv1 = int(jv1+djv1)
          jv2 = int(jv2+djv2)

          IV[k1] = iv1
          IV[k2] = iv2
          JV[k1] = jv1
          JV[k2] = jv2
          XX  = LON[JV,IV]
          YY  = LAT[JV,IV]
          xC = 0.25*np.sum(XX)
          yC = 0.25*np.sum(YY)
          XV, YV  = mblnr.lonlat2xy_wrtX0(XX, YY, xC, yC)
          x0c, y0c = mblnr.lonlat2xy_pnt(x0,y0, xC, yC)
          INp     = mmisc1.inpolygon_1pnt(x0c, y0c, XV, YV)

          if INp: break

          x1 = XV[k1]
          x2 = XV[k2]
          y1 = YV[k1]
          y2 = YV[k2]
          Dcentr = np.sign(mmisc1.orientation([x1,y1],[x2,y2],[xC,yC]))
          Dpnt0  = np.sign(mmisc1.orientation([x1,y1],[x2,y2],[x0c,y0c]))

          icc += 1
          if icc > 10: 
            raise Exception(f"Couldnot adjust side to include point #iter={icc}")


    
  ixx = np.array(IV).astype(int)
  jxx = np.array(JV).astype(int)
  
  return ixx, jxx

def fill_land(aa1,aa2,aa3,aa4,HH,A3d,JJ,II,Jocn,Iocn):
  """
    Fill all nans in 1D arrays - land points
    at least 1 array should be ocean pnt
  """
  la1  = len(np.where(np.isnan(aa1))[0])
  la2  = len(np.where(np.isnan(aa2))[0])
  la3  = len(np.where(np.isnan(aa3))[0])
  la4  = len(np.where(np.isnan(aa4))[0])
  kdmh = len(aa1)

  LL  = np.array([la1,la2,la3,la4])
  iL  = np.where(LL == kdmh)[0]
  inL = np.where(LL < kdmh)[0] 

# This should not happen - no land points
# assumed at least 1 array is nans
  if len(iL) == 0:
    return aa1, aa2, aa3, aa4

  AP  = np.column_stack((aa1,aa2,aa3,aa4))
# All land - find closes ocean point
  if len(iL) == 4:
#    raise Exception("All 4 1D arrays - land points, check Land masks")
    for ild in range(4):
      il1 = II[ild]
      jl1 = JJ[ild]    
      D   = np.square((Iocn-float(il1))**2 + (Jocn-float(jl1))**2)
      ixx = np.argmin(D)
      ioc = Iocn[ixx]
      joc = Jocn[ixx]
      AP[:, ild] = A3d[:,joc,ioc]

  APm = np.nanmean(AP, axis=1) 
  AP[:,iL] = APm[:,None]
 
  return AP[:,0].squeeze(), AP[:,1].squeeze(),\
         AP[:,2].squeeze(), AP[:,3].squeeze() 

def check_bottom(AA):
  """
    Make sure there are no NaNs in the 1D vertical data
    AA is 1D array
    If HYCOM land mask does not match MOM's --> all nans in the 1D profile
    from HYCOM
  """
  if AA[0] == np.nan:
    print(f"1D profile: all values are nans, Land mask mismatch")
    return 

  izb = np.argwhere(AA == np.nan)
  if len(izb) == 0:
    return

  AA[izb] = AA[min(izb)-1]

  return AA
  
def derive_TSprof_WOA23(seas, YR, Xp, Yp, grd=0.25, conv2pot=True):
  """
    Extract T & S profiles from WOA23
    for specified locations Xp, Yp
    Note: T - in situ T's
    conv2pot: Convert in situ T to potential 
    WOA:
    # season: 1-12 monthly, 13-winter (Jan-Mar), 14-spring (Apr-Jun), ...
    https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/bin/woa23.pl
    woa23_[DECA]_[v][tp]_[gr].nc - NetCDF format
    where:
    [DECA] - decade
    [v] - variable
    [tp] - time period
    [ft] - field type
    [gr] - grid
    e.g: woa23_95A4_t14_04.nc
  """
  import mod_swstate as msw

  if grd==0.25:
    cgrd=4
  woa='woa23'

  if YR > 1990 and YR < 2005:
    deca = "95A4"
  elif YR > 2005 and YR < 2014:
    deca = "A5B4"
  elif YR > 2014 and YR < 2023:
    deca = "B5C2"

  urlT = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/' + \
         f'temperature/netcdf/{deca}/0.25/'
  urlS = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/' + \
         f'salinity/netcdf/{deca}/0.25/'

  tfnm=f'{woa}_{deca}_t{seas:02d}_{cgrd:02d}.nc'
  sfnm=f'{woa}_{deca}_s{seas:02d}_{cgrd:02d}.nc'
#  print("Reading {1} from {0}".format(furl,varnm))
  
  tvar='t_an'
  svar='s_an'

  furl = os.path.join(urlT, tfnm)
  print(f"Extracting T prof from WOA23 decadal clim season={seas}")
  print(furl)
  nc   = ncFile(furl)
  T3d  = nc.variables[tvar][:].data.squeeze()
  T3d  = np.where(T3d>1.e10, np.nan, T3d)
  kdm, jdm, idm  = T3d.shape
 
  furl = os.path.join(urlS, sfnm)
  print(f"Extracting S prof from WOA23 decadal clim season={seas}")
  print(furl)
  nc   = ncFile(furl)
  S3d  = nc.variables[svar][:].data.squeeze()
  S3d  = np.where(S3d>1.e10, np.nan, S3d)

  ZZ   = nc.variables['depth'][:].data.squeeze()
  ZZ   = -abs(ZZ)
  lat  = nc.variables['lat'][:].data.squeeze()
  lon  = nc.variables['lon'][:].data.squeeze() 
  dlat = np.max(np.diff(lat))
  dlon = np.max(np.diff(lon))
 
  Npnts = len(Xp)
  print(f"N original locations={Npnts}, conv to potent temp={conv2pot}")

  if max(lon) <= 180. and min(lon) < 0.:
    Xp = np.where(Xp > 180., Xp-360., Xp)
 
  Iwoa = []
  Jwoa = []
  Tprf = np.zeros((kdm,1))
  Sprf = np.zeros((kdm,1))
  for kk in range(Npnts):
    x0 = Xp[kk]
    y0 = Yp[kk]

    dx = np.sqrt((lon-x0)**2)
    i0 = np.where(dx == min(dx))[0][0]
    dy = np.sqrt((lat-y0)**2)
    j0 = np.where(dy == min(dy))[0][0]
    chck = np.sqrt((x0-lon[i0])**2 + (y0-lat[j0])**2)
    if chck > dlat:
      print(f"ERR: WOA lon=lon[i0] lat=lat[j0]")
      raise Exception(f'ERR: no closest WOA grid pnt for x={x0} y={y0} kk={kk}')

    # check if this point already processed:
    if len(Iwoa) > 0:
      dI = int(np.min(np.sqrt((Iwoa-i0)**2 + (Jwoa-j0)**2)))
      if dI == 0: continue

    Iwoa.append(i0)
    Jwoa.append(j0)
    irc = len(Iwoa)
    tt  = T3d[:,j0,i0].squeeze()
    ss  = S3d[:,j0,i0].squeeze()
    if irc == 1:
      Tprf[:,0] = tt
      Sprf[:,0] = ss
    else:
      tt   = np.expand_dims(tt, axis=1)
      Tprf = np.append(Tprf, tt, axis=1)
      ss   = np.expand_dims(ss, axis=1)
      Sprf = np.append(Sprf, ss, axis=1)
    
  if conv2pot:
    pr_ref = 0.
    lat_obs = lat[Jwoa]
    for klr in range(kdm):
      t1d = Tprf[klr,:]
      s1d = Sprf[klr,:]
      inn = len(np.where(~np.isnan(t1d))[0])
      if inn == 0: break
      z0  = ZZ[klr]
      pr_db, pr_pa = msw.sw_press(z0, lat_obs)         
      Tpot = msw.sw_ptmp(s1d, t1d, pr_db, pr_ref)
      Tprf[klr,:] = Tpot

  return Tprf, Sprf, ZZ

def insitu2pot_3D(T3d, S3d, ZZ, LAT, z_ref=0, uref='m'):
  """
    Convert in situ T to potential with pressure 
      reference: z_ref either in m (depth) 
      or dbar (pressure)
    T, S, Z - 3D arrays, lat0 - local latitude (2D array)
    ZZ - can be 1D or 3D array
  """
  import mod_swstate as msw
  import mod_misc1 as mmisc

  kdm, jdm, idm = T3d.shape
  if len(ZZ.shape) == 1:
    Z3d = np.tile(ZZ, idm*jdm).reshape((idm,jdm,kdm))
    Z3d   = np.transpose(Z3d, (2, 1, 0))
  elif len(ZZ.shape) == 3:
    Z3d = ZZ
  else:
    raise Excpetion('ZZ array should be either 1D or 3D')

  Zref = np.zeros((jdm,idm)) 
  if uref == 'm':
    if abs(z_ref) < 1.e-3:
      prref_db = np.zeros((jdm,idm))
      prref_pa = np.zeros((jdm,idm))
    else:
      prref_db, prref_pa = msw.sw_press(Zref, LAT)
  else:
    prref_db = np.zeros((jdm,idm))
    prref_pa = np.zeros((jdm,idm))
 
  Tpot = np.zeros((kdm,jdm,idm))
  for klr in range(kdm):
    print(f'Converting T in situ --> T pot, layer {klr}')
    temp = T3d[klr,:].squeeze()
    sal  = S3d[klr,:].squeeze()
    z0   = Z3d[klr,:].squeeze()
    pr_db, pr_pa = msw.sw_press(z0, LAT)
    tp = msw.sw_ptmp(sal, temp, pr_db, prref_db)
    Tpot[klr,:] = tp

  return Tpot

def insitu2pot_1D(T1d, S1d, Z1d, lat0, z_ref=0, uref='m', printT=False):
  """
    Convert in situ T to potential with pressure reference: z_ref either in m (depth) 
    or dbar (pressure)
    T, S, Z - 1D arrays, 1 profile, lat0 - local latitude
  """
  import mod_swstate as msw
  import mod_misc1 as mmisc

  if uref == 'm':
    prref_db, prref_pa = msw.sw_press(z_ref, lat0)
  else:
    prref_db = z_ref

  Tpot = np.zeros((len(T1d)))
  for klr in range(len(T1d)):
    temp = T1d[klr]
    sal  = S1d[klr]
    z0   = Z1d[klr]
    if np.isnan(temp): continue
    pr_db, pr_pa = msw.sw_press(z0, lat0)
    Tpot[klr] = msw.sw_ptmp(sal, temp, pr_db, prref_db)

  if printT:
    mmisc.print_3col(Z1d,T1d,Tpot)

  return Tpot
 
def pot2insitu_1D(T1d, S1d, Z1d, lat0, z_ref=0, uref='m', printT=False):
  """
    Convert potential T to insitu with pressure reference for T potential: 
    z_ref either in m (depth) 
    or dbar (pressure)
    Z1d - depths (m) for computing in situ T
    T, S, Z - 1D arrays, 1 profile, lat0 - local latitude
  """
  import mod_swstate as msw
  import mod_misc1 as mmisc

  if uref == 'm':
    pr_db, _ = msw.sw_press(z_ref, lat0)
  else:
    pr_db = z_ref

  Tsitu = np.zeros((len(T1d)))
  for klr in range(len(T1d)):
    temp = T1d[klr]
    sal  = S1d[klr]
    z0   = Z1d[klr]
    if np.isnan(temp): continue
    prref_db, _ = msw.sw_press(z0, lat0)
    Tsitu[klr] = msw.sw_ptmp(sal, temp, pr_db, prref_db)

  if printT:
    mmisc.print_3col(Z1d,T1d,Tsitu)

  return Tsitu










