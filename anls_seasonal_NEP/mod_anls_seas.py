"""
  modules/ functions for analysis of seasonal runs
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import xarray
#import yaml
from yaml import safe_load
import importlib

PPTHN = []
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mmom6
import mod_misc1 as mmisc

def xsct_segments_woa(sctnm, lonW, latW, fyaml='paths_seasfcst.yaml', \
                      fyaml_paths='pypaths_gfdlpub.yaml'):
  """
    Find end points of the xsection segments
    similar to segments defined for NEP MOM6
    WOA grid is Mercator
  """
  with open(fyaml) as ff:
    pthseas = safe_load(ff)

  Is  = pthseas['ANLS_NEP'][sctnm]['II']
  Js  = pthseas['ANLS_NEP'][sctnm]['JJ']

  nlegs   = len(Is) - 1
  IJ      = np.zeros((nlegs+1,2))
  IJ[:,0] = Is
  IJ[:,1] = Js

  lonW = np.where(lonW < 0, lonW+360., lonW)
# Hgrid lon. lat:
  with open(fyaml_paths) as ff:
    gridfls = safe_load(ff)

  pthtopo    = gridfls['MOM6_NEP']['seasonal_fcst']['pthgrid']
  fgrid      = gridfls['MOM6_NEP']['seasonal_fcst']['fgrid']
  ftopo_mom  = gridfls["MOM6_NEP"]["seasonal_fcst"]["ftopo"]
  hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
  hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
  dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
  dfgrid_mom = os.path.join(pthtopo, fgrid)
  hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
  HH         = dstopo_nep['depth'].data
  HH         = np.where(HH < 1.e-20, np.nan, HH)
  HH         = -HH
  HH         = np.where(np.isnan(HH), 1., HH)
  DX, DY     = mmom6.dx_dy(hlon,hlat)
  SGMT       = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive',\
                                     check_pole=False)
  II         = SGMT.I_indx
  JJ         = SGMT.J_indx
  nLeg       = SGMT.Leg_number
  hLsgm1     = SGMT.half_Lsgm1
  hLsgm2     = SGMT.half_Lsgm2
  XX         = hlon[JJ,II]
  YY         = hlat[JJ,II]
  Hbtm       = HH[JJ,II]
  LSgm       = np.zeros((len(II)))  # total segment length = half1 + half2

  for ik in range(len(II)):
     LSgm[ik] = hLsgm1[ik] + hLsgm2[ik]

# Distance along the section
# normalize by the total distance
  Lsection = mmisc.dist_sphcrd(YY[-1],XX[-1],YY[0],XX[0]) # total length section, m
  Xdist = np.cumsum(LSgm)
  Xdist = Xdist-Xdist[0]
  Xdist = Xdist/Xdist[-1]*Lsection*1.e-3  # normalized, km

# Find end points of the segments for WOA:
# For latitudinal sections - simply follow the lat:
  nlat = sctnm.split("_")[1]
  if nlat[-1] == 'N':
    lat1 = float(nlat[:2])
    lat2 = lat1
    lon1 = XX[1]
    lon2 = XX[-1]
    Xs = np.array([lon1, lon2])
    Ys = np.array([lat1, lat2])
  else: 
# Not a latitud. section:
    Xs = np.zeros((len(Is)))
    Ys = np.zeros((len(Is)))
    for ii in range(len(Is)):
      i0 = Is[ii]
      j0 = Js[ii]
      Xs[ii] = hlon[j0,i0]
      Ys[ii] = hlat[j0,i0]    

# Find WOA indices:
  IsWOA = np.zeros((len(Xs)), dtype='int') - 999
  JsWOA = np.zeros((len(Xs)), dtype='int') - 999
  for ii in range(len(Xs)):
    x0 = Xs[ii]
    y0 = Ys[ii]
    dlat = abs(latW - y0)
    ilat = np.argmin(dlat)
    lat0 = latW[ilat]
    LAT  = np.zeros((len(lonW))) + lat0
    DD = mmisc.dist_sphcrd(LAT, lonW, y0, x0)
    ilon = np.argmin(DD)

    IsWOA[ii] = ilon
    JsWOA[ii] = ilat

# MOM6 info for plotting:
  indxsct   = np.arange(len(Xdist))
  dim_name  = "sect_indx"
  darr_btm  = xarray.DataArray(Hbtm, dims=(dim_name), \
            coords={dim_name: indxsct})
  darr_dist = xarray.DataArray(Xdist, dims=(dim_name), \
            coords={dim_name: indxsct})
  darr_lon  = xarray.DataArray(XX, dims=(dim_name), \
            coords={dim_name: indxsct})
  darr_lat  = xarray.DataArray(YY, dims=(dim_name), \
            coords={dim_name: indxsct})
  dsetBtm   = xarray.Dataset({
              "Hbtm_section": darr_btm, \
              "Dist_section": darr_dist, \
              "Lon_section": darr_lon, \
              "Lat_section": darr_lat
              }) 

  return IsWOA, JsWOA, dsetBtm

def season_decade_woa(YR,MM):
  """
    Find WOA season and year span for decadal averages
    File naming convention:
    woa23_[DECA]_[v][tp][ft][gr].[form_end] - all formats, except NetCDF
    woa23_[DECA]_[v][tp]_[gr].[form_end] - NetCDF format
    where:
    [DECA] - decade
    [v] - variable
    [tp] - time period
    [ft] - field type
    [gr] - grid
    [form_end] - file name extention
    Note: '.dat' - ASCII; '.csv' - comma separated value; 
          '.dbf', '.shp', '.shx' - ArcGIS shape files; '.nc' - netCDF files
  """
  SEAS = {"1" : 13,
          "2" : 13,
          "3" : 13, 
          "4" : 14,
          "5" : 14, 
          "6" : 14,
          "7" : 15, 
          "8" : 15,
          "9" : 15, 
          "10": 16,
          "11": 16, 
          "12": 16}

  DECA = np.array([[1955, 1964],
                   [1965, 1974],
                   [1975, 1984],
                   [1985, 1994],
                   [1995, 2004],
                   [2005, 2014],
                   [2015, 2022]])


  if MM > 12:
    seas = 0    # annual
  else:
    seas = SEAS[f"{MM}"]

  if YR > np.max(DECA) or YR < np.min(DECA):
    raise Exception(f"{YR} is not in the time range for WOA")
  ii = np.where((DECA[:,0]<=YR) & (DECA[:,1]>=YR))[0][0]
  yr1, yr2 = DECA[ii,:]
  yr1_end = f"{yr1}"[2:5]
  yr2_end = f"{yr2}"[2:5]
  if YR < 1995:
    decade = yr1_end+yr2_end
  elif YR >= 1995 and YR < 2005:
    decade = "95A4"
  elif YR >= 2005 and YR < 2015:
    decade = "A5B4"
  elif YR >= 2015 and YR < 2023:
    decade = "B5C2"

  return seas, decade, yr1, yr2


def calc_N2(T, S, ZZ, latW):
  """ 
  Input T, S - 2D or 3D arrays
#  Calculate N2 = -g/rho*d2(rho)/dz2 
#  Using in situ T and S - calculate rho relative to the
#  mid-grid depth - following Chelton, 1996
# Use the neutral density gradient method, Chelton et al., 1996

  Note for model output, T, S are in the mid-depths,
  then rho's are in the interfaces, surface rho is missing - BC

  """
  import mod_swstate
  #importlib.reload(mod_swstate)
  from mod_swstate import sw_press
  from mod_swstate import adiab_Tgrad
  from mod_swstate import sw_ptmp
  from mod_swstate import sw_dens0
  from mod_swstate import sw_smow
  from mod_swstate import sw_seck
  from mod_swstate import sw_dens

#  print('Calculating pressure at midpoints ...')
  kdm, idm = T.shape    
  grav   = 9.81
  rho0   = 1025.0
  Z_phi  = np.zeros(kdm-1)  # depths for e/function Phi
  N2     = np.zeros((kdm-1, idm))
  n2fill = 1.e-8     # missing values

  for kk in range(kdm-1):
#    print(' Layer {0}'.format(kk))
    z1 = ZZ[kk]
# To avoid negative P at the surface:
    if abs(z1) < 1.e-30: z1 = 0.
    z2 = ZZ[kk+1]
    t1 = T[kk,:].squeeze()
    t2 = T[kk+1,:].squeeze()
    s1 = S[kk,:].squeeze()
    s2 = S[kk+1,:].squeeze()
    Z1 = z1*np.ones((idm))
    Z2 = z2*np.ones((idm))
    p1_db, p1_pa = sw_press(Z1,latW)  # pressure upper interface
    p2_db, p2_pa = sw_press(Z2,latW)  # pressure bottom interface
  #
  # Find mid-point of the layers
    p0_db = 0.5*(p1_db+p2_db)
    z0    = 0.5*(z1+z2)
    Z_phi[kk] = z0
#    print(f'z1={z1:5.1f} z2={z2:5.2f} p1={np.nanmin(p1_db):5.1f}'  +\
#           f' p2={np.nanmin(p2_db):5.1f} z0={z0:5.1f}')
  #
  # Calculate rho(z1--->z0) with depth reference at midpoint
  # and rho(z2--->z0) at midpoint 
    t_z1z0 = sw_ptmp(s1,t1,p1_db,p0_db)
    t_z2z0 = sw_ptmp(s2,t2,p2_db,p0_db)
    rho_z1z0 = sw_dens(s1,t_z1z0,p0_db)
    rho_z2z0 = sw_dens(s2,t_z2z0,p0_db)

  # Calculate d(rho)/dz for z0 - center-difference
    drho_dz  = (rho_z1z0 - rho_z2z0)/(z1 - z2)
    N2z0     = -grav/rho0*drho_dz
    N2[kk,:] = N2z0
#    print(f'z1={z1} z2={z2} z0={z0}')
#    print(f'k={kk}, z={z0}, min/max N2: {np.nanmin(N2z0)}, {np.nanmax(N2z0)}')
#
# If N2 < 0 - density inversion happens in some ~homogeneous layers
# when parcels are brought down from z1 and up from z2
# replace with above N2 or below of surface layer
  print('Fixing N2<0')
  for kk in range(kdm-1):
    N2z0=N2[kk,:].squeeze()
    if kk > 0:
      dmm = N2[kk-1,:].squeeze()
      N2z0 = np.where(N2z0 < 0., dmm, N2z0)
    else:
      dmm = N2[kk+1,:].squeeze()
      N2z0 = np.where(N2z0 < 0., dmm, N2z0)

    N2z0     = np.where(N2z0<0, n2fill, N2z0)
    N2[kk,:] = N2z0

#    print(f'k={kk}, z={z0}, min/max N2: {np.nanmin(N2z0)}, {np.nanmax(N2z0)}')

  return N2, Z_phi

def solve_SturmLiouville(N2z, Hb0, Z_phi, latj, mode=1):
  """

  # Numerically Solve Sturm-Liouville e/value problem 

  """
  import mod_solver as msolv
  importlib.reload(msolv)
 
  omg = 7.29e-5 # Earth angular velocity
  #tic = timeit.default_timer()
  #ticR = timeit.default_timer()
#  print(f'Solving Strum-Liouville, Hbtm={Hb0:4.1f}, requested mode={mode}')

  if np.isnan(N2z[0]): 
    print(f'N2 profile is nan, land point? Skipping')
    return np.nan, np.nan, np.nan 

  k1 = np.where(np.isnan(N2z))[0]
  if k1.size:   
    kbtm = k1[0]-1
  else:
# Bottom deeper than the last layer
# Check if this is so
    kbtm = N2z.shape[0] - 1
    if abs(Hb0) < abs(Z_phi[-1]):
      raise Exception(f"Bottom={Hb0} last z={Z_phi[-1]} expected Z_phi[-1]>Hb0")

  zbtm = Hb0
  # Create Matrix A with Dk, Dk+1 for 2nd derivative of N2
  AA = msolv.form_mtrxA(Z_phi, kbtm, zbtm)

  # Form 1/N2*AA:
  # The matrix is cutoff by the bottom depth
  ka, na = AA.shape
  N2A = np.zeros((ka,na))
  for kk in range(ka):
  #  n2 = N2zF[kk]
    n2 = N2z[kk]  
    if n2 == 0:
      n2=1.e-12
    N2A[kk,:] = 1./n2*AA[kk,:]

# For now: use python eig function, for this
# form Matrix A unfolding
# the 3 -elemnts form and find eigenvalues/vectors
# W - eigenvalues, not sorted out by magnitude!
# V- corresponding eigenvectors (in columns)
    W, V = msolv.eig_unfoldA(kbtm, N2A)

# Calculate
# Choose lmbd1 = as min(abs(W)) and lmbd1<0
  absW = np.abs(W)
  Indx = np.argsort(absW)
  im   = Indx[mode-1]    # 1st mode
#  im2  = Indx[1]    # 2nd mode
#  im   = np.argmin(absW)
  if W[im] > 0.:
    print('ERROR: W[im] >0: im={0}, W[im]={1}, zbtm={4:6.1f}, ii={2}, jj={3}'.\
           format(im, W[im], ii, jj, zbtm))

#  latj = latW[ik]
  Rearth = 6371.e3  # Earth R
  fcor = 2.*omg*np.sin(np.deg2rad(latj))
  betta = (2.*omg*np.cos(np.deg2rad(latj)))/Rearth

  if abs(latj) >= 5.0:
    lmbd = np.sqrt(-1./(W[im]*fcor**2))
  else:
    lmbd = (-1./(4.*(betta**2)*(W[im])))**(1./4.)

  RsbNum = lmbd*1.e-3  # km
  Phi    = V[:,im]   # eigenfunction
  mu     = W[im]     # e/value ( <0)
  Cphs   = np.sqrt(-1./mu)  # m-th mode phase speed m/s

# Insert BCs' for Phi:
# Surface:
  Phi = np.insert(Phi, 0, 0.)
# Add phi at the bottom
  Phi = np.append(Phi, 0.)

  return Phi, RsbNum, Cphs 

def derive_dP_WOA(ZZ, A3d):
  """
  Estimate layer thicknesses for WOA Z layers
  Use any 3D field to find bottom
  """
  import mod_mom6 as mmom6
  kdm, jdm, idm = A3d.shape
  dP = A3d.copy()*0.
  ZM = mmom6.zz2zm(ZZ)
  ZM = np.append(ZM, ZZ[-1]) 

  for kk in range(kdm):
    a2d = A3d[kk,:].squeeze()
    if kk == 0: 
      dP[kk,:] = abs(ZM[0])
# Keep land mask to avoid land points
      J,I = np.where(np.isnan(a2d))
      dP[kk,J,I] = np.nan
    else:
      dP[kk,:] = abs(ZM[kk]-ZM[kk-1])
      J,I = np.where(np.isnan(a2d))
      dP[kk,J,I] = 0.
      dP[kk,:] = np.where(np.isnan(dP[kk-1,:]), np.nan, dP[kk,:])  # keep land mask

  return dP

def derive_dP_from_ZZtopo(ZZ, HH):
  """
  Derive layer thickness array from interface depths ZZ and bottom topography HH
  dP below bottom = 0.
  ZZ can be 1D or 3D, HH = 2D
  ZZ and HH use negative depths
  """
  if np.min(HH) > 0.:
    raise Exception("HH topography depths has to be <0")
  if np.min(ZZ) > 0:
    print(f"WARNING: ZZ > 0, converting to negative depths")
    ZZ = -ZZ

  ndim = len(ZZ.shape)
  jdm, idm = HH.shape
  kdm = ZZ.shape[0]-1
  if ndim == 1:
    Z3d   = np.tile(ZZ, idm*jdm).reshape((idm,jdm,kdm+1))
    Z3d   = np.transpose(Z3d, (2, 1, 0))
  elif ndim == 3:
    Z3d = ZZ

  dP3d = abs(np.diff(Z3d, axis=0))

# Make 0-thicknesses
  for kk in range(kdm):
    zup  = Z3d[kk,:].squeeze()
    zbtm = Z3d[kk+1,:].squeeze()
    dmm  = dP3d[kk,:].squeeze()
    dmm  = np.where(zup < HH, 0., dmm)
    dmm  = np.where( (zup > HH) & (zbtm < HH), abs(HH-zup), dmm)
    dP3d[kk,:] = dmm

  return dP3d

def find_closest_output(pthoutp, dnmb0, fld='oceanm'):
  """
    Find closest output file to given date
    MOM/HYCOM file naming assumed: fld_YYYY_DAY.nc
  """
  import os
  import mod_time as mtime

  if not os.path.isdir(pthoutp):
    print(f'not exist: {pthoutp}')
    return 0, 0, 0, 'none'

# Find 1st ouptut date:
  LF = os.listdir(pthoutp)
#  print(LF)
  if len(LF) == 0:
    raise Exception(f'No files in {pthoutp}')
#  YR0 = np.zeros((len(LF)))
#  JD0 = np.zeros((len(LF)))
  icc = -1
  for fls in LF:
    bsname = fls.split(".")[0]
    fld_name =  bsname.split("_")[0]
    if not fld_name == fld: continue
    yrf    = int(bsname.split("_")[1])
    jday   = int(bsname.split("_")[2])
    icc += 1
    if icc == 0:
      YR0 = np.array([yrf])
      JD0 = np.array([jday])
    else:
      YR0 = np.append(YR0,yrf)
      JD0 = np.append(JD0,jday)
#    print(f'icc={icc} yr={yrf} jday={jday}')

  if icc < 0:
    print(f'ERR: no files found {fld}* in {pthoutp}')
    print(LF)
    raise Exception(f'Check fld={fld} ???')
 
  DNMB  = mtime.jday2dnmb(YR0, JD0)
  D     = np.sqrt((DNMB-dnmb0)**2)
  ii    = np.argmin(D)
  year  = YR0[ii]
  jday  = JD0[ii]
#  print(f'Min ii={ii} year={year} jday={jday} dnmb={mtime.jday2dnmb(1994,58)}')
  dnmb  = mtime.jday2dnmb(year, jday)
  flname= LF[ii]

  return int(year), int(jday), dnmb, flname

def derive_ice_contour(AA, tz0=0.15, nmin=10):
  """
    Deduce coordinates of the ice edge contour, AA 2D field with 0 - 1 ice conc
    nmin - min # of points to keep the contour
  """
  plt.ioff()
  figA = plt.figure(10,figsize=(8,8))
  plt.clf()

  axA1 = plt.axes([0.1, 0.2, 0.7, 0.7])
  CS   = axA1.contour(AA,[tz0])

  axA1.axis('equal')
#  axA1.set_xlim([xl1,xl2])
#  axA1.set_ylim([yl1,yl2])

  SGS  = CS.allsegs[0]  # should be only 1 contoured value
  nsgs = len(SGS)
  
# Delete all closed contours and small contours
  CNTR = []
  for isg in range(nsgs):
    XY = SGS[isg]
    X  = XY[:,0]
    Y  = XY[:,1]

    if len(X) <= nmin: continue

    dEnd = np.sqrt((X[0]-X[-1])**2+(Y[0]-Y[-1])**2)
    if dEnd < 1.:
      continue

    CNTR.append(XY)
  
  plt.close(figA)

  plt.ion()

  return CNTR


