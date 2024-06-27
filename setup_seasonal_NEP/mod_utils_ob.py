"""
  Functions for extracting and writing
  OB file for regional MOM6
"""
import xarray
import os
import importlib
import numpy as np
import sys
import matplotlib.pyplot as plt
#import pickle
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
sys.path.append('./seasonal-workflow')
 
import mod_mom6 as mom6util
importlib.reload(mom6util)

def load_var(spear_dir, varnm, yr1=1993, yr2=2020, mstart=1, fzint=False):
  """
    Extract SPEAR monthly climatology
    fzint = True - include z layer interface depths into returned dataset
  """
  flname = f'spear_moclim.mstart{mstart:02d}.{yr1}-{yr2}.{varnm}.nc'
  dflinp = os.path.join(spear_dir, flname)

  if not varnm == 'ssh':
    spear_var = xarray.open_dataset(dflinp).rename({'z_l': 'z'})
  else:
    spear_var = xarray.open_dataset(dflinp)
    fzint = False

# Make sure missing/fill value was properly used
  assert np.abs(spear_var).max() < 1e5

  if fzint:
    zzi = spear_var.z_i
    dset_out = xarray.Dataset({f"{varnm}": spear_var[varnm], "z_interf": zzi})
  else:
    dset_out = xarray.Dataset({f"{varnm}":spear_var[varnm]})

  return dset_out



def add_valid_time(ds):
  months = ds.lead.values.astype(int) + mon
  dates = [dt.datetime(1993, m, 1) if m <= 12 else dt.datetime(1994, m-12, 1) for m in months]
  ds['time'] = (('lead', ), dates)
  ds = ds.swap_dims({'lead': 'time'})
  return ds


def add_coords(da):
  xcoord = [str(d) for d in da.coords if d in ['xh', 'xq']]
  ycoord = [str(d) for d in da.coords if d in ['yh', 'yq']]
  assert len(xcoord) == len(ycoord) == 1
  xc = xcoord[0]
  yc = ycoord[0]
  coord_tuple = (yc, xc)
  coord_sel = {xc: da[xc], yc: da[yc]}
  lat = next(static[v] for v in ['geolat', 'geolat_u', 'geolat_v'] if xc in static[v].coords and yc in static[v].coords)
  lon = next(static[v] for v in ['geolon', 'geolon_u', 'geolon_v'] if xc in static[v].coords and yc in static[v].coords)
  da['lat'] = (coord_tuple, lat.sel(**coord_sel).data)
  da['lon'] = (coord_tuple, lon.sel(**coord_sel).data)
  return da

def spear2momOB_gmapi(segments, lon_spear, lat_spear, dflout, fsave=True):
  """
   For interpolation find mapping indices: 4 grid points on SPEAR grid around
   a MOM grid point along the OB segment
   SPEAR lon/ lat - 2D array of u, v, or h-points
  """
  import mod_regmom as mrmom

  dset = xarray.Dataset()
  nsgm = len(segments)
  for ii in range(nsgm):
    isgm = ii+1
    dset_segm = segments[ii]
    xOB   = dset_segm.coords.lon.data
    yOB   = dset_segm.coords.lat.data
    npnts = len(xOB)

    INDX = np.zeros((npnts,4)).astype(int)-999
    JNDX = np.zeros((npnts,4)).astype(int)-999
    for ikk in range(npnts):
      x0 = xOB[ikk]
      y0 = yOB[ikk]
      ixx, jxx = mrmom.find_gridpnts_box(x0, y0, lon_spear, lat_spear, dhstep=0.8) 
      INDX[ikk,:] = ixx
      JNDX[ikk,:] = jxx

    darri = xarray.DataArray(INDX, dims=(f"npnts_{isgm:03d}","vertices"), \
                  coords={f"npnts_{isgm:03d}": np.arange(npnts), \
                  "vertices": np.arange(1,5)})
    darrj = xarray.DataArray(JNDX, dims=(f"npnts_{isgm:03d}","vertices"), \
                  coords={f"npnts_{isgm:03d}": np.arange(npnts), \
                  "vertices": np.arange(1,5)})

    if isgm == 1:
      dset = xarray.Dataset({f"indx_segm{isgm:03d}": darri, f"jndx_segm{isgm:03d}": darrj})
    else:
      dstmp = xarray.Dataset({f"indx_segm{isgm:03d}": darri, \
                  f"jndx_segm{isgm:03d}": darrj})
      dset = xarray.merge([dset, dstmp])
     
    f_plt = False
    if f_plt:
      plt.ion()

      fig1 = plt.figure(1,figsize=(9,8))
      plt.clf()
      ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])

      if min(xOB) > 0. and min(lon_spear) < 0.:
        lon_spear = np.where(lon_spear<0, lon_spear+360, lon_spear)

      ax1.plot(xOB, yOB,'r*')
      for ikk in range(npnts):
#      ax1.plot(x0,y0,'r*')
        ixx = INDX[ikk,:]
        jxx = JNDX[ikk,:]
        xv  = np.append(lon_spear[ixx], lon_spear[ixx[0]])
        yv  = np.append(lat_spear[jxx], lat_spear[jxx[0]])
        ax1.plot(xv, yv, 'o-')

      ax1.axis('scaled')
      ax1.set_title('mod_utils_ob.py: MOM6 NEP OB segment 1, SPEAR gmapi')

  if fsave:
    print(f"Saving gmapi indices --> {dflout}")
    dset.to_netcdf(dflout, format='NETCDF3_64BIT', engine='netcdf4')

  return

def derive_obsegm_lonlat(hgrid, segments, isgm):
  """
    Derive lon lat and construct data sets
  """
  nsgm   = isgm+1 
  dset_segm = segments[isgm]
  xOB    = dset_segm.coords.lon.data
  yOB    = dset_segm.coords.lat.data
  npnts  = len(xOB)
  nx     = dset_segm.nx
  ny     = dset_segm.ny
  bndary = dset_segm.border 

  if bndary == 'west':
    lon_segm = hgrid['x'].isel(nxp=0).data
    lat_segm = hgrid['y'].isel(nxp=0).data
  elif bndary == 'north':
    lon_segm = hgrid['x'].isel(nyp=-1).data
    lat_segm = hgrid['y'].isel(nyp=-1).data
  elif bndary == 'east':
    lon_segm = hgrid['x'].isel(nxp=-1).data
    lat_segm = hgrid['y'].isel(nxp=-1).data
  elif bndary == 'south':
    lon_segm = hgrid['x'].isel(nyp=0).data
    lat_segm = hgrid['y'].isel(nyp=0).data

  nxsgm = [x for x in range(0,nx)]
  nysgm = [x for x in range(0,ny)]
  xdim  = f'nx_segment_{nsgm:03d}'
  ydim  = f'ny_segment_{nsgm:03d}'
  if ny == 1:
    darr_lon = xarray.DataArray(lon_segm, dims=(xdim),
                              coords={xdim: nxsgm})
    darr_lat = xarray.DataArray(lat_segm, dims=(xdim),
                              coords={xdim: nxsgm})
  else:
    darr_lon = xarray.DataArray(lon_segm, dims=(ydim),
                              coords={ydim: nysgm})
    darr_lat = xarray.DataArray(lat_segm, dims=(ydim),
                              coords={ydim: nysgm})

  dset_segm = xarray.Dataset({f"lon_segment_{nsgm:03d}": darr_lon, \
           f"lat_segment_{nsgm:03d}": darr_lat})

  return dset_segm

 

def derive_obsegm_3D(hgrid, ds, segments, isgm, varnm, INDX, JNDX, icegrid=[], rot_angle=[]):
  """
    Derive data set for an OB segment = isgm for variable = varnm
    Scalar (h-point) or u/v 3D fields (+ time dim)  --> 2D vertical sections 
    ds - xarray of climatology field for 1 variable, subset to the region
    hgrid - topo/grid for the climatology field (MOM6 NEP grid)
    segments - xarray with information for OB segments created in boundary.Segment
    isgm - segment # [0, 1, ...]

    INDX, JNDX - n x 4 arrays of grid points (gmapi) for bilinear interpolation

    icegrid - spear static file with grid info, need rotation angles for
              U-vectors rotation onto MOM NEP grid from SPEAR
              alternatively, rotation angle rot_angle is provided
              angle is from J-axis to true North dir (math convention + is c/clockwise)

 Use bilinear interpolation 
 Points are mapped onto a reference square for interpolation then 
 Perform interpolation on reference rectangle, that is 
 similar to interp on actual rectangle
 The advantage of using a reference quadrilateral is that basis
 functions computed once and used for all grid points - saved lots of 
 computational time
     See p. 83-90, Gockenbach, "Understanding Finite-Element Methods" textbook

  """
  import mod_bilinear as mblnr
  importlib.reload(mblnr)

# Find basis functions for a reference rectangle:
  phi1,phi2,phi3,phi4 = mblnr.basisFn_RectRef()
  phi_basis           = np.array([phi1, phi2, phi3, phi4]).transpose() # basis funs in columns

  nsgm  = isgm+1 
  dset_segm = segments[isgm]
  xOB   = dset_segm.coords.lon.data
  yOB   = dset_segm.coords.lat.data
  npnts = len(xOB)
  nx    = dset_segm.nx
  ny    = dset_segm.ny

# Make sure that gmapi is for the right section:
  assert npnts==len(INDX), "INDX and OB segments mismatch in length"

# SPEAR grid coords:
  if varnm == 'uo':
    xSP   = ds.xq.data
    ySP   = ds.yh.data
  elif varnm == 'vo':
    xSP   = ds.xh.data
    ySP   = ds.yq.data
  else:
    xSP   = ds.xh.data
    ySP   = ds.yh.data
 
  dtm   = ds.time.data
  AA    = ds[f'{varnm}'].data
  ntime = len(dtm) 

# interface depths and dz
  zzi  = ds.z_interf.data
  dz   = np.diff(zzi)
  if ny == 1:
    dz2d = np.tile(dz,(nx,1)).transpose() 
  else:
    dz2d = np.tile(dz,(ny,1)).transpose() 

  dz3d = np.tile(dz2d,(ntime,1,1))

# Assuming 3D arrays and time dimenstion:
  ndim = len(AA.shape)
  if ndim < 4:
    raise Exception(f'4D array (time, depth, y, x) is assumed, check array dim')

# No vertical interpolation, so use n levels from input data:
  nzlev = AA.shape[1]
  Ai = np.zeros((ntime,nzlev,ny,nx))+1.e30

  for itime in range(ntime):
    print(f'Time = {itime+1}')
    A   = AA[itime,:]
# Fill missing values (bottom/land):
    dmm =  mom6util.fill_land3d(A)
    assert np.max(abs(dmm)) < 1.e20, "Bottom/land not filled correctly"
    assert not np.isnan(np.max(abs(dmm))), "Bottom/land not filled correctly"
    A3d = dmm.copy()

# Horizontal interpolation - bilinear
    for ikk in range(npnts):
      x0 = xOB[ikk]
      y0 = yOB[ikk]
      II = np.squeeze(INDX[ikk,:])
      JJ = np.squeeze(JNDX[ikk,:])
      II, JJ = mblnr.sort_gridcell_indx(II,JJ)
      ii1, jj1 = II[0], JJ[0]
      ii2, jj2 = II[1], JJ[1]
      ii3, jj3 = II[2], JJ[2]
      ii4, jj4 = II[3], JJ[3]

      xx = xSP[II]
      yy = ySP[JJ]

# Use cartesian coordinates for mapping
      if x0 < 0.: x0 = x0+360.
      xx = np.where(xx<0., xx+360., xx)
      XV, YV, x0c, y0c = mblnr.lonlat2xy_wrtX0(xx, yy, x0, y0)
      xht, yht = mblnr.map_x2xhat(XV, YV, x0c, y0c)   # cartesian coord

      aa1 = np.squeeze(A3d[:,jj1,ii1])
      aa2 = np.squeeze(A3d[:,jj2,ii2])
      aa3 = np.squeeze(A3d[:,jj3,ii3])
      aa4 = np.squeeze(A3d[:,jj4,ii4])

      HT  = np.array([aa1, aa2, aa3, aa4]).transpose()
      hintp  = mblnr.bilin_interp1D(phi1, phi2, phi3, phi4, xht, yht, HT)

# Check: abs. values of interpolated values <= original data
      mxHT = np.max(abs(HT), axis=1)
      dmx  = abs(hintp)/mxHT
      if max(dmx-1.) > 0.1:
        print(f"!!! segm{nsgm} {varnm} Min/Max test violated: ikk={ikk} dlt: {np.max(dmx)}")
        if max(dmx-1.) > 1.:
          raise Exception("MinMax test: Interp error Check interpolation")

      if ny == 1:
        Ai[itime, :, 0, ikk] = hintp
      else:
        Ai[itime, :, ikk, 0] = hintp
    
# Construct data set:
  if ny == 1:
    dz3d = np.expand_dims(dz3d, axis=2)
  else:
    dz3d = np.expand_dims(dz3d, axis=3)
  nxsgm = [x for x in range(0,nx)]
  nysgm = [x for x in range(0,ny)]
  stime = [x for x in range(1,ntime+1)]
  nzsgm = [x for x in range(0,nzlev)]
  sgmnm = f"segment_{nsgm:03d}"
  dim1  = "time"
  dim2  = f"nz_{sgmnm}"
  dim3  = f"ny_{sgmnm}"
  dim4  = f"nx_{sgmnm}"
  darr_var = xarray.DataArray(Ai, dims=(dim1, dim2, dim3, dim4),
             coords={dim1: stime, dim2: nzsgm, dim3: nysgm, dim4: nxsgm})
  darr_dz  = xarray.DataArray(dz3d, dims=(dim1, dim2, dim3, dim4),
             coords={dim1: stime, dim2: nzsgm, dim3: nysgm, dim4: nxsgm})
#  f"{varnm}_segment_{nsgm:03d}"
  dset_segm = xarray.Dataset({f"{varnm}_segment_{nsgm:03d}": darr_var, \
           f"dz_{varnm}_segment_{nsgm:03d}": darr_dz})

  return dset_segm
   
 
def derive_obsegm_ssh(hgrid, ds, segments, isgm, INDX, JNDX):
  """
    Derive data set for an OB segment = isgm for variable = varnm
    SSH fields at h-points 2D fields --> 1D OB segment 
    ds - xarray of climatology field for 1 variable, subset to the region
    hgrid - topo/grid for the climatology field
    segments - xarray with information for OB segments created in boundary.Segment
    isgm - segment # [0, 1, ...]

    INDX, JNDX - n x 4 arrays of grid points (gmapi) for bilinear interpolation

 Use bilinear interpolation 

  """
  import mod_bilinear as mblnr
  importlib.reload(mblnr)

# Find basis functions for a reference rectangle:
  phi1,phi2,phi3,phi4 = mblnr.basisFn_RectRef()
  phi_basis           = np.array([phi1, phi2, phi3, phi4]).transpose() # basis funs in columns

  nsgm  = isgm+1 
  dset_segm = segments[isgm]
  xOB   = dset_segm.coords.lon.data
  yOB   = dset_segm.coords.lat.data
  npnts = len(xOB)
  nx    = dset_segm.nx
  ny    = dset_segm.ny

# Make sure that gmapi is for the right section:
  assert npnts==len(INDX), "INDX and OB segments mismatch in length"

# SPEAR grid coords:
  xSP   = ds.xh.data
  ySP   = ds.yh.data
 
  dtm   = ds.time.data
  AA    = ds['ssh'].data
  ntime = len(dtm) 

# 2D array and time dimension:
  ndim = len(AA.shape)
  if ndim != 3:
    raise Exception(f'3D array (time, y, x) is assumed for ssh, check array dim')

  Ai = np.zeros((ntime,ny,nx))+1.e30

  for itime in range(ntime):
    print(f'Time = {itime+1}')
    A   = AA[itime,:]
# Fill missing values (bottom/land):
    dmm =  mom6util.fill_land3d(A)
    assert np.max(abs(dmm)) < 1.e20, "Bottom/land not filled correctly"
    assert not np.isnan(np.max(abs(dmm))), "Bottom/land not filled correctly"
    A2d = dmm.copy()

# Horizontal interpolation - bilinear
    for ikk in range(npnts):
      x0 = xOB[ikk]
      y0 = yOB[ikk]
      II = np.squeeze(INDX[ikk,:])
      JJ = np.squeeze(JNDX[ikk,:])
      II, JJ = mblnr.sort_gridcell_indx(II,JJ)
      ii1, jj1 = II[0], JJ[0]
      ii2, jj2 = II[1], JJ[1]
      ii3, jj3 = II[2], JJ[2]
      ii4, jj4 = II[3], JJ[3]

      xx = xSP[II]
      yy = ySP[JJ]

# Use cartesian coordinates for mapping
      if x0 < 0.: x0 = x0+360.
      xx = np.where(xx<0., xx+360., xx)
      XV, YV, x0c, y0c = mblnr.lonlat2xy_wrtX0(xx, yy, x0, y0)
      xht, yht = mblnr.map_x2xhat(XV, YV, x0c, y0c)   # cartesian coord

      aa1 = np.squeeze(A2d[jj1,ii1])
      aa2 = np.squeeze(A2d[jj2,ii2])
      aa3 = np.squeeze(A2d[jj3,ii3])
      aa4 = np.squeeze(A2d[jj4,ii4])

      HT  = np.array([aa1, aa2, aa3, aa4]).transpose()
      hintp  = mblnr.bilin_interp(phi1, phi2, phi3, phi4, xht, yht, HT)

# Check: abs. values of interpolated values <= original data
      mxHT = np.max(abs(HT))
      dmx  = abs(hintp)/mxHT
      if (dmx-1.) > 0.1:
        print(f"!!! segm{nsgm} ssh Min/Max test violated: ikk={ikk} dlt: {np.max(dmx)}")
        if (dmx-1.) > 1.:
          raise Exception("MinMax test: Interp error Check interpolation")

      if ny == 1:
        Ai[itime, 0, ikk] = hintp
      else:
        Ai[itime, ikk, 0] = hintp
    
# Construct data set:
  nxsgm = [x for x in range(0,nx)]
  nysgm = [x for x in range(0,ny)]
  stime = [x for x in range(1,ntime+1)]
  sgmnm = f"segment_{nsgm:03d}"
  dim1  = "time"
  dim2  = f"ny_{sgmnm}"
  dim3  = f"nx_{sgmnm}"
  darr_var = xarray.DataArray(Ai, dims=(dim1, dim2, dim3),
             coords={dim1: stime, dim2: nysgm, dim3: nxsgm})
  dset_segm = xarray.Dataset({f"zos_segment_{nsgm:03d}": darr_var})

  return dset_segm

def derive_obsegm_uv(hgrid, ds_uo, ds_vo, segments, isgm, theta_rot, varnm, INDX, JNDX):
  """
    Derive data set for an OB segment = isgm for U and V components
    Scalar (h-point) or u/v 3D fields (+ time dim)  --> 2D vertical sections 
    dsu / dsv - xarray of climatology field for U/V variable, subset to the region
    hgrid - topo/grid for the climatology field (MOM6 NEP grid)
    segments - xarray with information for OB segments created in boundary.Segment
    isgm - segment # [0, 1, ...]

    INDXu, JNDXu - n x 4 arrays of grid points (gmapi) for bilinear interpolation
                   for u-points on SPEAR
    similar for v-components INDXv, JNDXv

    icegrid - spear static file with grid info, need rotation angles for
              U-vectors rotation onto MOM NEP grid from SPEAR
              alternatively, rotation angle rot_angle is provided
              angle is from J-axis to true North dir (math convention + is c/clockwise)

   Use bilinear interpolation 

  """
  import mod_bilinear as mblnr
  importlib.reload(mblnr)

# Find basis functions for a reference rectangle:
  phi1,phi2,phi3,phi4 = mblnr.basisFn_RectRef()
  phi_basis           = np.array([phi1, phi2, phi3, phi4]).transpose() # basis funs in columns

  nsgm   = isgm+1 
  dset_segm = segments[isgm]
  xOB    = dset_segm.coords.lon.data
  yOB    = dset_segm.coords.lat.data
  npnts  = len(xOB)
  nx     = dset_segm.nx
  ny     = dset_segm.ny
  bndary = dset_segm.border  

# Make sure that gmapi is for the right section:
  assert npnts==len(INDX), "INDX-U and OB segments mismatch in length"
  assert npnts==len(JNDX), "JNDX-V and OB segments mismatch in length" 

# SPEAR grid coords:
# SPEAR grid coords:
  if varnm == 'u':
    xSP   = ds_uo.xq.data
    ySP   = ds_uo.yh.data
    ds    = ds_uo
  elif varnm == 'v':
    xSP   = ds_vo.xh.data
    ySP   = ds_vo.yq.data
    ds    = ds_vo

  dtm   = ds.time.data
  ntime = len(dtm) 
  U4D   = ds_uo['uo'].data
  V4D   = ds_vo['vo'].data

# interface depths and dz
  zzi  = ds.z_interf.data
  dz   = np.diff(zzi)
  if ny == 1:
    dz2d = np.tile(dz,(nx,1)).transpose() 
  else:
    dz2d = np.tile(dz,(ny,1)).transpose() 

  dz3d = np.tile(dz2d,(ntime,1,1))

# Assuming 3D arrays and time dimenstion:
  ndim = len(U4D.shape)
  if ndim < 4:
    raise Exception(f'4D array (time, depth, y, x) is assumed, check array dim')

# No vertical interpolation, so use n levels from input data:
  nzlev = U4D.shape[1]
  Ai = np.zeros((ntime,nzlev,ny,nx))+1.e30

# For vector rotation, make sure the angle sign matches math convention
# For NEP hgrid saves angle positive clockwise (from axis to true N/E)
# *(-1) make it in math sense: positive c/clockwise axis J/I to true N/E
# also can derive angles at the OB: theta=dset_segm.coords.angle.data
# use plot_rotangleNEP(hgrid, hmask, nfig=1) to plot angles for the whole domain
  if bndary == 'west':
    thetaNEP = -hgrid['angle_dx'].isel(nxp=0).data  
  elif bndary == 'north':
    thetaNEP = -hgrid['angle_dx'].isel(nyp=-1).data  
  elif bndary == 'east':
    thetaNEP = -hgrid['angle_dx'].isel(nxp=-1).data  
  elif bndary == 'south':
    thetaNEP = -hgrid['angle_dx'].isel(nyp=0).data  

  for itime in range(ntime):
    print(f'Time = {itime+1}')
    A   = U4D[itime,:]
# Fill missing values (bottom/land):
    dmm =  mom6util.fill_land3d(A)
    assert np.max(abs(dmm)) < 1.e20, "Bottom/land not filled correctly"
    assert not np.isnan(np.max(abs(dmm))), "Bottom/land not filled correctly"
    U3d = dmm.copy()

    A   = V4D[itime,:]
    dmm =  mom6util.fill_land3d(A)
    assert np.max(abs(dmm)) < 1.e20, "Bottom/land not filled correctly"
    assert not np.isnan(np.max(abs(dmm))), "Bottom/land not filled correctly"
    V3d = dmm.copy()

# Horizontal interpolation - bilinear
    for ikk in range(npnts):
      x0 = xOB[ikk]
      y0 = yOB[ikk]
      II = np.squeeze(INDX[ikk,:])
      JJ = np.squeeze(JNDX[ikk,:])
      II, JJ = mblnr.sort_gridcell_indx(II,JJ)
      ii1, jj1 = II[0], JJ[0]
      ii2, jj2 = II[1], JJ[1]
      ii3, jj3 = II[2], JJ[2]
      ii4, jj4 = II[3], JJ[3]

      xx = xSP[II]
      yy = ySP[JJ]

# Use cartesian coordinates for mapping
      if x0 < 0.: x0 = x0+360.
      xx = np.where(xx<0., xx+360., xx)
      XV, YV, x0c, y0c = mblnr.lonlat2xy_wrtX0(xx, yy, x0, y0)
      xht, yht = mblnr.map_x2xhat(XV, YV, x0c, y0c)   # cartesian coord

# Angle in radians!
      tht1 = theta_rot[jj1,ii1]
      tht2 = theta_rot[jj2,ii2]
      tht3 = theta_rot[jj3,ii3]
      tht4 = theta_rot[jj4,ii4]
      uu1  = np.squeeze(U3d[:,jj1,ii1])
      uu2  = np.squeeze(U3d[:,jj2,ii2])
      uu3  = np.squeeze(U3d[:,jj3,ii3])
      uu4  = np.squeeze(U3d[:,jj4,ii4])
      vv1  = np.squeeze(V3d[:,jj1,ii1])
      vv2  = np.squeeze(V3d[:,jj2,ii2])
      vv3  = np.squeeze(V3d[:,jj3,ii3])
      vv4  = np.squeeze(V3d[:,jj4,ii4])

# Rotate onto N/E grid
      uNE1, vNE1 = rotate_uv(uu1, vv1, tht1)
      uNE2, vNE2 = rotate_uv(uu2, vv2, tht2)
      uNE3, vNE3 = rotate_uv(uu3, vv3, tht3)
      uNE4, vNE4 = rotate_uv(uu4, vv4, tht4)

      f_plt = False
      if f_plt:
        u1 = uu1[0]*100
        v1 = vv1[0]*100
        ur = uNE1[0]*100
        vr = vNE1[0]*100
        plot_vectors(u1, v1, ur, vr, nfig=1)


# Interpolate U V components (on true N/E grid) to NEP OB segments
# then rotate onto NEP I/J grid
# NEP angle is J/I to N/E not in math sense but positive clockwise!!!
# i.e. this is N/E to J/I angle in math sense
      UT  = np.array([uu1, uu2, uu3, uu4]).transpose()
      VT  = np.array([vv1, vv2, vv3, vv4]).transpose()
      Uintp = mblnr.bilin_interp1D(phi1, phi2, phi3, phi4, xht, yht, UT)
      Vintp = mblnr.bilin_interp1D(phi1, phi2, phi3, phi4, xht, yht, VT)
    
      tht_nep = thetaNEP[ikk]
      UR, VR = rotate_uv(Uintp, Vintp, tht_nep, ij2NE=False)
      #plot_vectors(Uintp[0]*100, Vintp[0]*100, UR[0]*100, VR[0]*100, nfig=1)

# Check: abs. values of interpolated values <= original data
      mxHT = np.max(abs(UT), axis=1)
      dmx  = (abs(Uintp))/mxHT
      if max(dmx-1.) > 0.1:
        print(f"!!! segm{nsgm} {varnm} Min/Max test violated: ikk={ikk} dlt: {np.max(dmx)}")
        if max(dmx-1.) > 1.:
          raise Exception("MinMax test: Interp error Check interpolation")

      if varnm == 'u':
        hintp = UR
      else:
        hintp = VR

      if ny == 1:
        Ai[itime, :, 0, ikk] = hintp
      else:
        Ai[itime, :, ikk, 0] = hintp
    
# Construct data set:
  if ny == 1:
    dz3d = np.expand_dims(dz3d, axis=2)
  else:
    dz3d = np.expand_dims(dz3d, axis=3)
  nxsgm = [x for x in range(0,nx)]
  nysgm = [x for x in range(0,ny)]
  stime = [x for x in range(1,ntime+1)]
  nzsgm = [x for x in range(0,nzlev)]
  sgmnm = f"segment_{nsgm:03d}"
  dim1  = "time"
  dim2  = f"nz_{sgmnm}"
  dim3  = f"ny_{sgmnm}"
  dim4  = f"nx_{sgmnm}"
  darr_var = xarray.DataArray(Ai, dims=(dim1, dim2, dim3, dim4),
             coords={dim1: stime, dim2: nzsgm, dim3: nysgm, dim4: nxsgm})
  darr_dz  = xarray.DataArray(dz3d, dims=(dim1, dim2, dim3, dim4),
             coords={dim1: stime, dim2: nzsgm, dim3: nysgm, dim4: nxsgm})
#  f"{varnm}_segment_{nsgm:03d}"
  dset_segm = xarray.Dataset({f"{varnm}_segment_{nsgm:03d}": darr_var, \
           f"dz_{varnm}_segment_{nsgm:03d}": darr_dz})

  return dset_segm
  
def rotate_uv(uu, vv, tht, ij2NE=True):
  """
    Rotate OLD axes to NEW axes by angle tht (positive in math sense), i.e.
    rotate vector to (-tht) vector U={uu,vv} ---> Ur={ur,vr}

    angle tht (in radians): angle between I/J axis and N/E true directions 
    positive / negative in math sense
    e.g. J-axis to true North angle tht > 0 
         means positive rotation of J to N
    ij2NE=False:         to rotate back from true North/East to I/J axes

    Use rotation matrix for (-tht):  cos(tht) sin(tht)
                                    -sin(tht) cos(tht)
  """
  if ij2NE:
    ur = uu*np.cos(tht) + vv*np.sin(tht)
    vr = uu*(-np.sin(tht)) + vv*np.cos(tht)
  else:
    ur = uu*np.cos(tht) - vv*np.sin(tht)
    vr = uu*np.sin(tht) + vv*np.cos(tht)

  return ur, vr

def get_rotangle(icegrid_spear, fconfig, grid_spear):
  """
    Derive angle between the model IJ-axes and true North-East
    from ice.static.nc file
    saved are components of the rotation matrix for vector rotation from
    the model I/J grid onto true N/E axes
    U_trueNE = R*U
    Derive angle as:
    arcsin(SINROT) ---> angle between J-axis and N (positive in math sense, J to N)
    arccos(COSROT) ---> angle between I-axis and E (sign of angle is not preserved)
  """
  from yaml import safe_load
  import matplotlib.pyplot as plt
  from matplotlib import colormaps
  from mod_utils_fig import bottom_text

  with open(fconfig) as ff:
    config = safe_load(ff)
  spear_dir = os.path.join(config['filesystem']['spear_month_ens'], 'monthly_clim')

  # Domain lon/lat boundaries:
  lonW = config['domain']['west_lon']
  lonE = config['domain']['east_lon']
  latN = config['domain']['north_lat']
  latS = config['domain']['south_lat']

  LON = icegrid_spear['GEOLON'].data
  LAT = icegrid_spear['GEOLAT'].data
  if lonW > np.max(LON):
    lonW = lonW-360.
  elif lonW < np.min(LON):
    lonW = lonW+360.

  if lonE > np.max(LON):
    lonE = lonE-360.
  elif lonE < np.min(LON):
    lonE = lonE+360.

  lon_slice = slice(lonW, lonE)
  lat_slice = slice(latS, latN)

  dset    = icegrid_spear.sel(yT=lat_slice, xT=lon_slice) 
  cosrot  = dset['COSROT'].data
  sinrot  = dset['SINROT'].data
  lon_sub = dset['GEOLON'].data
  lat_sub = dset['GEOLAT'].data
  Lmask   = grid_spear['wet'].sel(yh=lat_slice, xh=lon_slice).data

  theta_I2E = np.arccos(cosrot) # does not preserve the sign of the angle
  theta_J2N = np.arcsin(sinrot) 
  rad2deg   = 180./np.pi
  theta     = theta_J2N

  f_chck = False
  if f_chck:
    plt.ion()
    fig = plt.figure(1,figsize=(9,8))
    plt.clf()
    #ax1.cla()
    ax1 = plt.axes([0.1, 0.55, 0.8, 0.32])
    im=ax1.pcolormesh(theta_I2E*rad2deg, cmap='seismic', vmin=-60, vmax=60)
    ax1.contour(Lmask, [0.01], colors=[(0,0,0)])
  #  ax1.contour(cosrot, [x/10 for x in range(0,10)])
    ax1.contour(cosrot, [0.99], colors=[(0,0,1)])
#    ax1.contour(lon_sub,[x for x in range(-205,-100,5)], linestyles='solid', 
#                linewidths=0.5, colors=[(0.7,0.7,0.7)])
    ax1.contour(lat_sub,[x for x in range(10,88,5)], linestyles='solid', 
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
    ax1.axis('scaled')
    if lon_slice is None: ax1.set_ylim([150, 319])

    ax2 = plt.axes([0.91, 0.55, 0.03, 0.32])
    fig.colorbar(im, cax=ax2)
    ax1.set_title('Angle I to East from COSROT')

    ax3 = plt.axes([0.1, 0.08, 0.8, 0.32])
    im2=ax3.pcolormesh(theta_J2N*rad2deg, cmap='seismic', vmin=-60, vmax=60)
    ax3.contour(Lmask, [0.01], colors=[(0,0,0)])
    ax3.contour(sinrot, [0.01], colors=[(0,0,1)])
#    ax3.contour(lon_sub,[x for x in range(-205,-100,5)], linestyles='solid', 
#                linewidths=0.5, colors=[(0.7,0.7,0.7)])
    ax3.contour(lon_sub,[x for x in range(-300,-1,10)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
    ax3.contour(lon_sub,[x for x in range(0,100,10)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
    ax3.axis('scaled')
    if lon_slice is None: ax3.set_ylim([150, 319])

    ax4 = plt.axes([0.91, 0.08, 0.03, 0.32])
    fig.colorbar(im2, cax=ax4)
    ax3.set_title('Angle J to North from SINROT')

    btx = 'mod_utils_ob.py'
    bottom_text(btx)


  return theta, cosrot, sinrot 

def plot_rotangle(grid_spear, icegrid_spear, lon_slice=None, lat_slice=None):
  """
    Plot rot angles (cos/sine) from ice grid file
  """
  import matplotlib.pyplot as plt
  from matplotlib import colormaps
  from mod_utils_fig import bottom_text

  if not lon_slice is None:
    Lmask = grid_spear['wet'].sel(yh=lat_slice, xh=lon_slice).data
    cosrot = icegrid_spear['COSROT'].sel(yT=lat_slice, xT=lon_slice).data
    sinrot = icegrid_spear['SINROT'].sel(yT=lat_slice, xT=lon_slice).data
    lon_sub = icegrid_spear['GEOLON'].sel(yT=lat_slice, xT=lon_slice).data
    lat_sub = icegrid_spear['GEOLAT'].sel(yT=lat_slice, xT=lon_slice).data
  else:
    Lmask = grid_spear['wet'].data
    cosrot = icegrid_spear['COSROT'].data
    sinrot = icegrid_spear['SINROT'].data
    lon_sub = icegrid_spear['GEOLON'].data
    lat_sub = icegrid_spear['GEOLAT'].data


  theta = np.arcsin(sinrot)*180./np.pi  # degrees

  plt.ion()
  fig = plt.figure(1,figsize=(9,8))
  plt.clf()
  #ax1.cla()
  ax1 = plt.axes([0.1, 0.55, 0.8, 0.32])
  im=ax1.pcolormesh(cosrot, cmap='plasma', vmin=0.0, vmax=1)
  ax1.contour(Lmask, [0.01], colors=[(0,0,0)])
#  ax1.contour(cosrot, [x/10 for x in range(0,10)])
  ax1.contour(cosrot, [0.99], colors=[(0,0,1)])
  ax1.contour(lat_sub,[x for x in range(10,88,5)], linestyles='solid',     
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
  ax1.axis('scaled')
  if lon_slice is None: ax1.set_ylim([150, 319])

  ax2 = plt.axes([0.91, 0.55, 0.03, 0.32])
  fig.colorbar(im, cax=ax2)
  ax1.set_title('COSROT')


  ax3 = plt.axes([0.1, 0.08, 0.8, 0.32])
  im2=ax3.pcolormesh(sinrot, cmap='seismic', vmin=-1, vmax=1)
  ax3.contour(Lmask, [0.01], colors=[(0,0,0)])
  ax3.contour(sinrot, [0.01], colors=[(0,0,1)])
  ax3.contour(lon_sub,[x for x in range(-300,-1,10)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
  ax3.contour(lon_sub,[x for x in range(0,100,10)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
  ax3.axis('scaled')
  if lon_slice is None: ax3.set_ylim([150, 319])

  ax4 = plt.axes([0.91, 0.08, 0.03, 0.32])
  fig.colorbar(im2, cax=ax4)
  ax3.set_title('SINROT')

  btx = 'mod_utils_ob: plot_rotangle'
  bottom_text(btx)

  fig2 = plt.figure(2,figsize=(9,8))
  plt.clf()
  ax21 = plt.axes([0.1, 0.4, 0.8, 0.5])
  im2  = ax21.pcolormesh(theta, cmap='turbo', vmin=-90, vmax=90)
  ax21.contour(Lmask, [0.01], colors=[(0,0,0)])
  ax21.contour(sinrot, [0.01], colors=[(0,0,1)])
  ax21.contour(lat_sub,[x for x in range(10,88,5)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
  ax21.contour(lon_sub,[x for x in range(-300,-1,10)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
  ax21.contour(lon_sub,[x for x in range(0,100,10)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
  ax21.axis('scaled')
  if lon_slice is None: ax21.set_ylim([150, 319])

  ax22 = plt.axes([0.91, 0.4, 0.03, 0.5])
  fig2.colorbar(im2, cax=ax22)
  ax21.set_title('Angle: Axes to true N/E, degrees')

  bottom_text(btx)

def plot_rotangleNEP(hgrid, hmask, nfig=1):
  """
    Plot rot angles for NEP domain
  """
  import matplotlib.pyplot as plt
  from matplotlib import colormaps
  from mod_utils_fig import bottom_text

  Lmask = hmask['mask'].data
  theta = hgrid['angle_dx'].data*180./np.pi  # degrees
  ny_mask, nx_mask = Lmask.shape
  xmask = [x for x in range(0,2*nx_mask,2)]
  ymask = [x for x in range(0,2*ny_mask,2)]
  lon_nep = hgrid['x'].data
  lat_nep = hgrid['y'].data

  plt.ion()
  fig2 = plt.figure(nfig,figsize=(9,8))
  plt.clf()
  ax21 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im2  = ax21.pcolormesh(theta, cmap='turbo', vmin=0, vmax=90)
#  ax21.contour(Lmask, [0.01], colors=[(0,0,0)])
  ax21.contour(lat_nep,[x for x in range(5,88,5)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
  ax21.contour(lon_nep,[x for x in range(120,359,10)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
  ax21.axis('scaled')
#  if lon_slice is None: ax21.set_ylim([150, 319])

  ax22 = plt.axes([0.91, 0.4, 0.03, 0.5])
  fig2.colorbar(im2, cax=ax22)
  ax21.set_title('NEP Angle: Axes to true N/E, degrees')

  btx = 'mod_utils_ob.py: plot_rotangleNEP'
  bottom_text(btx)

def plot_vectors(u1, v1, ur, vr, nfig=1):
  """
    Quick plot 2 vectors to check rotation
  """
  plt.ion()
  fig = plt.figure(nfig,figsize=(9,8))
  fig.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  ax1.quiver(u1, v1, angles='xy', scale_units='xy', scale=1)
  ax1.quiver(ur, vr, angles='xy', scale_units='xy', scale=1, color='r')
  vmax = np.max(abs(np.array([u1,v1,ur,vr])))
  ax1.set_xlim([-vmax, vmax])
  ax1.set_ylim([-vmax, vmax])
  ax1.set_title('Original: black, Rotated: Red')

  return
 
