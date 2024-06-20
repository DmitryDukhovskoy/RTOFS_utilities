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

def load_var(spear_dir, varnm, yr1=1993, yr2=2020, mstart=1, fzint=True):
  flname = f'spear_moclim.mstart{mstart:02d}.{yr1}-{yr2}.{varnm}.nc'
  dflinp = os.path.join(spear_dir, flname)

  if not varnm == 'ssh':
    spear_var = xarray.open_dataset(dflinp).rename({'z_l': 'z'})
  else:
    spear_var = xarray.open_dataset(dflinp)

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

#def interp2line(AF,

def derive_obsegment(hgrid, ds, segments, isgm, varnm, INDX, JNDX):
  """
    Derive data set for an OB segment = isgm for variable = varnm
    ds - xarray of climatology field for 1 variable, subset to the region
    hgrid - topo/grid for the climatology field
    segments - xarray with information for OB segments created in boundary.Segment
    isgm - segment # [0, 1, ...]

    INDX, JNDX - n x 4 arrays of grid points (gmapi) for bilinear interpolation

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
   
  
