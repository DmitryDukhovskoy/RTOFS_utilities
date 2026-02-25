"""
  Interpolated ice fields from GFSv16 from regular Mercator grid --> mesh025 bipolar grid
  for N and S hemispheres bounded by 50 lat. N/S


"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
from yaml import safe_load
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, cm
import argparse
import time

# Append custom module paths
PPTHN = None
if 'PPTHN' not in locals() or PPTHN is None:
  cwd = os.getcwd()
  parts = cwd.split(os.sep)
  if 'python' in parts:
    idx = parts.index('python')
    PPTHN = os.sep + os.path.join(*parts[:idx + 1])
  else:
    raise RuntimeError("Directory 'python' not found in current working directory path.")

sys.path.extend([
    os.path.join(PPTHN, 'MyPython', 'hycom_utils'),
    os.path.join(PPTHN, 'MyPython', 'draw_map'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_regmom as mrmom

YR = 2025
MM = 12
DD = 31
FH = 0    # forecast hours
regn = 'global'
fload = 0

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help=f"year of GFSv16 run, default={YR}", type=int)
parser.add_argument("--mm", help=f"month of GFSv16 run, default={MM}", type=int)
parser.add_argument("--dd", help=f"motnh day of GFSv16 run, default={DD}", type=int)
parser.add_argument("--fhr", help=f"forecasts hour, default={FH}", type=int)
parser.add_argument("--fload", help=f"Load and finish tmp file 1=yes, 0=no, default=0", 
                               choices=[0,1], type=int)
args = parser.parse_args()

YR    = args.yr    if args.yr    else YR
MM    = args.mm    if args.mm    else MM
DD    = args.dd    if args.dd    else DD 
FH    = args.fhr   if args.fhr   is not None else FH 
fload = args.fload if args.fload is not None else fload

load_tmp = bool(fload)

syst_info = os.uname()
machine = syst_info.nodename

if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
elif 'an' in machine:
  print("Running on PPAN node:", machine)
  node_nm = "ppan"
else:
  print("Unknown machine:", machine)

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdmM, idmM = HH.shape   # mesh025

# Derive geo-coordinates from GFSv16 Polar Coordinates
pthgfs  = f'/gpfs/f6/gfs-cpu/world-shared/Lydia.B.Stefanova/forDmitry/GFSv16_icesubset/{YR}{MM:02d}{DD:02d}'
pthdump = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/gmapi_mapping2mesh025'

flgfs = f'gfsv16_icesubset.{YR}{MM:02d}{DD:02d}.fh{FH:03d}.nc'  
dflgfs = os.path.join(pthgfs,flgfs)
with xarray.open_dataset(dflgfs) as ds_gfs:
  lon1d = ds_gfs['longitude'].values
  lat1d = ds_gfs['latitude'].values

#idm_org = len(lon1d)
#jdm_org = len(lat1d)

# Add an extra row below 90N to avoid errors during findex search
# It also a good idea to avoid 90N lat to avoid errors in spehrical distance calculation
# Replace it with 89.999 
lat1d_org = lat1d.copy()
lat1d_north = np.insert(lat1d_org, -1, lat1d_org[-1] - 0.01)

# Add 2 halo-columns at 0 amd 360 longitudes 
# Create extended LON arrays for grid points near western bndry (x0 ~0) and east bndry (x0~360)
nhalo_i = 2
nhalo_j = 1
lon1d_org = lon1d.copy()
lon_west1  = lon1d_org[-1] - 360.
lon_west2  = lon1d_org[-2] - 360.
lon_east1  = lon1d_org[0]  + 360.
lon_east2  = lon1d_org[1]  + 360.
lon1d_west = np.insert(lon1d_org,  0, lon_west1)
lon1d_west = np.insert(lon1d_west, 0, lon_west2)
lon1d_east = np.append(lon1d_org,  lon_east1)
lon1d_east = np.append(lon1d_east, lon_east2)

LON, LAT = np.meshgrid(lon1d_org,lat1d_org)
LON_west, LAT_west = np.meshgrid(lon1d_west, lat1d_north)
LON_east, LAT_east = np.meshgrid(lon1d_east, lat1d_north)

jdm_org, idm_org = LON.shape
jdm_ext, idm_ext = LON_west.shape

# Find lat bounds of GFSv16 data:
# MOM6 points should be inside the GFSv16 domain 
# to be able to find 4-vertice of the bounding box for interpolation
lat_min = np.min(LAT)
lat_max = np.max(LAT)

ignore_north_lim = lat_max >= 90.
if ignore_north_lim:
  print(f"WARN: Indices north of northernmost lat={np.max(LAT):.4f} will be searched")

row_min = np.min(hlat, axis=1)  # min lat in each row
row_max = np.max(hlat, axis=1)  # max lat in each row
jS = np.argmax(row_min >= lat_min)
jE = len(row_max) - np.argmax(row_max[::-1] <= lat_max) - 1 # note reverse indexing for [::-1]

npnts_tot = (jE+1-jS)*(idmM)

print('Start searching gmapii ...')

fgmapi  = f'GFSv16_reggrid2mesh025_gmapi_{idmM}x{jdmM}_SNpoles.nc'
INDX_tmp = JNDX_tmp = IMOM_tmp = JMOM_tmp = None
prcs_pnts = None
if load_tmp:
  ftmp = os.path.splitext(fgmapi)[0] + '_tmp.nc'
  dftmp = os.path.join(pthdump, ftmp)

  if not os.path.isfile(dftmp):
    raise Exception(f"tmp gmapi file not found {dftmp}")

  with xarray.open_dataset(dftmp) as ds_tmp:
    INDX_tmp = ds_tmp['gmapi_i'].values
    JNDX_tmp = ds_tmp['gmapi_j'].values
    IMOM_tmp = ds_tmp['mom_indx'].values
    JMOM_tmp = ds_tmp['mom_jndx'].values
    jdim0    = ds_tmp.sizes['jdim']
    idim0    = ds_tmp.sizes['idim']

  assert jdm_org == jdim0 and idm_org == idim0, f"2D dims in {ftmp} does not match {jdm_org} x {idm_org}" 
    
  prcs_pnts = set(zip(IMOM_tmp.tolist(), JMOM_tmp.tolist()))


def convert_ihalo(nhalo_i, halo_grid, ixx0, idm_org):
  """
    Convert ixx halo indices back to the original array:
    for "west halos" 
        0         1    ...   nhalo_i-1 nhalo_i nhalo_i+1 ...           nhalo_i + idm-1  Halo indices
    -nhalo_i   -nhalo_i+1 ...    -1      0        1      2   3   4  ...  idm-1          Grid indices
         *        *               *      *        *      *   *   *       *
     <--        west halo       -->       <--  grid -->
             
    subtract indices of the west halos
    2 halo points on west/east are assumed
    West halos --> last original columns
  """
  if halo_grid == 'west':
    ixx = ixx0 - nhalo_i
    for inh in range(1,nhalo_i+1):
      ixx[ixx == -inh] = idm_org - inh

  elif halo_grid == 'east':
    # Convert ixx halo indices back to the original array:
    # East halos --> 1st/2nd original column
    ixx = ixx0.copy()
    for inh in range(nhalo_i):
      ixx[ixx0 == idm_org+inh] = inh

  return ixx

icc = -1
latS = -50.
latN = 50.
INDX = None
JNDX = None
IMOM = []
JMOM = []
start_time = time.perf_counter()
end_time = start_time
for ii in range(idmM):
  for jj in range(jS,jE+1):
    pnt_nmb = ii*(jE-jS+1) + jj + 1
    if pnt_nmb == 1 or pnt_nmb%50000 == 0:
      dt_step = (time.perf_counter() - end_time) / 60.
      end_time = time.perf_counter()
      dlt_time = (end_time - start_time) / 60.
      prcnt_done = float(pnt_nmb)/npnts_tot * 100.
      print(f'icc={icc}  Elapsed time={dlt_time:.1f} min, dt_step={dt_step:.3f} min, {prcnt_done:.2f}% done ...')

    if HH[jj,ii] >= 0:
      continue

    x0 = hlon[jj,ii]
    y0 = hlat[jj,ii]

    # Skip if outside the GFS domain
    if y0 < lat_min or y0 > lat_max:
      continue

    # Skip if outside polar regions:
    if y0 > latS and y0 < latN:
      continue 

    icc += 1

    # skip if already processed in tmp file:
    if load_tmp and (ii, jj) in prcs_pnts:
      idx = np.where((IMOM_tmp == ii) & (JMOM_tmp == jj))[0][0]
      ixx = INDX_tmp[idx]
      jxx = JNDX_tmp[idx]

    else:
      halo_grid = 'west'
      ixx0, jxx = mrmom.find_gridpnts_box(x0, y0, LON_west, LAT_west, dhstep=1., ignore_north_lim=True)

      if len(ixx0) == 0:
        # Try east halo grid:
        halo_grid = 'east'
        ixx0, jxx = mrmom.find_gridpnts_box(x0, y0, LON_east, LAT_east, dhstep=1., ignore_north_lim=True)

      ixx = convert_ihalo(nhalo_i, halo_grid, ixx0, idm_org)

      if len(ixx)==0 or len(jxx)==0:
        print(f"Box not found: ii={ii} jj={jj} x0={x0:.3f} y0={y0:.3f}")
        continue

    # Check that the ghost row is not in the box:
    if np.max(jxx) == jdm_org:
      jxx[jxx == jdm_org] = jdm_org - 1    

    ixx = np.expand_dims(ixx, axis=0)
    jxx = np.expand_dims(jxx, axis=0)

    if icc == 0:
      INDX = ixx.copy()
      JNDX = jxx.copy()
    else:
      INDX = np.append(INDX, ixx, axis=0)
      JNDX = np.append(JNDX, jxx, axis=0)

    IMOM.append(ii)
    JMOM.append(jj)

IMOM = np.array(IMOM)
JMOM = np.array(JMOM)

# GFSv16 coord dimensions:
jdim, idim = LON.shape

npnts = len(IMOM)
darr_imom = xarray.DataArray(IMOM, dims=("npoints"),\
                   coords={"npoints": np.arange(npnts)})
darr_jmom = xarray.DataArray(JMOM, dims=("npoints"),\
                   coords={"npoints": np.arange(npnts)})
darr_indx = xarray.DataArray(INDX, dims=("npoints","nvert"),\
                   coords={"npoints": np.arange(npnts),\
                           "nvert": np.arange(4)})
darr_jndx = xarray.DataArray(JNDX, dims=("npoints","nvert"),\
                   coords={"npoints": np.arange(npnts),\
                           "nvert": np.arange(4)})
darr_lon = xarray.DataArray(LON, dims=("jdim","idim"),\
                  coords={"jdim": np.arange(jdim),\
                          "idim": np.arange(idim)})
darr_lat = xarray.DataArray(LAT, dims=("jdim","idim"),\
                  coords={"jdim": np.arange(jdim),\
                          "idim": np.arange(idim)})

dset = xarray.Dataset({"mom_indx": darr_imom, \
                       "mom_jndx": darr_jmom, \
                       "gmapi_i": darr_indx,\
                       "gmapi_j": darr_jndx,\
                       "longit":  darr_lon,\
                       "latit":   darr_lat})

dset['mom_indx'].attrs['long_name'] = 'MOM6 grid I indices corresponding gmapi'
dset['mom_jndx'].attrs['long_name'] = 'MOM6 grid J indices corresponding gmapi'
dset['gmapi_i'].attrs['long_name'] = 'I indices GFSv16 grid for interpolation'
dset['gmapi_j'].attrs['long_name'] = 'J indices GFSv16 grid for interpolation'
dset['longit'].attrs['long_name']  = 'Longitudes derived from GFSv16 polar grid'
dset['latit'].attrs['long_name']   = 'Latitudes derived from GFSv16 polar grid'

# Global attributes
dset.attrs['title']       = 'Grid mapping between GFSv16 global regular Mercator grid and tripolar mesh025 grid'
dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
dset.attrs['source']      = 'get_gmapi_gfs16_reggrid_to_mesh025.py'
dset.attrs['history']     = 'GFSv16 NRT daily ice concentration, polar stereogr converted to geogr. coordinates'
dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
dset.attrs['region']      = f'S and N polar regions south of {latS:.1f} and north of {latN:.1f}'
dset.attrs['Grid_idm_jdm'] = f'{idmM}x{jdmM}'

dfgmapi = os.path.join(pthdump, fgmapi)

print(f'Saving gmapi --> {dfgmapi}')
dset.to_netcdf(dfgmapi, format='NETCDF4', engine='netcdf4')

f_chck = False
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:
  clrmp = mclrmps.colormap_conc()
  rmin = 0.
  rmax = 1.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  LON1 = LON.copy()
  LON2 = LON.copy()
  LON1 = np.where(LON1 < -175, np.nan, LON1)
  LON2 = np.where(LON2 > 172, np.nan, LON2)
  LON3 = np.where(LON < 0, LON+360., LON)
  LON3 = np.where(LON3 > 350., np.nan, LON3)
  lon_cntr1 = [x for x in range(-170,0,10)]  # grey -180:0
  lon_cntr2 = [x for x in range(10,178,10)]  # blue: 0 to 180 E
  if regn == 'south':
    lat_cntr = [x for x in range(-85,-20,5)]
  else:
    lat_cntr = [x for x in range(50,89,5)]

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.12, 0.12, 0.8, 0.8])

  if f_xy: 
    ax1.pcolormesh(XX,YY,C2d, cmap=clrmp, vmin=rmin, vmax=rmax)
    #ax1.invert_yaxis()
    # plot on XX,YY:
    cs = ax1.contour(XX,YY,LON1, lon_cntr1, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LON2, lon_cntr2, linestyles='solid', colors=[(0.,0.5,0.9)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LAT, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LON1,[0], linestyles='solid', colors=[(1,0.,0.)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LAT,[-75], linestyles='solid', colors=[(1,0.,0.)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LON3,[180], linestyles='solid', colors=[(0.6,0.,0.9)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    ax1.axis('scaled')
  else:
    # Plot on grid:
    ax1.pcolormesh(C2d, cmap=clrmp, vmin=rmin, vmax=rmax)
    cs = ax1.contour(LON1, lon_cntr1, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs2 = ax1.contour(LON2, lon_cntr2, linestyles='solid', colors=[(0.,0.5,0.9)], linewidths=1)
    ax1.clabel(cs2, inline=True, fontsize=10, fmt="%.1f")
    cs3 = ax1.contour(LON1,[0], linestyles='solid', colors=[(1,0.,0.)], linewidths=2)
    ax1.clabel(cs3, inline=True, fontsize=12, fmt="%.1f")
    ax1.contour(LON3,[180], linestyles='solid', colors=[(0.6,0.,0.9)], linewidths=1)
    cs = ax1.contour(LAT, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    ax1.axis('scaled')
    ax1.invert_yaxis() 
    ax1.set_ylabel('Inverted j index')
    ax1.set_xlabel('i index')

  ax1.set_title('Derived lon/lat from GFSv16 polar sterogr. projection')
  btx = 'get_gmapi_gfs16_reggrid_to_mesh025.py'
  bottom_text(btx) 


