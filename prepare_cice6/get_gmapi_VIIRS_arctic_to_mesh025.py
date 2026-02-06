"""
  Derive mapping indices gmapi to interpolate VIIRS high-res data --> mesh025 grid

  VIIRS ice thickness
  on native grid   
  daily and 4-day composites 
                    
  ice and snow depth, Arctic

https://coastwatch.noaa.gov/cwn/products/viirs-sea-ice-concentration-ice-thickness-ice-surface-temperature.html

Use NSIDC Polar Stereogr Proj:
semi_major_axis: 6378137.0
inverse_flattening: 298.257223563
straight_vertical_longitude_from_pole: -45.0
latitude_of_projection_origin: 90.0
standard_parallel: 70.0
false_easting: 0.0
false_northing: 0.0
srid: EPSG:3413
proj4text: +proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
crs_wkt: PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3413"]]


  Due to very large VIIRS grid, recommended to run several serial jobs saving temporary files with
  gmapi indices, using get_gmapi_subVIIRS_arctic_to_mesh025.py  e.g.
  run get_gmapi_subVIIRS_arctic_to_mesh025.py --iS0 0 --iE0 200
  run get_gmapi_subVIIRS_arctic_to_mesh025.py --iS0 201 --iE0 400
  ...

  then use this script to combine all saved pieces and creating the final netcdf

"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray
from copy import copy
import matplotlib.colors as colors
from yaml import safe_load
from mpl_toolkits.basemap import Basemap, cm
import argparse
import re    # regular expression module to search strings etc

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

import mod_time as mtime
from mod_utils_fig import bottom_text
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_regmom as mrmom
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

regn = 'north'
iS0 = 0
parser = argparse.ArgumentParser()
parser.add_argument("--usetmp", help="Load temporary INDX files, continue from saved (1=yes, 0=no)",
                                  choices=[0,1], required=True, type=int)
args = parser.parse_args()

stmp = args.usetmp if args.usetmp is not None else None
use_tmp = bool(stmp)

syst_info = os.uname()
machine = syst_info.nodename

if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
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

jdm, idm = hlon.shape


def extract_iS(fname):
  mnmb = re.search(r"_iS(\d+)_iE(\d+)", fname)  # extract numbers into 2 groups
  if mnmb is None:
    raise ValueError(f"Cannot parse iS/iE from {fname}")
  return int(mnmb.group(1))

def combine_tmp(varnm):
  """
    Combine saved tmp files with INDX, JNDX
  """
  # TMP files:
  pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
  pthtmp  = os.path.join(pthdata,'gmapi_NSIDC','tmp')

  fll_patt = os.path.join(pthtmp,f"VIIRS_gmapi_tmp_{varnm}_iS*_iE*.npy")
  files_tmp = glob.glob(fll_patt)

  if len(files_tmp) == 0:
    #raise FileNotFoundError(f"No TMP files found for {varnm}")
    print(f"No tmp files found for {varnm}")
    return []

  # Sort files by starting iS index
  files_tmp = sorted(files_tmp, key=extract_iS)

  blocks = []
  for ftmp in files_tmp:
    print(f"Reading {ftmp}")
    blocks.append(np.load(ftmp))

  # Combine appropriately
  if varnm in ("IMOM", "JMOM"):
    return np.concatenate(blocks)
  else:        
    # INDX, JNDX - numpy arrays
    return np.vstack(blocks)



def check_cntr(icc, ii, jj, iS0, jS0, iE0, jE0, ichk=200):
  if icc % ichk == 0:
    Npnts = (iE0 - iS0) * (jE0 - jS0)
    nIs_row = iE0 - iS0
    nJs_col = jE0 - jS0
    # Assuming loop is by i-indx first
    ncols_done = ii - iS0
    npnts_done = ncols_done * nJs_col + (jj - jS0 + 1)
    print(f' icc={icc} ii={ii} jj={jj} {(float(npnts_done)/Npnts*100.):.2f}% done ...')

# Get grid data using OpenDAP
YR=2021
MM=9
DD=15
dnmb = mtime.datenum([YR,MM,DD])
jday = int(mtime.date2jday([YR,MM,DD]))
#url = f"https://www.star.nesdis.noaa.gov/thredds/dodsC/CoastWatch/VIIRS/npp/IceThick/FourDaySector/NP06"
url = f"https://www.star.nesdis.noaa.gov/thredds/dodsC/IceThickVIIRSnppSectorFourDayNP06/{YR}/"
flviirs = f"VXSACW_B{YR}{jday-3:03d}_B{YR}{jday:03d}_H4_NP06_edgemask_IceThickness.nc"
dflviirs = os.path.join(url, flviirs)

print(f"Reading {dflviirs}")
with xarray.open_dataset(dflviirs) as ds_ices:
  LON = ds_ices['longitude'].values
  LAT = ds_ices['latitude'].values

# Convert to -180, 180:
LON = (LON + 180.) % 360. - 180.


# Find lat bounds of CryoSat data:
# MOM6 points should be inside the CryoSat domain 
# to be able to find 4-vertice of the bounding box for interpolation
lat_min = np.min(LAT)
lat_min = max([50.,lat_min])
lat_max = 90. # override np.max(LAT) for Polar stereogr. projection, the code should work
                # for correctly finding 4 points around hlat>np.max(LAT) but only
                # for Polar stereogr. projection by grabbing points over the N. Pole

ignore_north_lim = lat_max >= 90.  # make this true to avoid a discont. line across the Arctic Ocean
if ignore_north_lim:
  print(f"WARN: Indices north of northernmost lat={np.max(LAT):.4f} will be searched")

row_min = np.min(hlat, axis=1)
row_max = np.max(hlat, axis=1)
jS0 = np.argmax(row_min >= lat_min)
jE = len(row_max) - np.argmax(row_max[::-1] <= lat_max) - 1
#jE = np.argmin(row_max <= lat_max) - 1

INDX = None
JNDX = None
INDX_list = []
JNDX_list = []
IMOM = []
JMOM = []
icc = -1
iS = 0 
jS = jS0
iE = idm-1
Npnts = (iE+1) * ((jE+1) - jS)

print(f"Running gmapi i-indx={iS}:{iE}, j-indx={jS}:{jE}")

# Read saved tmp file if exist, and update indx arrays
if use_tmp:
  IMOM = combine_tmp('IMOM')
  JMOM = combine_tmp('JMOM')
  INDX = combine_tmp('INDX')
  JNDX = combine_tmp('JNDX')
  INDX_list = list(INDX)
  JNDX_list = list(JNDX)
  IMOM = list(IMOM)
  JMOM = list(JMOM)
  assert len(IMOM) == len(JMOM), f"IMOM and JMOM different lengths"
  assert INDX.shape[0] == JNDX.shape[0], f"INDX and JNDX different lengths"

# Set of processed points for checking:
pnts_done = set(zip(JMOM,IMOM))

Npnts = (iE + 1) * (jE + 1 - jS)
print(f"TMP updated: Running gmapi i-indx={iS}:{iE}, j-indx={jS}:{jE}")

if not use_tmp:
  # Very slow !!!
  for ii in range(iS, iE+1):
    if ii >= idm:
      continue
    for jj in range(jS, jE+1):
      if jj >= jdm:
        continue

      if (jj,ii) in pnts_done:
        icc += 1 
        check_cntr(icc, ii, jj, iS, jS, iE, jE, ichk=5000)
        continue

      if HH[jj,ii] >= 0:
        continue
      #ii, jj = mutil.find_indx_lonlat(358.047,89.816, hlon,hlat)
      x0 = hlon[jj,ii]
      y0 = hlat[jj,ii]
      if y0 < lat_min or y0 > lat_max:
        continue

      ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=0.1, ignore_north_lim=ignore_north_lim)
      ixx = np.asarray(ixx).reshape(-1)
      jxx = np.asarray(jxx).reshape(-1)

      check_cntr(icc, ii, jj, iS, jS, iE, jE, ichk=1000)
      if ixx.size == 0 or jxx.size == 0:
        continue

      icc += 1

      INDX_list.append(ixx.copy())
      JNDX_list.append(jxx.copy())
      IMOM.append(ii)
      JMOM.append(jj)

print("Main loop finished")

IMOM = np.array(IMOM)
JMOM = np.array(JMOM)
INDX = np.asarray(INDX_list)   # shape (N, 4)
JNDX = np.asarray(JNDX_list)

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
dset['gmapi_i'].attrs['long_name'] = 'I indices VIIRS grid for interpolation'
dset['gmapi_j'].attrs['long_name'] = 'J indices VIIRS grid for interpolation'
dset['longit'].attrs['long_name']  = 'Longitudes VIIRS polar grid'
dset['latit'].attrs['long_name']   = 'Latitudes VIIRS polar grid'

# Global attributes
dset.attrs['title']       = 'Grid mapping between VIIRS on NSIDC Polar Stereographic North WGS84 and MOM6 mesh025 grids'
dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
dset.attrs['source']      = 'get_gmapi_VIIRS_arctic_to_mesh025.py'
dset.attrs['history']     = 'Arctic sea ice thickness from VIIRS 750m L3 product, NOAA/NESDIS/STAR/SOCD'
dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
dset.attrs['region']      = regn

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump = os.path.join(pthdata,'gmapi_NSIDC')
fgmapi  = f'VIIRS_ithkn_MOM6_gmapi_{idm}x{jdm}_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)

print(f'Saving gmapi --> {dfgmapi}')
dset.to_netcdf(dfgmapi, format='NETCDF4', engine='netcdf4')


f_chck = False
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.12, 0.12, 0.8, 0.8])

  # Plot on grid:
  #img = ax1.pcolormesh(C2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  # Check longitudes:

  JD, ID = zip(*pnts_done)

  ax1.cla()
  ax1.contour(HH, [0], linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
  ax1.plot(ID, JD, '.', markersize=1)


  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  if rmin < 0:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  else:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')

  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'get_gmapi_CryoSat_to_mesh025.py'
  bottom_text(btx) 



