"""
  Find mesh025 closest index for GLORYS grid: Artic or Antarctic
  using nearest neighbor

  Find mapping GLORYS(Jg,Ig) <--> mesh025(Jm,Im)

  This is 2-way mapping, can be use for "interpolating"
  from / to mesh025 to/from GLORYS

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray as xr
import pandas as pd
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, cm
from yaml import safe_load
import argparse
from pathlib import Path
import urllib.request


#ROOT = Path(__file__).resolve().parent

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

import mod_glorys as mglr


parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="Region: north or south", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--debug", help="Run in debug mode to check cKDTRee method vs brute force nearest dist",
                    choices=[0,1], type=int, default=0)
args = parser.parse_args()

run_debug = args.debug == 1
regn = args.regn
f_save = True

regions = {
    "north": ("Arctic", 65.0),
    "south": ("Antarctic", -60.0),
}
regn_name, lat0 = regions[regn]

dgr2rad = np.pi/180.

pthindx = '/archive/Dmitry.Dukhovskoy/data/remap_indx'

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

DIRS = {
  "pthithkn" : config_predictor["linregr"]["pthithkn"],
  "pthiconc" : config_predictor["linregr"]["pthiconc"],
  "pthsst"   : config_predictor["linregr"]["pthsst"],
  "pthssh"   : config_predictor["linregr"]["pthssh"],
  "ptht2m"   : config_predictor["linregr"]["ptht2m"].format(regn_name=regn_name),
  "pthout"   : config_predictor["linregr"]["pthout"],
  }

# Read mesh025 grid
pthdata = '/archive/Dmitry.Dukhovskoy/data'
pthice    = os.path.join(pthdata, 'ithkn_clim_combined')
fliceout  = 'ithkn_mnthclim_cryo_avhrr_ices_1440x1080_north.nc'
dfliceout = os.path.join(pthice,fliceout)

with xr.open_dataset(dfliceout) as dsice:
  LON_m25 = dsice['lon'].data
  LAT_m25 = dsice['lat'].data

# Subset grid for faster search algorithm
jlat0 = None
if regn == 'north':
  lat_row = np.nanmax(LAT_m25, axis=1)
  jlat0 = np.where(lat_row < lat0)[0][-1]
  subLON_m25 = LON_m25[jlat0:,:]
  subLAT_m25 = LAT_m25[jlat0:,:]
elif regn == 'south':
  lat_row = np.nanmin(LAT_m25, axis=1)
  jlat0 = np.where(lat_row > lat0)[0][0]
  subLON_m25 = LON_m25[:jlat0+1,:]
  subLAT_m25 = LAT_m25[:jlat0+1,:]


# Convert to spherical coord:
Z_m25 = np.sin(subLAT_m25 * dgr2rad)
X_m25 = np.cos(subLAT_m25 * dgr2rad) * np.cos(subLON_m25 * dgr2rad)
Y_m25 = np.cos(subLAT_m25 * dgr2rad) * np.sin(subLON_m25 * dgr2rad)


# Read GLORYS grid
YR = 1993
pthice = f"/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick/{YR}"
rdate = f"{YR*10000+100+1}"
 
dflice = mglr.find_file(rdate, pthice)
if not os.path.isfile(dflice):
  raise RuntimeError(f"File not found: {dflice}")

with xr.open_dataset(dflice) as dsice:
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values

# 0 <= lon < 360
LON = (LON + 360) % 360
hlon, hlat = np.meshgrid(LON, LAT)

# Land mask:
LMsk = None
pthssh = os.path.join(DIRS["pthssh"],f"{YR}")
dflssh = mglr.find_file(rdate, pthssh)
with xr.open_dataset(dflssh) as dszos:
  SSH = dszos['zos'].isel(time=0).values.squeeze()

LMsk = np.where(np.isfinite(SSH),1,0)

# Define domain:
if regn == 'north':
  DOMAIN = (hlat > lat0) & (LMsk == 1)
elif regn == 'south':
  DOMAIN = (hlat < lat0) & (LMsk == 1)

J, I = np.where(DOMAIN)
Npnts = np.shape(J)[0]


def save_netcdf(JM25, IM25, JGLR, IGLR, LONm025, LATm025, LONGLR, LATGLR, dfgmapi, regn):
  """
    NetCDF writer
  """
  # Target indices on mesh025 grid
  JM25 = np.asarray(JM25, dtype=int)
  IM25 = np.asarray(IM25, dtype=int)
  # GLORYS indices where field to be interpolated:
  JGLR = np.asarray(JGLR, dtype=int)
  IGLR = np.asarray(IGLR, dtype=int)

  Npnts = len(JGLR)
  assert len(IM25) == Npnts, f"Nmb of mapped points on mesh025 {len(IM25)} does not match GLORYS={Npnts}"

  # dimensions of mesh025 grid
  jdimE, idimE = LATm025.shape
  # Subset lon/lat from mesh025 for closest neighbour points:
  lon_m025 = LONm025[JM25,IM25]
  lat_m025 = LATm025[JM25,IM25]

  # GLORYS coord:
  jdimG = LATGLR.shape[0]
  idimG = LONGLR.shape[0]

  pnts_arr = np.arange(Npnts)

  darr_lonM = xr.DataArray(lon_m025, dims="npoints",\
                           coords={"npoints": pnts_arr})
  darr_latM = xr.DataArray(lat_m025, dims="npoints",\
                           coords={"npoints": pnts_arr})
  darr_lonG = xr.DataArray(LONGLR, dims="idimG",\
                           coords={"idimG": np.arange(idimG)})
  darr_latG = xr.DataArray(LATGLR, dims="jdimG",\
                           coords={"jdimG": np.arange(jdimG)})
  darr_im25 = xr.DataArray(IM25, dims="npoints",\
                           coords={"npoints": pnts_arr})
  darr_jm25 = xr.DataArray(JM25, dims="npoints",\
                           coords={"npoints": pnts_arr})
  darr_jglr = xr.DataArray(JGLR, dims="npoints",\
                           coords={"npoints": pnts_arr})
  darr_iglr = xr.DataArray(IGLR, dims="npoints",\
                           coords={"npoints": pnts_arr})

  dset = xr.Dataset({"mesh025_indx":   darr_im25, \
                     "mesh025_jndx":   darr_jm25, \
                     "glorys_indx":    darr_iglr, \
                     "glorys_jndx":    darr_jglr, \
                     "mesh025_longit": darr_lonM, \
                     "mesh025_latit":  darr_latM, \
                     "glorys_longit":  darr_lonG, \
                     "glorys_latit":   darr_latG})

  dset['mesh025_indx'].attrs['long_name']   = 'UFS mesh 0.25-degree grid I index'
  dset['mesh025_jndx'].attrs['long_name']   = 'UFS mesh 0.25-degree grid J index'
  dset['glorys_indx'].attrs['long_name']    = 'GLORYS grid I index'
  dset['glorys_jndx'].attrs['long_name']    = 'GLORYS grid J index'
  dset['mesh025_longit'].attrs['long_name'] = 'UFS mesh 0.25-degree longitudes at I,J points'
  dset["mesh025_longit"].attrs["units"]     = "degrees_east"
  dset['mesh025_latit'].attrs['long_name']  = 'UFS mesh 0.25-degree latitudes at I,J points'
  dset["mesh025_latit"].attrs["units"]      = "degrees_north"
  dset['glorys_longit'].attrs['long_name']  = 'GLORYS longitudes'
  dset['glorys_latit'].attrs['long_name']   = 'GLORYS latitudes'
  dset["glorys_longit"].attrs["units"]      = "degrees_east"
  dset["glorys_latit"].attrs["units"]       = "degrees_north"

  # Global attributes:
  dset.attrs['title']    = 'Grid mapping from/to UFS 0.25-dgr grid to/from GLORYS, {regn}, closest neighbour'
  dset.attrs['source']   = 'find_remap_indx_GLORYS_to_mesh025.py'
  dset.attrs['region']   = regn

  print(f'Saving gmapi --> {dfgmapi}')
  dset.to_netcdf(dfgmapi, format='NETCDF4', engine='netcdf4')

# Use cKDTree for quicker search of the closest points:
from scipy.spatial import cKDTree

xyz_m25 = np.column_stack((
  X_m25.ravel(),
  Y_m25.ravel(),
  Z_m25.ravel()
))

tree = cKDTree(xyz_m25)

JGL = []    # GLORYS global index j where fields need to be interpolated
IGL = []    # GLORYS global index i
JM25 = []    # mesh025 j index
IM25 = []    # mesh025 i index 
cntr = 0
dmin0 = 0.5   # max allowed distance (degrees) to the closest point
print("Searching closest neighbours ....")
for jj, ii in zip(J,I):
  cntr += 1

  if cntr%25000 == 0:
    print(f' cntr={cntr} {cntr/Npnts*100:.2f}% done ...')

  if LMsk[jj, ii] == 0:
    JGL.append(jj)
    IGL.append(ii)
    JM25.append(-999)
    IM25.append(-999)
    continue

  # GLORYS:
  x0 = hlon[jj, ii]
  y0 = hlat[jj, ii]

  x0 = (x0 + 360) % 360.

  # Unit vector in spherical coord
  x0rad = np.deg2rad(x0)
  y0rad = np.deg2rad(y0)

  x0plr = np.cos(y0rad) * np.cos(x0rad)
  y0plr = np.cos(y0rad) * np.sin(x0rad)
  z0plr = np.sin(y0rad)

  dist, idx = tree.query((x0plr, y0plr, z0plr))
  jmin, imin = np.unravel_index(idx, X_m25.shape)

  # Brute force, slow - computes N_m25_points per 1 glorys point
  # Dot product with every mesh025 grid point
  if run_debug:
    DOT = (
      X_m25 * x0plr +
      Y_m25 * y0plr +
      Z_m25 * z0plr
    )

    # For unit vector in sph. coord, the min arc distance d = R*Theta
    # u * v = |u||v|cos(theta), minimize d --> maximize dot product, i.e. theta = 0
    jbf, ibf = np.unravel_index(np.argmax(DOT), DOT.shape)

    # Angular distance (degrees)
    dmin_degr = np.degrees(np.arccos(np.clip(DOT[jbf, ibf], -1., 1.)))

    # Check: convert distance in radian to degrees:
    assert dmin_degr < dmin0, f"Min distance {Dmin_degr:.5f} > {dmin0:.5f} dgr"

    # Check KDTree and distance:
    dot_tree = (
      X_m25[jmin, imin] * x0plr +
      Y_m25[jmin, imin] * y0plr +
      Z_m25[jmin, imin] * z0plr
    )

    dot_bf = DOT[jbf, ibf]

    assert np.isclose(dot_tree, dot_bf, atol=1.e-14)

    if (jmin != jbf) or (imin != ibf):
      print(f"jj={jj}, ii={ii}: Different indices but can be 2 equally distant points:")
      print(f"jmin / imin = {jmin}/{imin}, jbf / ibf = {jbf} / {ibf}")
      print(f"KDTree dot = {dot_tree:.16f}, Brute   dot = {dot_bf:.16f}")

  # Convert subset indices back to global:
  jMmin = jmin + (jlat0 if regn == "north" else 0)
  iMmin = imin

  # Check final indices:
  x0m25 = (LON_m25[jMmin, iMmin] + 360) % 360  
  y0m25 = LAT_m25[jMmin, iMmin]
  dlon = abs((x0m25 - x0 + 180.0) % 360.0 - 180.0)
  dlat = abs(y0m25 - y0)
  dlt_dgr = max(dlon, dlat)
  km1dgr = 111.
  dlt_km = np.sqrt((np.cos(y0m25 * dgr2rad)*km1dgr*dlon)**2 + (km1dgr*dlat)**2)
  dlt_dgr_max = 1.  # does not work well near pole - close points may have big long. difference
  dlt_km_max = dmin0 * km1dgr

  if dlt_km > dlt_km_max:
    print(f"jj={jj}, ii={ii}: dlt degrees btw closest pnt and x0,y0 too big={dlt_dgr:.7f}")
    print(f"distance btw the closest point and x0, y0 = {dlt_km}")
    print(f"x0={x0}, y0={y0} Found pnt: x={x0m25}, y={y0m25}")
    raise RuntimeError(f"dlt_km {dlt_km} > {dlt_km_max}")

  # Regitser
  JM25.append(jMmin)    # mesh025
  IM25.append(iMmin)    # mesh025
  JGL.append(jj)        # GLORYS
  IGL.append(ii)        # GLORYS

if f_save:
  flout = f"gmapi_closenghb_mesh025_to_GLORYS_{regn}.nc"
  dflout = os.path.join(pthindx, flout)

  save_netcdf(JM25, IM25, JGL, IGL, LON_m25, LAT_m25, LON, LAT, dflout, regn)




