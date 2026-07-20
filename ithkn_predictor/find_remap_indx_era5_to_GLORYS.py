"""
  Find ERA5 closest index for GLORYS grid: Artic or Antarctic
  using nearest neighbor

  Downloaded every 7-day daily mean 2m SAT for specified region from ERA5 website
  https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics?tab=download



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

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="Region: north or south", choices=['north','south'], 
                    required=True, type=str)
parser.add_argument("--debug", help="Run in debug mode to check cKDTRee method vs brute force nearest dist",
                    choices=[0,1], type=int, default=0)
args = parser.parse_args()
run_debug = args.debug == 1
regn = args.regn

pthindx = '/archive/Dmitry.Dukhovskoy/data/remap_indx'

f_save = True


if regn == 'north':
  #lat0 = 60
  pthera5 = '/archive/Dmitry.Dukhovskoy/data/ERA5/Arctic'
  flbase = 'era5_2mTemp_daily7day_Arctic_'
else:
  lat0 = -55

dgr2rad = np.pi/180.

# Read ERA5 grid:
YR=2000
flera5 = f"{flbase}{YR}.nc"
dflera5 = os.path.join(pthera5,flera5)

with xr.open_dataset(dflera5) as ds:
  LON_era = ds.longitude.values
  LAT_era = ds.latitude.values

LON_era = (LON_era + 360) % 360

if regn == 'north':
  lat0 = np.min(LAT_era)
elif regn == 'south':
  lat0 = np.max(LAT_era)

LONR_era, LATR_era = np.meshgrid(LON_era * dgr2rad, LAT_era * dgr2rad)

# Convert to spherical coord:
Z_era = np.sin(LATR_era)
X_era = np.cos(LATR_era) * np.cos(LONR_era)
Y_era = np.cos(LATR_era) * np.sin(LONR_era)



# Read GLORYS grid:
pthice = f"/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick/{YR}"
rdate = f"{YR*10000+100+1}"

# Find file:
def find_file(rdate, pthice):
  from pathlib import Path
  try:
    dflice = next(
        Path(pthice).glob(
            f"*_mean_{rdate}_R*.nc"
        )
    )
    print(f"Found file: {dflice}")
    return dflice
  except StopIteration:
    print(f"No file found for {rdate} in {pthice}")


#dflice = os.path.join(pthice, flice)
dflice = find_file(rdate, pthice)
if not os.path.isfile(dflice):
  raise RuntimeError(f"File not found: {dflice}")

with xr.open_dataset(dflice) as dsice:
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values

# 0 <= lon < 360
LON = (LON + 360) % 360
hlon, hlat = np.meshgrid(LON, LAT)

LMsk = None
use_lmask = True
if use_lmask:
  pthssh = f"/uda/Global_Ocean_Physics_Reanalysis/global/daily/zos/{YR}"
  dflssh = find_file(rdate, pthssh)
  with xr.open_dataset(dflssh) as dszos:
    SSH = dszos['zos'].isel(time=0).values.squeeze()

  LMsk = np.where(np.isfinite(SSH),1,0)

# Regional domain:
if regn == 'north':
  Ireg = hlat >= lat0
elif regn == 'south':
  Ireg = hlat <= lat0

J, I = np.where(Ireg)
Npnts = np.shape(J)[0]


def save_netcdf(JERA, IERA, JGLR, IGLR, LONERA, LATERA, LONGLR, LATGLR, dfgmapi, regn):
  """
    NetCDF writer
  """
  # Target indices on ERA subset grid
  JERA = np.asarray(JERA, dtype=int)
  IERA = np.asarray(IERA, dtype=int)
  # GLORYS indices where field to be interpolated:
  JGLR = np.asarray(JGLR, dtype=int)
  IGLR = np.asarray(IGLR, dtype=int)

  Npnts = len(JGLR)

  # ERA5 coord:
  jdimE = LATERA.shape[0]
  idimE = LONERA.shape[0]

  # GLORYS coord:
  jdimG = LATGLR.shape[0]
  idimG = LONGLR.shape[0]

  darr_lonE = xr.DataArray(LONERA, dims="idimE",\
                           coords={"idimE": np.arange(idimE)})
  darr_latE = xr.DataArray(LATERA, dims="jdimE",\
                           coords={"jdimE": np.arange(jdimE)})
  darr_lonG = xr.DataArray(LONGLR, dims="idimG",\
                           coords={"idimG": np.arange(idimG)})
  darr_latG = xr.DataArray(LATGLR, dims="jdimG",\
                           coords={"jdimG": np.arange(jdimG)})
  darr_iera = xr.DataArray(IERA, dims="npoints",\
                           coords={"npoints": np.arange(Npnts)})
  darr_jera = xr.DataArray(JERA, dims="npoints",\
                           coords={"npoints": np.arange(Npnts)})
  darr_jglr = xr.DataArray(JGLR, dims="npoints",\
                           coords={"npoints": np.arange(Npnts)})
  darr_iglr = xr.DataArray(IGLR, dims="npoints",\
                           coords={"npoints": np.arange(Npnts)})

  dset = xr.Dataset({"era_indx":      darr_iera, \
                     "era_jndx":      darr_jera, \
                     "glorys_indx":   darr_iglr, \
                     "glorys_jndx":   darr_jglr, \
                     "era_longit":    darr_lonE, \
                     "era_latit":     darr_latE, \
                     "glorys_longit": darr_lonG, \
                     "glorys_latit":  darr_latG})


  dset['era_indx'].attrs['long_name']      = 'ERA5 grid I index'
  dset['era_jndx'].attrs['long_name']      = 'ERA5 grid J index'
  dset['glorys_indx'].attrs['long_name']   = 'GLORYS grid I index'
  dset['glorys_jndx'].attrs['long_name']   = 'GLORYS grid J index'
  dset['era_longit'].attrs['long_name']    = 'ERA5 longitudes'
  dset["era_longit"].attrs["units"]        = "degrees_east"
  dset['era_latit'].attrs['long_name']     = 'ERA5 latitudes'
  dset["era_latit"].attrs["units"]         = "degrees_north"
  dset['glorys_longit'].attrs['long_name'] = 'GLORYS longitudes'
  dset['glorys_latit'].attrs['long_name']  = 'GLORYS latitudes'
  dset["glorys_longit"].attrs["units"]     = "degrees_east"
  dset["glorys_latit"].attrs["units"]      = "degrees_north"

  # Global attributes:
  dset.attrs['title']    = 'Grid mapping from ERA5 --> GLORYS, closest neighbour'
  dset.attrs['source']   = 'find_remap_indx_era5_to_GLORYS.py'
  dset.attrs['region']   = regn 

  print(f'Saving gmapi --> {dfgmapi}')
  dset.to_netcdf(dfgmapi, format='NETCDF4', engine='netcdf4')


# Use cKDTree for quicker search of the closest points:
from scipy.spatial import cKDTree

xyz_era = np.column_stack((
  X_era.ravel(),
  Y_era.ravel(),
  Z_era.ravel()
))

tree = cKDTree(xyz_era)

JGL = []    # GLORYS global index j where fields need to be interpolated
IGL = []    # GLORYS global index i
JER = []    # ERA5 j index of the saved domain (typically not global)
IER = []    # ERA5 i index of the saved domain
cntr = 0
dmin0 = 0.5   # max allowed distance (degrees) to the closest point
for jj, ii in zip(J,I):
  cntr += 1

  if cntr%50000 == 0:
    print(f' cntr={cntr} {cntr/Npnts*100:.2f}% done ...')

  if use_lmask and LMsk[jj, ii] == 0:
    JGL.append(jj)
    IGL.append(ii)
    JER.append(-999)
    IER.append(-999)
    continue

  x0 = hlon[jj, ii]
  y0 = hlat[jj, ii]

  # Make sure longitude convention matches ERA5
  if x0 < 0.:
    x0 += 360.

  # Unit vector in spherical coord
  x0rad = np.deg2rad(x0)
  y0rad = np.deg2rad(y0)

  x0plr = np.cos(y0rad) * np.cos(x0rad)
  y0plr = np.cos(y0rad) * np.sin(x0rad)
  z0plr = np.sin(y0rad)

  dist, idx = tree.query((x0plr, y0plr, z0plr))
  jmin, imin = np.unravel_index(idx, X_era.shape)

  # Brute force, slow - computes N_Era_points per 1 glorys point
  # Dot product with every ERA5 grid point
  if run_debug:
    DOT = (
      X_era * x0plr +
      Y_era * y0plr +
      Z_era * z0plr
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
      X_era[jmin, imin] * x0plr +
      Y_era[jmin, imin] * y0plr +
      Z_era[jmin, imin] * z0plr
    )

    dot_bf = DOT[jbf, ibf]

    assert np.isclose(dot_tree, dot_bf, atol=1.e-14)

    if (jmin != jbf) or (imin != ibf):
      print("Different indices but")
      print(f"KDTree dot = {dot_tree:.16f}")
      print(f"Brute   dot = {dot_bf:.16f}")


  # Regitser
  JER.append(jmin)
  IER.append(imin) 
  JGL.append(jj)
  IGL.append(ii)

if f_save:
  flout = f"gmapi_ERA5_to_GLORYS_{regn}.nc"
  dflout = os.path.join(pthindx, flout)

  save_netcdf(JER, IER, JGL, IGL, LON_era, LAT_era, LON, LAT, dflout, regn)



