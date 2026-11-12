"""
Find ERA5 closest index for mesh025 SFS ocean grid: Artic or Antarctic
using nearest neighbor
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
from yaml import safe_load

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

from mod_utils_fig import bottom_text, colorbar_horiz
import mod_time as mtime
import mod_mom6 as mmom6

def save_netcdf(JERA, IERA, JGLR, IGLR, LONERA, LATERA, LONGLR, LATGLR, dfgmapi):
  """
  NetCDF writer
  """
  # Target indices on ERA subset grid
  JERA = np.asarray(JERA, dtype=int)
  IERA = np.asarray(IERA, dtype=int)
  # mesh025 indices where field to be interpolated:
  JGLR = np.asarray(JGLR, dtype=int)
  IGLR = np.asarray(IGLR, dtype=int)

  Npnts = len(JGLR)

  assert len(JERA) == len(IERA) == len(JGLR) == len(IGLR), \
       "Inconsistent number of mapping points"

  # ERA5 coord:
  jdimE = LATERA.shape[0]
  idimE = LONERA.shape[0]

  # mesh025 coord:
  jdimG, idimG = LATGLR.shape

  darr_lonE = xr.DataArray(LONERA, dims="idimE",\
                           coords={"idimE": np.arange(idimE)})
  darr_latE = xr.DataArray(LATERA, dims="jdimE",\
                           coords={"jdimE": np.arange(jdimE)})
  darr_lonG = xr.DataArray(
      LONGLR, dims=("jdimG", "idimG"),
      coords={"jdimG": np.arange(jdimG),
              "idimG": np.arange(idimG)}
  )

  darr_latG = xr.DataArray(
      LATGLR, dims=("jdimG", "idimG"),
      coords={"jdimG": np.arange(jdimG),
              "idimG": np.arange(idimG)}
  )
  coords = {"npoints": np.arange(Npnts)}

  darr_iera = xr.DataArray(IERA, dims="npoints", coords=coords)
  darr_jera = xr.DataArray(JERA, dims="npoints", coords=coords)
  darr_iglr = xr.DataArray(IGLR, dims="npoints", coords=coords)
  darr_jglr = xr.DataArray(JGLR, dims="npoints", coords=coords)
  
  dset = xr.Dataset({"era_indx":      darr_iera, \
                     "era_jndx":      darr_jera, \
                     "mesh025_indx":   darr_iglr, \
                     "mesh025_jndx":   darr_jglr, \
                     "era_longit":    darr_lonE, \
                     "era_latit":     darr_latE, \
                     "mesh025_longit": darr_lonG, \
                     "mesh025_latit":  darr_latG})


  dset['era_indx'].attrs['long_name']      = 'ERA5 grid I index'
  dset['era_jndx'].attrs['long_name']      = 'ERA5 grid J index'
  dset['mesh025_indx'].attrs['long_name']   = 'mesh025 grid I index'
  dset['mesh025_jndx'].attrs['long_name']   = 'mesh025 grid J index'
  dset['era_longit'].attrs['long_name']    = 'ERA5 longitudes'
  dset["era_longit"].attrs["units"]        = "degrees_east"
  dset['era_latit'].attrs['long_name']     = 'ERA5 latitudes'
  dset["era_latit"].attrs["units"]         = "degrees_north"
  dset['mesh025_longit'].attrs['long_name'] = 'mesh025 longitudes'
  dset['mesh025_latit'].attrs['long_name']  = 'mesh025 latitudes'
  dset["mesh025_longit"].attrs["units"]     = "degrees_east"
  dset["mesh025_latit"].attrs["units"]      = "degrees_north"

  # Global attributes:
  dset.attrs['title']  = 'ERA5-to-UFS mesh025 grid mapping using nearest neighbors'
  dset.attrs['source'] = 'find_remap_indx_era5_to_mesh025.py'
  dset.attrs['index_base'] = 0

  print(f'Saving gmapi --> {dfgmapi}')
  dset.to_netcdf(dfgmapi, format='NETCDF4', engine='netcdf4')




parser = argparse.ArgumentParser()
parser.add_argument(
  "--debug", 
  help="Run in debug mode to check cKDTRee method vs brute force nearest dist",
  choices=[0,1], 
  type=int, 
  default=0
)
args = parser.parse_args()

run_debug = args.debug == 1
regn  = "global"
f_save = True

rdate = 20250101

dnmb = mtime.rdate2datenum(rdate)
YR, MM, DD = mtime.datevec(dnmb)[:3]

fyaml = "config_MLithkn_oras5.yaml"
with open(fyaml) as ff:
  config_ml = safe_load(ff)

pthdata   = config_ml["MOM6"]["pthdata"]  # root dir for processed data
#pthindx = os.path.join(pthdata,'remap_indx')
pthdinx = config_ml["ERA5"]["pthgmapi"]
pthera5 = config_ml["ERA5"]["pthera"].format(YR=YR)
flbase  = 'ERA5_forcing_'

if regn == 'north':
  lat0 = 50
elif regn == 'south':
  lat0 = -50
elif regn == 'global':
  lat0 = 90

dgr2rad = np.pi/180.

# Read ERA5 grid
flera5 = f"{flbase}{YR}-{MM:02d}-{DD:02d}.nc"
dflera5 = os.path.join(pthera5, flera5)

with xr.open_dataset(dflera5) as ds:
  LON_era = ds["longitude"].values
  LAT_era = ds["latitude"].values

LON_era = (LON_era + 360) % 360

if lat0 is None:
  if regn == 'north':
    lat0 = np.min(LAT_era)
  elif regn == 'south':
    lat0 = np.max(LAT_era)
  elif regn == 'global':
    lat0 = 90.

# Radians:
LONR_era, LATR_era = np.meshgrid(LON_era * dgr2rad, LAT_era * dgr2rad)

# Convert to spherical coord:
Z_era = np.sin(LATR_era)
X_era = np.cos(LATR_era) * np.cos(LONR_era)
Y_era = np.cos(LATR_era) * np.sin(LONR_era)

# Read mesh025 grid:
pthgrid = config_ml["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xr.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, 1., -HH)

jdm, idm = HH.shape
LMsk = HH < 0
use_lmask = True

# Regional domain:
if regn == 'north':
  Ireg = hlat >= lat0
elif regn == 'south' or regn == 'global':
  Ireg = hlat <= lat0

J, I = np.where(Ireg)
Npnts = np.shape(J)[0]


# Use cKDTree for quicker search of the closest points:
from scipy.spatial import cKDTree

xyz_era = np.column_stack((
  X_era.ravel(),
  Y_era.ravel(),
  Z_era.ravel()
))

tree = cKDTree(xyz_era)

JGL = []    # UFS mesh025 global index j where fields need to be interpolated
IGL = []    # UFS mesh025 global index i
JER = []    # ERA5 j index of the saved domain (typically not global)
IER = []    # ERA5 i index of the saved domain
cntr = 0
dmin0 = 0.5   # max allowed distance (degrees) to the closest point

for jj, ii in zip(J,I):
  cntr += 1

  if cntr%50000 == 0:
    print(f' cntr={cntr} {cntr/Npnts*100:.2f}% done ...')

  if use_lmask and not LMsk[jj, ii]:
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
      print(f"Different indices: jmin={jmin} jbf={jbf}, imin={imin} ibf={ibf} but")
      print(f"KDTree dot = {dot_tree:.16f}  " 
            f"Brute   dot = {dot_bf:.16f}")


  # Regitser
  JER.append(jmin)
  IER.append(imin)
  JGL.append(jj)
  IGL.append(ii)

if f_save:
  flout = f"gmapi_ERA5_to_UFSmesh025_{regn}.nc"
  dflout = os.path.join(pthindx, flout)

  save_netcdf(JER, IER, JGL, IGL, LON_era, LAT_era, hlon, hlat, dflout)




