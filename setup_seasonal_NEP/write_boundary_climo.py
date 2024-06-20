"""
  Write nc files with SPEAR climatology OBC
  for seasonal forecasts

  Not finished, need to transfer NEP obc files from gaea to PPANLS
  

"""
import numpy as np
import xarray
from utils import modulo, smooth_climatology
import argparse
from pathlib import Path
from yaml import safe_load

#parser = argparse.ArgumentParser()
#parser.add_argument('-c', '--config')
#args = parser.parse_args()
#with open(args.config, 'r') as file: 
#    config = safe_load(file)

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)


first_year = config['climatology']['first_year']
last_year = config['climatology']['last_year']                   
pathin = Path(config['filesystem']['open_boundary_files'])
pathout = Path(config['filesystem']['model_input_data']) / 'boundary' / f'climatology_{first_year}_{last_year}'
pathout.mkdir(exist_ok=True, parents=True)
write_boundary(first_year, last_year, pathin, pathout)
  


def write_boundary(ystart, yend, pathin, pathout):
  for var in ['zos', 'thetao', 'so', 'uv']:
    for segment in [1, 2, 3]:
      print(f'{var} {segment}')
      boundary = xarray.open_dataset(pathin / f'{var}_{segment:03d}.nc')
      boundary = boundary.sel(time=slice(str(ystart), str(yend)))
      # To be sure
      assert int(boundary['time.year'].min()) == ystart
      assert int(boundary['time.year'].max()) == yend
      if var == 'uv':
        vardata = boundary[[f'u_segment_{segment:03d}', f'v_segment_{segment:03d}']]
      else:
        vardata = boundary[f'{var}_segment_{segment:03d}']
      ave = vardata.groupby('time.dayofyear').mean('time').sel(dayofyear=slice(1, 365))
      smoothed = smooth_climatology(ave).rename({'dayofyear': 'time'})

      encoding = {
          'time': dict(_FillValue=1.0e20),
          f'lon_segment_{segment:03d}': dict(dtype='float64', _FillValue=1.0e20),
          f'lat_segment_{segment:03d}': dict(dtype='float64', _FillValue=1.0e20),
          f'{var}_segment_{segment:03d}': dict(_FillValue=1.0e20),
      }

      if var == 'zos':
        # zos doesn't have z coordinates to worry about
        res = smoothed.to_dataset()
      else:
        # z coordinates don't really vary in time. use the first coord and expand over time.
        # do it for both u and v if it is a velocity file.
        if var == 'uv':
          z = boundary[[f'dz_u_segment_{segment:03d}', f'dz_v_segment_{segment:03d}']].isel(time=0).drop('time').expand_dims(time=365)
          encoding = {
              'time': dict(_FillValue=1.0e20),
              f'lon_segment_{segment:03d}': dict(dtype='float64', _FillValue=1.0e20),
              f'lat_segment_{segment:03d}': dict(dtype='float64', _FillValue=1.0e20),
              f'u_segment_{segment:03d}': dict(_FillValue=1.0e20),
              f'v_segment_{segment:03d}': dict(_FillValue=1.0e20)
          }
        else:
          z = boundary[f'dz_{var}_segment_{segment:03d}'].isel(time=0).drop('time').expand_dims(time=365)

        z['time'] = smoothed['time']
        res = xarray.merge([smoothed, z])

      for coord in ['lat', 'lon']:
        fullcoord = f'{coord}_segment_{segment:03d}'
        if fullcoord not in res:
          res[fullcoord] = boundary[fullcoord]

      res = modulo(res)
      res.to_netcdf(
          pathout / f'{var}_c_{segment:01d}.nc',
          format='NETCDF3_64BIT',
          engine='netcdf4',
          encoding=encoding,
          unlimited_dims='time'
      )


