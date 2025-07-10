"""
  COnvert hdf4 grid file from SMOS/SMAP ice thickness website
  to netcdf


Huntemann, M., Heygster, G., Kaleschke, L., Krumpen, T., Mäkynen, M., and Drusch, M.: Empirical sea ice thickness retrieval during the freeze-up period from SMOS high incident angle observations, The Cryosphere, 8, 439-451, doi:10.5194/tc-8-439-2014, 2014.

"""
import os
import numpy as np
import xarray
from pyhdf.SD import SD, SDC
from netCDF4 import Dataset

# Input HDF4 file
grid_res = 12500  # 6250 or 12500
pthgrid = '/work/Dmitry.Dukhovskoy/data/SMOS_SMAP_thin_ice'
flnm_hdf = f'LongitudeLatitudeGrid-n{grid_res}-Arctic.hdf'
flnm_nc = f'LongitudeLatitudeGrid-n{grid_res}-Arctic.nc'
dfl_hdf = os.path.join(pthgrid,flnm_hdf)
dfl_nc  = os.path.join(pthgrid,flnm_nc)

# Open HDF4 file
hdf = SD(dfl_hdf, SDC.READ)
datasets = hdf.datasets()

#
datasets = hdf.datasets()
print("Datasets in file:")
for varname, info in datasets.items():
    print(f"{varname} shape: {info[1]}, type: {info[0]}")

darr_lon = None
darr_lat = None
for varname in datasets.keys():
  print(f"Processing dataset: {varname}")
  A2d = hdf.select(varname)[:]
  jdim,idim = A2d.shape
  X = np.arange(idim)
  Y = np.arange(jdim)
  if varname.lower().startswith('longit'):
    darr_lon = xarray.DataArray(A2d, dims=('Y', 'X'), coords={'Y': Y, 'X': X})
  elif varname.lower().startswith('latit'):
    darr_lat = xarray.DataArray(A2d, dims=('Y', 'X'), coords={'Y': Y, 'X': X})

dset_lonlat = xarray.Dataset({
  'Longitudes': darr_lon,
  'Latitudes': darr_lat
})

dset_lonlat.attrs["history"] = (
 f"Converted from hfd4 sis2_relax/convert_SMOS_lonlat_hdf4_to_netcdf.py"
)

print(f'Saving Lon/lat --> {dfl_nc}')
dset_lonlat.to_netcdf(
     dfl_nc,
     format='NETCDF4',
     engine='netcdf4'
)

print(f"Completed: {dfl_nc}")

