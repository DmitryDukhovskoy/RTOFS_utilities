"""
  OI SST high resolution fields
  https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
# OpenDap to PSL data does not work on PPAN
# it works on Gaea
#
# So 2 options: 
# (1) read SST on Gaea and extract NEP domain - save netcdf --> PPAN
#
# (2) copy files from PSL Linux via Niagara
# [ddukhovskoy@linux256 noaa.oisst.v2.highres]$ pwd
#/Datasets/noaa.oisst.v2.highres
# ---> Niagara untrusted ---> Gaea --- gcp ---> PPAN

 Use this script on Gaea to extract NEP domain
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
from netCDF4 import Dataset as ncFile
import importlib
import xarray
import yaml
from yaml import safe_load

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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_anls_seas as manseas

YR = 1994

def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  dmm = np.copy(dmm0)
  return dmm


# OpenDap to PSL data does not work on PPAN
#urlBase = 'http://psl.noaa.gov/thredds/dodsC/Datasets/noaa.oisst.v2.highres/'
#urlT    = f'sst.day.mean.{YR}.nc'
#furl = os.path.join(urlBase,urlT)
#T2d  = read_field(furl,'SST')
#ds = xarray.open_dataset(furl, chunks={})
#ds_sst = xarray.open_dataset(furl)

pthsst = '/work/Dmitry.Dukhovskoy/data/OISST/'
flsst  = os.path.join(pthsst,f'sst.day.mean.{YR}.nc')
print(f'Subsetting {YR} {flsst}')
ds_sst = xarray.open_dataset(flsst)
lon1d  = ds_sst['lon'].data
lat1d  = ds_sst['lat'].data

xlim1 = 155.
xlim2 = 260.
ylim1 = 10.
ylim2 = 84.

lon1d = np.where(lon1d<0., lon1d+360., lon1d)
DD    = np.abs(lon1d - xlim1)
ii1   = np.argmin(DD)
DD    = np.abs(lon1d - xlim2)
ii2   = np.argmin(DD)
DD    = np.abs(lat1d - ylim1)
jj1   = np.argmin(DD)
DD    = np.abs(lat1d - ylim2)
jj2   = np.argmin(DD)

SST   = ds_sst['sst'].data[:,jj1:jj2+1, ii1:ii2+1].squeeze()
lon   = lon1d[ii1:ii2+1]
lat   = lat1d[jj1:jj2+1]
idm   = len(lon)
jdm   = len(lat)
TM    = ds_sst['time'].data
ntm   = len(TM)  
 
darr_lon = xarray.DataArray(lon, dims=('lon'), coords={'lon': lon})
darr_lat = xarray.DataArray(lat, dims=('lat'), coords={'lat': lat})
darr_tm  = xarray.DataArray(TM, dims=('time'), coords={'time': TM})
darr_sst = xarray.DataArray(SST, dims=('time','lat','lon'), 
                            coords={'time': TM, 'lat': lat, 'lon': lon})

dset_sub = xarray.Dataset({'time': TM, 'LON': darr_lon, 'LAT': darr_lat, 'sst': darr_sst})
dset_sub.attrs["title"] = 'NOAA Daily Optimum Interpolation Sea Surface Temperature'
dset_sub.attrs["References"] = 'https://www.psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html'
dset_sub.attrs['Info'] = 'Subset for Northeast Pacific domain'
dset_sub.attrs['code'] = '/home/Dmitry.Dukhovskoy/python/setup_seasonal_NEP/extract_OISST_NEPdomain.py'

dset_sub['sst'].attrs['long_name'] = 'Daily Mean Sea Surface Temperature'
dset_sub['sst'].attrs['units'] = 'degC'


flsst_out = os.path.join(pthsst, f'oisst_dayily_NEPsubset_{YR}.nc')
print(f'Saving ---> {flsst_out}')
dset_sub.to_netcdf(flsst_out,
                   format='NETCDF3_64BIT',
                   engine='netcdf4',
                   unlimited_dims='time')








