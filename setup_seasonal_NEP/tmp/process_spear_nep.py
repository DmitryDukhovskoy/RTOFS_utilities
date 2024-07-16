"""
  Script to extract monthly or daily SPEAR data for individual ensemble members
  Original version written by Andrew C. Ross
  Modified for my needs

  For many fields, better to use sbatch:
  write_spear_atmos_job.sh 
  # Usage: sbatch write_spear_atmos_job.sh YEAR MONTH ENSEMBLE CONFIG

"""
import datetime as dt
from glob import glob
import numpy as np
from os import path
import xarray
from yaml import safe_load
from pathlib import Path, PurePath

import argparse
parser = argparse.ArgumentParser()
#parser.add_argument('-d', '--domain')
#parser.add_argument('-f', '--freq')
parser.add_argument('-v', '--var')  # 
#parser.add_argument('-e', '--ensemble')
#parser.add_argument('-c', '--config', default=None)
#args = parser.parse_args()

# Top level path to all SPEAR medium reforecast data on archive
ROOT = Path('/archive') / 'l1j' / 'spear_med' / 'rf_hist' / 'fcst' / 's_j11_OTA_IceAtmRes_L33'
root = ROOT

freq    = 'daily' # daily or monthly
fconfig = 'config_nep.yaml'
ens     = 
with open(fconfig) as ff:
  config = safe_load(ff)
  xslice = (config['domain']['west_lon'], config['domain']['east_lon'])
  yslice = (config['domain']['south_lat'], config['domain']['north_lat'])

def prepro(ds):
    ds['start'] = (('start', ), [dt.datetime(int(ds['time.year'][0]), int(ds['time.month'][0]), 1)])
    ds['lead'] = (('time', ), np.arange(len(ds['time'])))
    ds = ds.swap_dims({'time': 'lead'})
    ds = ds.rename({'time': 'valid_time'})
    # convert calendar from julian to normal gregorian
    # do it here while valid_time is 1D
    ds['valid_time'] = (['lead'], ds['valid_time'].to_index().to_datetimeindex())
    return ds


def slice_ds(ds, xslice, yslice):
    if xslice is None and yslice is None:
        return ds
    
    slice_dict = {}
    if xslice is not None:
        x = xslice.split(',')
        for xcoord in ['xh', 'xq', 'xT']:
            if xcoord in ds.coords:
                slice_dict.update({xcoord: slice(float(x[0]), float(x[1]))})
    if yslice is not None:
        y = yslice.split(',')
        for ycoord in ['yh', 'yq', 'yT']:
            if ycoord in ds.coords:
                slice_dict.update({ycoord: slice(float(y[0]), float(y[1]))})
    return ds.sel(**slice_dict)


def process_monthly(root, domain, var, ens=None, xslice=None, yslice=None):
    files = sorted(glob(path.join(root, f'{domain}.*-*.{var}.nc')))
    processed = xarray.open_mfdataset(files, preprocess=prepro, combine='nested', concat_dim='start')[var]
    processed = slice_ds(processed, xslice, yslice)

    if ens != 'pp_ensemble':
        fname = path.join(root, f'{domain}.monthly_mean.ens_{int(ens):02d}.{var}.nc')
    else:
        fname = path.join(root, f'{domain}.monthly_mean.ensmean.{var}.nc')
    processed.to_netcdf(fname)
    print(fname)


def process_daily(root, domain, var, ens=None, xslice=None, yslice=None):
    files = sorted(glob(path.join(root, f'{domain}.*-*.{var}.nc')))
    processed = xarray.open_mfdataset(files, preprocess=prepro, combine='nested', concat_dim='start')[var]
    processed = slice_ds(processed, xslice, yslice)

    if ens != 'pp_ensemble':
        fname = path.join(root, f'{domain}.daily_mean.ens_{int(ens):02d}.{var}.nc')
    else:
        fname = path.join(root, f'{domain}.daily_mean.ensmean.{var}.nc')
    processed.to_netcdf(fname)
    print(fname)


print(f"Processing {freq} {args.var} {ens}")
if freq == 'monthly':
  process_monthly(root, args.domain, args.var, ens=args.ensemble, xslice=xslice, yslice=yslice)
elif freq == 'daily':
  process_daily(root, args.domain, args.var, ens=args.ensemble, xslice=xslice, yslice=yslice)


