"""
  Extract domain from SPEAR daily SSH fields
  individual ensembles
  save by variables

  I use pp_ens - postprocessed files separated by variables

  To run the script use either this script directly or for multiple years/ months:
  scritps/seasonal_fcst/subset_spear_ocean.sh

"""
import numpy as np
import os
from pathlib import Path
import importlib
import spear_path
from spear_path import get_spear_paths
import subprocess
import xarray
import argparse
from yaml import safe_load

pthtmp = '/work/Dmitry.Dukhovskoy/tmp/'
YR     = 2020
mstart = 6
ens    = 10
varnm  = 'ssh'
#varnm  = 'SSH'
regn   = 'NEP'
tmpdir = os.path.join(pthtmp,f"{YR}",f"ens{ens:02d}")

prefix = 'ocean_daily'

#flspear = f"{YR}{ens:02d}01.ocean_z_month.nc"
flspear = spear_path.get_spear_file(YR, mstart, prefix, 'daily', varnm, fstr=True)
finpnc  = os.path.join(tmpdir,flspear)
fconfig = 'config_nep.yaml'

with open(fconfig) as ff:
  config = safe_load(ff)

# Check lon convention:
LON = xarray.open_dataset(finpnc).variables['xh'].data 

lonW = config['domain']['west_lon']
lonE = config['domain']['east_lon']
if lonW > np.max(LON):
  lonW = lonW-360.
elif lonW < np.min(LON):
  lonW = lonW+360.

if lonE > np.max(LON):
  lonE = lonE-360.
elif lonE < np.min(LON):
  lonE = lonE+360.

# Note conversion from [-180, 180] to [0, 360]
#lon_slice = slice(config['domain']['west_lon'] % 360, config['domain']['east_lon'] % 360)
lon_slice = slice(lonW, lonE)
lat_slice = slice(config['domain']['south_lat'], config['domain']['north_lat'])
#outp_dir = Path(config['filesystem']['spear_month_ens']) / f"{YR}" /f"ens{ens:02d}"
outp_dir = os.path.join(config['filesystem']['spear_month_ens'], f"{YR}",f"ens{ens:02d}")
subprocess.run([f'mkdir -pv {outp_dir}'], shell=True)

print(f"Extracting SPEAR ocean {varnm}  for NEP, {YR}, ensemble={ens}")
dset = (xarray.open_dataset(finpnc).sel(yh=lat_slice, xh=lon_slice))
#    dset = dset.rename({'xh': 'lon'})
#    dset = dset.rename({'yh': 'lat'})

#dset = dset.swap_dims({'time': 'lead'})
#dset = dset.rename({'time': 'valid_time'})
# convert calendar from julian to normal gregorian
#dset['valid_time'] = (['lead'], dset['valid_time'].to_index().to_datetimeindex())

foutp = os.path.join(outp_dir, f"{regn}_spear_{YR}{mstart:02d}.ssh_daily.nc")

print(f'Saving --> {foutp}')
dset.to_netcdf(foutp)


