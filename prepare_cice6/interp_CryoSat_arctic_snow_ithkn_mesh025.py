"""
  Interpolate CryoSat hsnow or ice thickn. monthly fileds 2018-2021
  winter months only

  NSIDC data
  Monthly fields for three Arctic growth seasons (October to April) from 2018 to 2021.
  Sahra Kacimi and Ron Kwok

  Data are on polar stereographic coordinates

  gmapi indices: get_gmapi_CryoSat_arcticNSIDC_to_mesh025.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib  
import xarray
from copy import copy
import matplotlib.colors as colors 
from yaml import safe_load
from mpl_toolkits.basemap import Basemap, cm
import pandas as pd
import argparse
                   
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
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

fsave = 1
regn = 'north'

parser = argparse.ArgumentParser()
parser.add_argument("--field", help=f"Field to interpolate: snow thkciness or ice thickn",
                    choices=['sndpth','ithkn'], required=True, type=str)
parser.add_argument("--fsave", help=f"Save interp. monthly fields netcdf, default={fsave}", 
                    choices=[0,1], type=int)
args = parser.parse_args()
  
field_name = args.field if args.field else None
fsave = args.fsave if args.fsave is not None else fsave

save_nc = fsave == 1
save_tmp = False

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

jdm, idm = HH.shape
LMsk = np.where(HH<0, 1, 0)

# Mask out not needed latitudes:
LMsk = np.where(hlat < 50, 0, LMsk)

# Get gmapi 4 NSIDC grid points for interpolation
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump  = os.path.join(pthdata,"gmapi_NSIDC")
fgmapi  = f'CryoSat_NSIDC_MOM6_gmapi_1440x1080_north.nc'
dfgmapi = os.path.join(pthdump, fgmapi)
print(f'Loading gmapi --> {dfgmapi}')

with xarray.open_dataset(dfgmapi) as dgmapi:
  IMOM = dgmapi['mom_indx'].data
  JMOM = dgmapi['mom_jndx'].data
  INDX = dgmapi['gmapi_i'].data
  JNDX = dgmapi['gmapi_j'].data

with xarray.open_dataset(dfgmapi) as dgmapi: 
  LON = dgmapi['longit'].data
  LAT = dgmapi['latit'].data

icc = 0
A3d = np.zeros((jdm,idm))
print(f"Saving netcdf: {save_nc}")

# Find N of records in snow/ice file
# Note only winter months exist
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthsnow = os.path.join(pthdata,'CryoSat_arctic_ice_snow_thkn')
flsnow = 'NSIDC-0773_SD-THK_25km_IS2-CS2_ArcticGrowthSeasons2018-2021_v01.nc'
dflsnow = os.path.join(pthsnow, flsnow)

with xarray.open_dataset(dflsnow) as dsn:
  time = dsn.time
years = time.dt.year.values
months = time.dt.month.values
nrec = len(months)

if field_name == 'sndpth':
  varnm = 'sd'
elif field_name == 'ithkn':
  varnm = 'thk'

for irec in range(nrec):
  YR = years[irec]
  MM = months[irec]
  print(f"Processing {YR}/{MM:02d} ...")

  if field_name == 'sndpth':
    attr_str = 'snow depth'
    fliceout = f'hsnow_NSIDC_CryoSat_arctic_interp_mesh025_{jdm}x{idm}_{YR}{MM:02d}.nc'
  elif field_name == 'ithkn':
    attr_str = 'ice thickness'
    fliceout = f'ithkn_NSIDC_CryoSat_arctic_interp_mesh025_{jdm}x{idm}_{YR}{MM:02d}.nc'

  pthnsidc = os.path.join(pthdata,'CryoSat_arctic_ice_snow_thkn/interp_NSIDC_monthly')
  dfliceout = os.path.join(pthnsidc,fliceout)

  # Skip already saved files:
  if os.path.isfile(dfliceout): 
    print(f"File exists: {dfliceout}, skipping ...\n")
    continue

  cff_m = None
  with xarray.open_dataset(dflsnow) as dsn:
    AA = dsn[varnm].isel(time=irec).data.squeeze()
    units = dsn[varnm].attrs.get('units', None)
    if units == 'cm':
      cff_m =0.01      # cm --> m
    elif units == 'm':
      cff_m = 1.

  AA = np.where(AA > 9999., np.nan, AA) * cff_m  # cm --> m  
 
  HSint = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)

  # Fill the N.Pole hole:
  HSint = mrmom.fill_npole(HSint, hlon, hlat, HH, Rpole=2.)

  A3d = np.where(HH>=0, np.nan, HSint)
  A3d = np.expand_dims(A3d, axis=0) 

  time_out = np.array([np.datetime64(f"{YR:04d}-{MM:02d}-01", "ns")])
  darr_cice = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
                     coords={"time": time_out,\
                             "jdim": np.arange(jdm),\
                             "idim": np.arange(idm)})
  if field_name == 'sndpth':
    dset = xarray.Dataset({"snow_depth": darr_cice})
    dset['snow_depth'].attrs['long_name']='snow depth on ice'
    dset['snow_depth'].attrs['units']='m'
  elif field_name == 'ithkn':
    dset = xarray.Dataset({"ice_thkn": darr_cice})
    dset['ice_thkn'].attrs['long_name']='ice thickness'
    dset['ice_thkn'].attrs['units']='m'

  dset["time"].attrs = {
       "long_name": "time"
  }

  # Add global attributes:
  dset.attrs['title']       = f'Arctic {attr_str} from ICESat-2 and CryoSat-2' 
  dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
  dset.attrs['source']      = 'interp_CryoSat_arctic_snow_ithkn_mesh025.py'
  dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
  dset.attrs['region']      = 'north'

  print(f'Dumping interpolated {field_name} --> {dfliceout}\n')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')



