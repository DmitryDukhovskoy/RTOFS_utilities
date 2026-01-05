"""
  Interpolate snow climatology from EWG snow data
  to MOM6/CICE6 mesh025 grid
  2 options: snow depth, liquid water equivalent

  gmapi indices: get_gmapi_EWG_snow_Arctic_to_mesh025.py

  From EWG 
  https://nsidc.org/data/search#keywords=Arctic+snow+climatology/sortKeys=score,,desc/facetFilters=%257B%257D/pageNumber=1/itemsPerPage=25


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
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

fsave = 1
regn = 'north'

parser = argparse.ArgumentParser()
parser.add_argument("--snfld", help=f"snow field to interpolate",
                    choices=['sndpth','swe'], required=True, type=str)
parser.add_argument("--fsave", help=f"Save final dataset with all days as netcdf, default={fsave}", 
                    choices=[0,1], type=int)
args = parser.parse_args()
  
snfld = args.snfld if args.snfld else None
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
fgmapi  = f'EWG_Atlas_MOM6_gmapi_1440x1080_north.nc'
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

# Save temporary fields 
if save_tmp:
  tmp_dir = os.path.join(pthdata,'NRT_NOAA_NSIDC_seaconc','tmp')
  os.makedirs(tmp_dir, exist_ok=True)

icc = 0
A3d = np.zeros((12,jdm,idm))
print(f"Saving netcdf at the end: {save_nc}")

for imo in range(1,13):
  print(f"Processing month {imo} ...")

  pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
  pthsnow = os.path.join(pthdata,'Warren_snow_clim_EWG_atlas/DATA/GRIDDED_FIELDS/SNOW_DEPTH')
  if snfld == 'sndpth':
    flsnow = f'snow_depth.{imo}.1954_1991.dat'
  elif snfld == 'swe':
    flsnow = f'swe.{imo}.1954_1991.dat'

  dflnm = os.path.join(pthsnow,flsnow)
  print(f"Reading {dflnm}")

  data_hsnow = np.loadtxt(dflnm)

  # sanity check
  assert data_hsnow.size == 23 * 23, f"Check file length expected {23*23} lines"

  AA = data_hsnow.reshape((23, 23))
  AA = np.where(AA > 9999., np.nan, AA) * 0.01  # cm --> m  
 
  HSint = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
  HSint = np.where(HH>=0, np.nan, HSint)
  A3d[imo-1,:,:] = HSint

time_months = [x for x in range(1,13)]
if not save_nc:
  print(f"Final netcdf is not saved, save_nc={save_nc}")
 
else:
  darr_cice = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
                     coords={"time": time_months,\
                             "jdim": np.arange(jdm),\
                             "idim": np.arange(idm)})
  if snfld == 'sndpth':
    dset = xarray.Dataset({"snow_depth": darr_cice})
    dset['snow_depth'].attrs['long_name']='snow depth'
    dset['snow_depth'].attrs['units']='m'
  elif snfld == 'swe':
    dset = xarray.Dataset({"swe": darr_cice})
    dset['swe'].attrs['long_name']='snow water equivalent'
    dset['swe'].attrs['units']='m'

  # Add global attributes:
  dset.attrs['title']       = 'hsnow EWG Atlas clim interpolated onto mash025 grid'
  dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
  dset.attrs['source']      = 'interp_EWG_snow_mesh025.py'
  dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
  dset.attrs['region']      = regn
  dset.attrs['Grid_idm_jdm'] = f'{idm}x{jdm}'

  if snfld == 'sndpth':
    fliceout = f'hsnow_EWGatlas_interp_mesh025_{jdm}x{idm}_{regn}.nc'
  elif snfld == 'swe':
    fliceout = f'swe_EWGatlas_interp_mesh025_{jdm}x{idm}_{regn}.nc'


  pthnsidc = os.path.join(pthdata,'Warren_snow_clim_EWG_atlas/snow_clim_interp')
  dfliceout = os.path.join(pthnsidc,fliceout)
  print(f'Dumping interpolated ice conc --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')



