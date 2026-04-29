"""
  Interpolate RTOFS CICE4 initial fields
  with corrected north pole seam 
  onto 0.25 grid

  see:
  get_gmapi_RTOFS_to_mesh025.py
  cice_rtofs/correct_RTOFS_NPole_seam.py

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray
import matplotlib.colors as colors 
from yaml import safe_load
from mpl_toolkits.basemap import Basemap, cm
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

regn = 'global'
init = 20250704
init_hr = 0
fhr  = 0      # RTOFS f/cast hours
field_name = 'ithkn' 

parser = argparse.ArgumentParser()
parser.add_argument("--field", help=f"Field to interpolate",
                    choices=['hsnow','rhosn','ithkn','iconc'], required=True, type=str)
parser.add_argument("--init", help="RTOFS init date", 
                    choices=[20250704], required=True, type=int)
args = parser.parse_args()
  
field_name = args.field 
init_date = args.init 


dnmbI = mtime.rdate2datenum(init_date*100+init_hr)
YR, MM, DD = mtime.datevec(dnmbI)[:3]

if dnmbI < mtime.datenum([2025,8,1]):
  rtofs_vers = "2.4"
else:
  rtofs_vers = "2.5"

print(f"\nInterpolating {field_name} RTOFS CICE4 {init_date} --> mesh025\n")

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

# Convert to negative depths:
if np.nanmin(HH) > -1e-6:
  HH = np.where(HH < 1.e-6, np.nan, HH) # assuming land ~0
  HH = -HH
  HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape

# Mask out not needed latitudes:
skip_mask = (hlat > -50) & (hlat < 50)
skip_mask &= (HH >= 0)
LMsk = np.ones_like(HH)
LMsk[skip_mask] = 0


print(f"Interpolating  {field_name} {init_date}")
A3d = np.zeros((jdm,idm))

# Get gmapi 4 NSIDC grid points for interpolation
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump  = os.path.join(pthdata,"gmapi_NSIDC")
fgmapi  = 'RTOFS_CICE4_gmapi_1440x1080_global.nc'
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


if field_name == 'ithkn':
  varnm = 'hi'
  str_long = 'ice thickness grid cell mean'
elif field_name == 'hsnow':
  varnm = 'hs'
  str_long = 'grid cell mean snow thickness over sea ice'
elif field_name == 'iconc':
  varnm = 'aice'
  str_long = 'aggregated ice partial area'
else:
  raise Exception(f"Unrecognized variable {field_name}")

def write_nc(A2d, time_out, field_name, dfliceout, varnm, rtofs_vers):
  model_str = f'RTOFSv{rtofs_vers} CICE4' 
  # Dump netcdf:
  darr_cice = xarray.DataArray(A2d, dims=("time","jdim","idim"),\
                     coords={"time": time_out,\
                             "jdim": np.arange(jdm),\
                             "idim": np.arange(idm)})
  if field_name == 'hsnow':
    dset = xarray.Dataset({"snow_depth": darr_cice})
    dset['snow_depth'].attrs['long_name'] = 'snow depth on ice'
    dset['snow_depth'].attrs['units'] = 'm'
  elif field_name == 'ithkn':
    dset = xarray.Dataset({"ice_thkn": darr_cice})
    dset['ice_thkn'].attrs['long_name'] = 'ice thickness'
    dset['ice_thkn'].attrs['units'] = 'm'
  elif field_name == 'iconc':
    dset = xarray.Dataset({"ice_conc": darr_cice})
    dset['ice_conc'].attrs['long_name'] = 'ice partial area'
    dset['ice_conc'].attrs['units'] = 'fraction m2_ice/m2_cell'

  dset["time"].attrs = {
       "long_name": "time"
  }

  # Add global attributes:
  dset.attrs['title']       = f'{varnm} on mesh025 grid from RTOFSv{rtofs_vers} CICE4 initial fields {init_date}' 
  dset.attrs['institution'] = 'NOAA NWS OMD'
  dset.attrs['source']      = 'interp_RTOFS_CICE4_ice_hsnow_mesh025.py'
  dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
  dset.attrs['region']      = 'global'

  print(f'Dumping interpolated {field_name} --> {dfliceout}\n')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')

  return


# Input RTOFS fields:
pthfld = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/cice4_fcast/{init_date}"
if fhr == 0:
  flname = "rtofs_glo.t00z.n00.cice_inst.seam_crct.nc"
else:
  flname = f"rtofs_glo.t00z.f{fhr}.cice_inst.seam_crct.nc"
dflname = os.path.join(pthfld, flname)

# Output file name:
pthintrp = os.path.join(pthfld,'interp_mesh025')
os.makedirs(pthintrp, exist_ok=True)
fliceout = f"{field_name}_RTOFSv{rtofs_vers}_{init_date}_{jdm}x{idm}.nc"
dfliceout = os.path.join(pthintrp, fliceout)

print(f"Reading {varnm} from {dflname}")

with xarray.open_dataset(dflname) as ds_ices:
  A2d = ds_ices[varnm].values.squeeze()
  units = ds_ices[varnm].attrs.get("units", None)

if units == 'cm' or units == 'centimeters':
  cff2m = 100.
elif units == 'm' or units == 'meters':
  cff2m = 1.
elif units == '1':
  cff2m = 1.
elif units is None:
  print(f"WARN: units not found for {varnm}, assumed meters")
  cff2m = 1.

AA = np.where(A2d > 1.e30, np.nan, A2d) * cff2m  # cm --> m if needed 

# Inpterolation
AAi = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
A2d = AAi.copy()


BDmask = HH >= 0
A2d[np.isnan(A2d)] = 0.
A2d[BDmask] = np.nan

time_out = np.array([np.datetime64(f"{YR:04d}-{MM:02d}-15", "ns")])
A2d = np.expand_dims(A2d, axis=0) 
write_nc(A2d, time_out, field_name, dfliceout, varnm, rtofs_vers)

  
