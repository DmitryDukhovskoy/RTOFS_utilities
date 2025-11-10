"""
  Interpolate NASA AMSR &  SSM/I snow thickness in Antarctic
  to MOM6/CICE6 mesh025 grid  

  Since NASA and NOAA NSIDC seem to use same South polar proejction, 
  use gmapi indices derived for NSIDC: get_gmapi_NSIDC_to_mesh025.py

  Monthly hsnow clim data derived (on GFDL/PPAN ):
  derive_hsnow_monthly_AMSR_Antarctic.py

  daily snow data are on GFDL/PPAN:
  /work/Dmitry.Dukhovskoy/data/snow_nasa/{YR}


  NSIDC fields from 
  https://noaadata.apps.nsidc.org/NOAA/G02202_V6/north/daily/2025/

  Derive monhtly clim of snow thickness in Antarctica 
  derived from AMSR 19 and 37 GHz microwave brightness temperatures
  https://earth.gsfc.nasa.gov/cryo/data/antarctic-snow-depth-sea-ice
  available data: 1993-2008
  daily data


  Code to derived monthly SSM/I fileds are on PPAN:
  derive_hsnow_monthly_SSMI_Antarctic.py

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
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

regn = 'south'  # only south region has been done so far
yrS = 1998
yrE = 2007

parser = argparse.ArgumentParser()
parser.add_argument(f"--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument(f"--yrs", help=f"start year of derived monthly clim, default={yrS}", type=int)
parser.add_argument(f"--yre", help=f"end year of derived monthly clim, default={yrE}", type=int)
args = parser.parse_args()
  
regn = args.regn if args.regn else regn
yrS = args.yrs if args.yrs else yrS
yrE = args.yre if args.yre else yrE

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
if regn == 'south':
  LMsk = np.where(hlat > -55, 0, LMsk)
else:
  LMsk = np.where(hlat < 50, 0, LMsk)

# Get gmapi 4 NSIDC grid points for interpolation
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump  = os.path.join(pthdata,"gmapi_NSIDC")
fgmapi  = f'SSMI_hsnow_MOM6_gmapi_{idm}x{jdm}_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)
print(f'Loading gmapi --> {dfgmapi}')

dgmapi = xarray.open_dataset(dfgmapi)
IMOM = dgmapi['mom_indx'].data
JMOM = dgmapi['mom_jndx'].data
INDX = dgmapi['gmapi_i'].data
JNDX = dgmapi['gmapi_j'].data

# monthly snow depth, Antarctica:
pthsnow = os.path.join(pthdata,'snow_nasa','monthly_clim')
flsnow = f'SSMI_Antarctic_hsnow_month_clim_{yrS}_{yrE}_316x332.nc'
dflsnow = os.path.join(pthsnow, flsnow)

print(f"Loading {dflsnow}")
with xarray.open_dataset(dflsnow) as ds_snow:
  LON = ds_snow['lon'].data
  LAT = ds_snow['lat'].data
  HSNOW = ds_snow['snow_depth'].data

jdim, idim = HH.shape

icc = 0
A3d = np.zeros((12,jdm,idm))
for imonth in range(1,13):
  print(f"Processing month {imonth} ...")

  AA = HSNOW[imonth-1,:].squeeze()
  CIint = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
  CIint = np.where(HH>=0, np.nan, CIint)
  A3d[imonth-1,:,:] = CIint

A3d = A3d.astype('float32')
LON = LON.astype('float32')
LAT = LAT.astype('float32')
JD = np.arange(jdim, dtype='int32')
ID = np.arange(idim, dtype='int32')
time_months = np.arange(1,13, dtype='int32')

darr_hs = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
                   coords={"time": time_months,\
                           "jdim": JD,\
                           "idim": ID,})
darr_lon = xarray.DataArray(hlon, dims=("jdim","idim"),
                 coords={"jdim": JD,\
                         "idim": ID,})
darr_lat = xarray.DataArray(hlat, dims=("jdim","idim"),
                 coords={"jdim": JD,\
                         "idim": ID,})

dset_hs = xarray.Dataset({
  "snow_depth": darr_hs,
  "lon": darr_lon,
  "lat": darr_lat,
})
dset_hs['time'].attrs.update({
  "long_name": "months"
})
dset_hs['snow_depth'].attrs.update({
  "long_name": "snow depth on ice",
  "units": "cm",
})
dset_hs['lon'].attrs.update({
  "long_name": "Longitudes",
  "units": "degrees_east",
})
dset_hs['lat'].attrs.update({
  "long_name": "Latitudes",
  "units": "degrees_north",
})

dset_hs.attrs.update({
    "title": f"Snow depth climatology from SSM/I ({yrS}–{yrE}) interpolated onto mesh025 grid",
    "institution": "NOAA NWS NCEP MDC",
    "source": "interp_SSMI_hsnow_monthly_antarct_mesh025.py",
    "region": regn,
})

fliceout = f'SSMI_hsnow_mnthclim_{yrS}_{yrE}_mesh025_{idm}x{jdm}_{regn}.nc'
dfliceout = os.path.join(pthsnow,fliceout)
print(f'Dumping interpolated snow depth --> {dfliceout}')
dset_hs.to_netcdf(dfliceout, 
        encoding={var: {'_FillValue': 1e30} for var in dset_hs.data_vars},
        format='NETCDF4')


