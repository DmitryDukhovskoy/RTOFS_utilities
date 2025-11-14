"""
  Interpolate montly clim 
  Derived monthly climatology of ice thickness in Antarctica region from
  Gridded estimates of Antarctic sea ice physical properties derived from 
  CryoSat-2 Baseline-D SAR and SARIn data spanning July 2010 through August 2021. 
  Data are generated using the CryoSat-2 Waveform-Fitting method for Antarctic sea ice (CS2WFA).

  Fons, S., Kurtz, N., & Bagnardi, M. (2022). 
  Antarctic Sea Ice Thickness Estimates from CryoSat-2: 2010-2021 (0.1.1) [Data set]. 
  Zenodo. https://doi.org/10.5281/zenodo.7327711

  gmapi indices:
  get_gmapi_CryoSat_to_mesh025.py

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
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)  

regn = 'south'
yrS  = 2011
yrE  = 2020
    
syst_info = os.uname()
machine = syst_info.nodename
        
if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
else:
  print("Unknown machine:", machine)

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
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

jdim, idim = HH.shape
LMsk = np.where(HH<0, 1, 0)

# Mask out not needed latitudes:
if regn == 'south':
  LMsk = np.where(hlat > -55, 0, LMsk)
else:
  LMsk = np.where(hlat < 50, 0, LMsk)

# Get gmapi 4 NSIDC grid points for interpolation
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump  = os.path.join(pthdata,"gmapi_NSIDC")
fgmapi  = f'CryoSat_MOM6_gmapi_{idim}x{jdim}_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)
print(f'Loading gmapi --> {dfgmapi}')

dgmapi = xarray.open_dataset(dfgmapi)
IMOM = dgmapi['mom_indx'].data
JMOM = dgmapi['mom_jndx'].data
INDX = dgmapi['gmapi_i'].data
JNDX = dgmapi['gmapi_j'].data

# Monthly ice thickness, Antarctica, original grid:
pth0   = os.path.join(pthdata, 'CryoSat2_antarctic_ice_snow_thkn')
pthice = os.path.join(pth0,'clim')
fhice  = f'CryoSat_ithkn_mnthly_clim_316x332_{regn}.nc'
dfhice = os.path.join(pthice, fhice)

print(f'Reading ice thickn climatology {dfhice}')
with xarray.open_dataset(dfhice) as ds_hice:
  LON = ds_hice['lon'].data
  LAT = ds_hice['lat'].data
  HICE = ds_hice['ice_thickness'].data

icc = 0
A3d = np.zeros((12,jdim,idim))
for imonth in range(1,13):
  print(f"Processing month {imonth} ...")

  AA = HICE[imonth-1,:].squeeze()
  CIint = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat, land_mask=True)
  #CIint = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
  #CIint = np.where(np.isnan(CIint), 0., CIint)
  CIint = np.where(HH>=0, np.nan, CIint)
  A3d[imonth-1,:,:] = CIint

A3d = A3d.astype('float32')
LON = LON.astype('float32')
LAT = LAT.astype('float32')
JD = np.arange(jdim, dtype='int32')
ID = np.arange(idim, dtype='int32')
time_months = np.arange(1,13, dtype='int32')

darr_hi = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
                   coords={"time": time_months,\
                           "jdim": JD,\
                           "idim": ID,})
darr_lon = xarray.DataArray(hlon, dims=("jdim","idim"),
                 coords={"jdim": JD,\
                         "idim": ID,})
darr_lat = xarray.DataArray(hlat, dims=("jdim","idim"),
                 coords={"jdim": JD,\
                         "idim": ID,})

dset_hi = xarray.Dataset({
  "ice_thkn": darr_hi,
  "lon": darr_lon,
  "lat": darr_lat,
})
dset_hi['time'].attrs.update({
  "long_name": "months"
})
dset_hi['ice_thkn'].attrs.update({
  "long_name": "sea ice thickness",
  "units": "m",
})
dset_hi['lon'].attrs.update({
  "long_name": "Longitudes",
  "units": "degrees_east",
})
dset_hi['lat'].attrs.update({
  "long_name": "Latitudes",
  "units": "degrees_north",
})

dset_hi.attrs.update({
  "title": f"Antarctic Sea Ice Thickness Estimates from CryoSat-2 monthly clim interpolated to mesh025 grid",
  "info": f"Monthly climatology fields on native grid years: {yrS}-{yrE}",
  "info2": "https://zenodo.org/records/7327711",
  "institution": "NOAA NWS NCEP MDC",
  "source": "interp_CryoSat_ithkn_antarct_mesh025.py",
  "contact": "dmitry.dukhovskoy@noaa.gov",
  "region": regn,
})

fliceout = f'CryoSat_hice_mnthclim_{yrS}_{yrE}_mesh025_{idim}x{jdim}_{regn}.nc'
dfliceout = os.path.join(pthice,fliceout)
print(f'Dumping interpolated ice thickness --> {dfliceout}')
dset_hi.to_netcdf(dfliceout,
        encoding={var: {'_FillValue': 1e30} for var in dset_hi.data_vars},
        format='NETCDF4')




