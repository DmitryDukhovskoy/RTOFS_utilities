"""
  Interpolate CryoSat ice thickn. monthly fileds 

  NSIDC Arctic snow depth and sea ice thickness from ICESat-2 and CryoSat-2
  CryoSat-2 monthly mean sea ice thickness for 2011 to 2013 on EASE100 grid
  12 months data

  gmapi indices: get_gmapi_CryoSat_arcticNSIDC_EASE100_to_mesh025.py
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

regn = 'north'
field_name = 'ithkn' 
tmpf = 1  # for climatology, save temporary monthly and start from last saved

parser = argparse.ArgumentParser()
parser.add_argument("--tmpf", choices=[0,1],
                    help=f"1: Save, start from last processed field, default={tmpf}", type=int)
args = parser.parse_args()

tmpf = args.tmpf if args.tmpf is not None else tmpf

save_tmp = tmpf == 1

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
LMsk = np.where(HH<0, 1, 0)

# Mask out not needed latitudes:
LMsk = np.where(hlat < 50, 0, LMsk)

print(f"NSIDC EASE100 monthly clim {field_name}  interpolating")
A3d = np.zeros((12,jdm,idm))

# Get gmapi 4 NSIDC grid points for interpolation
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump  = os.path.join(pthdata,"gmapi_NSIDC")
fgmapi  = f'CryoSat_NSIDC_EASE100_MOM6_gmapi_1440x1080_north.nc'
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

# Save temporary clim fields 
if save_tmp:
  tmp_dir = os.path.join(pthdata,'CryoSat_AWI_arctic_ithkn','tmp')
  os.makedirs(tmp_dir, exist_ok=True)

# Original data:
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthice  = os.path.join(pthdata, 'CryoSat_NSIDC_arctic_ithkn')
flice   = 'cryosat_seaice_thickness_mean_month_2011to2013.nc'
dflice  = os.path.join(pthice,flice)

assert os.path.isfile(dflice), f"Does not exist: {dflice}"
    
varnm = 'thick'

for imo in range(12):
  MM = imo + 1

  if save_tmp:
    tmp_file = os.path.join(tmp_dir, f"tmp_ithkn_NSIDC_EASE100_{MM:02d}")
    # Skip if already processed
    tmp_file_npy = f"{tmp_file}.npy"
    if os.path.exists(tmp_file_npy):
      print(f"Skipping month {MM}: already computed")
      A2d = np.load(tmp_file_npy)
      A3d[imo,:,:] = A2d
      continue

  print(f"Processing {MM:02d}")

  with xarray.open_dataset(dflice) as dsn:
    AA = dsn[varnm].isel(time=imo).data.squeeze()
    units = dsn[varnm].attrs.get('unit', None)
    if units == 'cm':
      cff_m =0.01      # cm --> m
    elif units == 'm' or units == 'kg m-3':
      cff_m = 1.
    else:
      raise Exception(f"Unrecognized units {units}")
  
    AA = np.where(AA > 1.e30, np.nan, AA) * cff_m  # cm --> m  

    HSint = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)

    # Fill the N.Pole hole:
    HSint = mrmom.fill_npole(HSint, hlon, hlat, HH, Rpole=2.)
    A2d = np.where(HH>=0, np.nan, HSint)
    A2d = np.expand_dims(A2d, axis=0) 
    A3d[imo,:,:] = A2d

    if save_tmp:
      print(f"Saving temporary --> {tmp_file}")
      np.save(tmp_file, A2d)


fliceout = f"{field_name}_CryoSat_arcticNSIDC_EASE100_mnthclim.nc"
pthnsidc = os.path.join(pthdata,'CryoSat_NSIDC_arctic_ithkn','clim')
dfliceout = os.path.join(pthnsidc,fliceout)
   
time_out = np.arange(1,13)
jdim, idim = HH.shape

darr_cice = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
                   coords={"time": time_out,\
                           "jdim": np.arange(jdm),\
                           "idim": np.arange(idm)})

dset = xarray.Dataset({"ice_thkn": darr_cice})
dset['ice_thkn'].attrs['long_name'] = 'ice thickness'
dset['ice_thkn'].attrs['units'] = 'm'

dset["time"].attrs = {
       "long_name": "time"
  }

# Add global attributes:
dset.attrs['title']       = f'Arctic ice thickness clim from 2011-2013 CryoSat2 on EASE100 grid' 
dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
dset.attrs['source']      = 'interp_CryoSat_arcticNSIDC_EASE100_ithkn_mesh025.py'
dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
dset.attrs['region']      = 'north'

print(f'Dumping interpolated {field_name} --> {dfliceout}\n')
dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')



f_chck = False
if f_chck:
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  clrmp.set_under(color=[1,1,1]) 

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.1, 0.15, 0.8, 0.8])
      
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l', ax=ax1)
  xh, yh = m(hlon, hlat)

  m.drawparallels(np.arange(60, 90, 5), labels=[1,0,0,0])
  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,1])
  m.drawcoastlines()

  AP = HSint.copy()
  AP[HH >= 0] = np.nan   # land
  AP[np.isnan(AP) & (HH < 0)] = -1.  # ocean
  img = m.pcolormesh(xh,yh, AP, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"{varnm}, CryoSat NSIDC EASE100 {MM:02d}")

  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = 'interp_CryoSat_arcticNSIDC_EASE100_ithkn_mesh025.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)





