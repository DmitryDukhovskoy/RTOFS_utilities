"""
  Interpolate CryoSat hsnow or ice thickn. monthly fileds 
  winter months only

  When averaging, 0-thickness is ignored, i.e. only hice > 0 is considered for 
  computing the multi-year mean to track mean ice thickness in the grid cell 
  only when it presents there


  Summer ice and snow depth, densities, Arctic

  https://zenodo.org/records/18004849
  
  Monthly gridded summer Arctic sea ice thickness from ICESat-2, v1
  Creators
  Petty, Alek Aaron (Producer)1, 2
  ORCID icon
  Description
  Monthly gridded summer Arctic sea ice thickness from ICESat-2. 
  Produced by combining Release 006 ATL10 freeboards with SnowModel-LG snow loading, 
  processed as in the IS2SITMOGR4 winter Arctic thickness dataset (https://nsidc.org/data/IS2SITMOGR4). 
    
ICESat-2 (Ice, Cloud, and land Elevation Satellite-2) uses a space-based laser altimeter to measure surface elevation changes—especially ice sheets and sea ice—down to centimeter accuracy.

  
  get_gmapi_ICESat2_arctic_to_mesh025.py

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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_regmom as mrmom 
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

regn = 'north'
YRS = 2019
YRE = 2021
field_name = 'ithkn' 

parser = argparse.ArgumentParser()
parser.add_argument("--field", help=f"Field to interpolate, default={field_name}",
                    choices=['hsnow','rhosn','ithkn','iconc'], type=str)
args = parser.parse_args()
  
field_name = args.field if args.field else field_name

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

print(f"Interpolating monthly {field_name} for {YRS}-{YRE}")
A3d = np.zeros((jdm,idm))

# Get gmapi 4 NSIDC grid points for interpolation
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump  = os.path.join(pthdata,"gmapi_NSIDC")
fgmapi  = f'ICESat2_IS2SITMOGR4_MOM6_gmapi_{idm}x{jdm}_{regn}.nc'
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


# Original data:
# There are many versions of ice thickness and snow depth estimates in the file
# ice_thickness_sm_e5_int: Monthly mean gridded and smoothed/interpolated sea ice thickness 
#                          calculated using redistributed SnowModel-LG snow loading 
#                          with ERA5 forcing (Liston et al., 2021, 10.5067/27A0P5M6LZBI, SM) 
#                          and fixed ice density (916 kg/m3)
#
# snow_depth_sm_e5_int:   Monthly mean gridded and smoothed/interpolated redistributed 
#                         SnowModel-LG with ERA5 forcing snow depths
#
# snow_density_sm_e5:     Monthly mean gridded SnowModel-LG with ERA5 forcing (Liston et al., 2021, 
#                         snow density. Data currently available up to July 2021 on the NSIDC
#

if field_name == 'ithkn':
  varnm = 'ice_thickness_sm_e5_int'
elif field_name == 'hsnow':
  varnm = 'snow_depth_sm_e5_int'
elif field_name == 'rhosn':
  varnm = 'snow_density_sm_e5'
elif field_name == 'iconc':
  varnm = 'sea_ice_conc'
else:
  raise Exception(f"Unrecognized variable {field_name}")

def write_nc(A2d, time_out, field_name, dfliceout, varnm):
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
  elif field_name == 'rhosn':
    dset = xarray.Dataset({"snow_dens": darr_cice})
    dset['snow_dens'].attrs['long_name'] = 'snow density'
    dset['snow_dens'].attrs['units'] = 'kg / m^3'

  dset["time"].attrs = {
       "long_name": "time"
  }

  # Add global attributes:
  dset.attrs['title']       = f'{varnm} on mesh025 grid from monthly gridded summer Arctic sea ice thickness from ICESat-2, v2' 
  dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
  dset.attrs['source']      = 'interp_ICESat2_icesnow_summer_arctic_mesh025.py'
  dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
  dset.attrs['region']      = 'north'

  print(f'Dumping interpolated {field_name} --> {dfliceout}\n')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')

  return


for YR in range(YRS,YRE+1):
  for MM in range(5,9):
    if YR == 2021 and MM > 7:
      print(f"Last record: 2021/07 skipping ...")
      continue

    # Output file name:
    pthintrp, _  = mc6util.pathfname_icesnow_mesh025(fyaml, node_nm, 'ICESat2_mnth_interp')
    fliceout = f"{field_name}_ICESat2_arctic_{YR}{MM:02d}_{jdm}x{idm}.nc"
    dfliceout = os.path.join(pthintrp, fliceout)
    if os.path.isfile(dfliceout):
      print(f"  ===> Already created {fliceout}, skipping ...")
      continue

    pthfld, flname = mc6util.pathfname_icesnow_mesh025(fyaml, node_nm, 'ICESat2_orig', YR=YR, MM=MM)
    dflname = os.path.join(pthfld, flname)
    print(f"Reading {varnm} from {dflname}")

    with xarray.open_dataset(dflname) as ds_ices:
      A2d = ds_ices[varnm].values.squeeze()
      units = ds_ices[varnm].attrs.get("units", None)

    if units == 'cm' or units == 'centimeters':
      cff2m = 100.
    elif units == 'm' or units == 'meters':
      cff2m = 1.

    AA = np.where(A2d > 1.e30, np.nan, A2d) * cff2m  # cm --> m if needed 
 
    # Inpterolation
    AAi = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)

    # Fill the N.Pole hole:
    #HSint = mrmom.fill_npole(AAi, hlon, hlat, HH, Rpole=2.)
    #A2d = HSint.copy()
    A2d = AAi.copy()

    BDmask = HH >= 0
    A2d[np.isnan(A2d)] = 0.
    A2d[BDmask] = np.nan

    time_out = np.array([np.datetime64(f"{YR:04d}-{MM:02d}-15", "ns")])
    A2d = np.expand_dims(A2d, axis=0) 
    write_nc(A2d, time_out, field_name, dfliceout, varnm)

  
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
  xL, yL = m(LON, LAT)

  m.drawparallels(np.arange(60, 90, 5), labels=[1,0,0,0])
  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,1])
  m.drawcoastlines()

  img = m.pcolormesh(xh,yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"{varnm}, CryoSat AWI, {YR}/{MM:02d}")

  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = 'interp_CryoSat_arcticAWI_iceflds_mesh025.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)





