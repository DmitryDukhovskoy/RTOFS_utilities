"""
  Use  interpolated ICESat2 ice thickn. monthly fileds 
  summer 05-08 months only

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
fld_name = 'ithkn' 

parser = argparse.ArgumentParser()
parser.add_argument("--field", help=f"Field to interpolate, default={fld_name}",
                    choices=['hsnow','ithkn'], type=str)
args = parser.parse_args()
  
fld_name = args.field if args.field else fld_name

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
    

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
jdim, idim = hlon.shape


if fld_name == 'ithkn':
  varnm = 'ice_thkn'
elif fld_name == 'hsnow':
  varnm = 'snow_depth'

# Available months: May - August
A3d = np.zeros((4,jdim,idim))
time_out = []
imo = -1
for MM in range(5,9):

  ASUM = np.zeros((jdim,idim))
  count_ice = np.zeros((jdim,idim))
  for YR in range(YRS,YRE+1):
    print(f"Processing {YR}/{MM}")
    if YR == 2021 and MM > 7:
      print(f"Last record: 2021/07 skipping ...")
      continue

    pthintrp = os.path.join(pthdata,'ICESat2_arctic_summer_ithkn_hsnow','interp_mesh025')
    fliceout = f"{fld_name}_ICESat2_arctic_{YR}{MM:02d}_{jdim}x{idim}.nc"
    dflice = os.path.join(pthintrp, fliceout)

    print(f"Reading {dflice}")

    with xarray.open_dataset(dflice) as ds_ices:
      A2d = ds_ices[varnm].values.squeeze()
      units = ds_ices[varnm].attrs.get("units", None)


    if units == 'cm' or units == 'centimeters':
      A2d = A2d * 0.01   # cm ---> m

    # Count & average only non-zero thicknesses
    # Do not count 0 thickn. this will bias ice thickness 
    ice_grid = np.isfinite(A2d) & (A2d > 0.0)
    ASUM[ice_grid] += A2d[ice_grid]
    count_ice[ice_grid] += 1

  # Average:
  Aavrg = np.divide(ASUM, count_ice, out=np.zeros_like(ASUM), where = count_ice > 0)
  Aavrg[np.isnan(A2d)] = np.nan 

  imo += 1
  A3d[imo,:,:] = Aavrg
  time_out.append(MM) 

# Dump netcdf:
darr_cice = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
                   coords={"time": time_out,\
                           "jdim": np.arange(jdim),\
                           "idim": np.arange(idim)})
if fld_name == 'hsnow':
  dset = xarray.Dataset({"snow_depth": darr_cice})
  dset['snow_depth'].attrs['long_name'] = 'snow depth on ice'
  dset['snow_depth'].attrs['units'] = 'm'
elif fld_name == 'ithkn':
  dset = xarray.Dataset({"ice_thkn": darr_cice})
  dset['ice_thkn'].attrs['long_name'] = 'ice thickness'
  dset['ice_thkn'].attrs['units'] = 'm'

dset["time"].attrs = {
     "long_name": "time",
     "units": "Month",
}

# Add global attributes:
dset.attrs['title']       = f'{varnm} on mesh025 grid from monthly summer Arctic data ICESat-2, v2 {YRS}-{YRE}' 
dset.attrs['institution'] = 'NOAA NWS NCEP EMC'
dset.attrs['source']      = 'derive_ithkn_hsnow_summer_clim_ICESat2_arctic.py'
dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
dset.attrs['region']      = 'north'

pthclm = os.path.join(pthdata,'ICESat2_arctic_summer_ithkn_hsnow','clim')
flclm = f"{fld_name}_ICESat2_arctic_mnth_{YRS}-{YRE}_{jdim}x{idim}.nc"
dflclm = os.path.join(pthclm, flclm)

print(f'Dumping interpolated {fld_name} --> {dflclm}\n')
dset.to_netcdf(dflclm, format='NETCDF4', engine='netcdf4')




f_chck = False
if f_chck:
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  clrmp.set_under(color=[1,1,1]) 

  MP = 8
  TM = np.array(time_out)
  Aclm = A3d[TM == MP,:,:].squeeze()

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.1, 0.15, 0.8, 0.8])
      
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l', ax=ax1)
  xh, yh = m(hlon, hlat)

  m.drawparallels(np.arange(60, 90, 5), labels=[1,0,0,0])
  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,1])
  m.drawcoastlines()

  img = m.pcolormesh(xh,yh, Aclm, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"{varnm}, ICESat-2 mnthly clim {MP}")

  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = 'derive_ithkn_hsnow_summer_clim_ICESat2_arctic.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)





