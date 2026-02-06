"""
  Interpolate CryoSat hsnow or ice thickn. monthly fileds 
  winter months only

  When averaging, 0-thickness is ignored, i.e. only hice > 0 is considered for 
  computing the multi-year mean to track mean ice thickness in the grid cell 
  only when it presents there

  AWI L4 gridded 25 km 
  Arctic snow depth, density, and sea ice thickness, freeboard etc 
   from CryoSat-2
  https://data.seaiceportal.de/relaunch/thickness.php?lang=en

  Data are on polar stereographic coordinates

  gmapi indices: get_gmapi_CryoSat_arcticAWI_to_mesh025.py
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
YRS = 2015
YRE = 2024
field_name = 'ithkn' 
avrg = 1  # 1 - derive climatology and interp, 0 - interpolate monthly data by years and save
tmpf = 1  # for climatology, save temporary monthly and start from last saved

parser = argparse.ArgumentParser()
parser.add_argument("--field", help=f"Field to interpolate, default={field_name}",
                    choices=['hsnow','rhosn','ithkn','iconc','rhoice'], type=str)
parser.add_argument("--yrs", help=f"Start year of data, default={YRS}", type=int)
parser.add_argument("--yre", help=f"End year of data, default={YRE}", type=int)
parser.add_argument("--avrg", help=f"1 = average to produce climatology, default={avrg}",
                   choices=[0,1], type=int)
parser.add_argument("--tmpf", choices=[0,1],
                    help=f"1: Save, start from last processed field, default={tmpf}", type=int)
args = parser.parse_args()
  
field_name = args.field if args.field else field_name
YRS = args.yrs if args.yrs else YRS
YRE = args.yre if args.yre else YRE
avrg = args.avrg if args.avrg else avrg
tmpf = args.tmpf if args.tmpf is not None else tmpf

mnthly_clim = avrg == 1  # True - produce monthly climatologies by averaging over years YRS-YRE and intrp
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

if mnthly_clim:
  print(f"Deriving monthly clim {field_name} for {YRS}-{YRE} and interpolating")
  A3d = np.zeros((12,jdm,idm))
else:
  print(f"Interpolating monthly {field_name} for {YRS}-{YRE}")
  A3d = np.zeros((jdm,idm))

# Get gmapi 4 NSIDC grid points for interpolation
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump  = os.path.join(pthdata,"gmapi_NSIDC")
fgmapi  = 'CryoSat_AWI_EASE2_MOM6_gmapi_1440x1080_north.nc'
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
pthice  = os.path.join(pthdata, 'CryoSat_AWI_arctic_ithkn')
fname0 = 'awi-siral-l3c-sithick-cryosat2-rep-nh_25km_ease2'
fsfx   = 'fv2p6'       

match field_name:
  case 'ithkn':
    varnm = 'sea_ice_thickness'
  case 'iconc':
    varnm = 'sea_ice_concentration'
  case 'rhoice':
    varnm = 'sea_ice_density'
  case 'rhosn':
    varnm = 'snow_density'
  case 'hsnow':
    varnm = 'snow_depth'

def write_nc(A3d, time_out, field_name, dfliceout, attr_str ):
  darr_cice = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
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
  elif field_name == 'rhoice':
    dset = xarray.Dataset({"ice_dens": darr_cice})
    dset['ice_dens'].attrs['long_name'] = 'ice density'
    dset['ice_dens'].attrs['units'] = 'kg / m^3'
  elif field_name == 'rhosn':
    dset = xarray.Dataset({"snow_dens": darr_cice})
    dset['snow_dens'].attrs['long_name'] = 'snow density'
    dset['snow_dens'].attrs['units'] = 'kg / m^3'

  dset["time"].attrs = {
       "long_name": "time"
  }

  # Add global attributes:
  dset.attrs['title']       = f'Arctic {attr_str} from CryoSat2 AWI L4 product' 
  dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
  dset.attrs['source']      = 'interp_CryoSat_arcticAWI_iceflds_mesh025.py'
  dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
  dset.attrs['region']      = 'north'

  print(f'Dumping interpolated {field_name} --> {dfliceout}\n')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')

  return

mnth_keep = []   # some months are missing for ice thickness
for imo in range(12):
  MM = imo + 1

  if MM > 4 and MM < 10:
    # only winter months
    continue

  mnth_keep.append(MM)
  if save_tmp and mnthly_clim:
    tmp_file = os.path.join(tmp_dir, f"tmp_{field_name}_AWIinterp_{MM:02d}")
    # Skip if already processed
    tmp_file_npy = f"{tmp_file}.npy"
    if os.path.exists(tmp_file_npy):
      print(f"Skipping month {MM}: already computed")
      A2d = np.load(tmp_file_npy)
      A3d[imo,:,:] = A2d
      continue

  nyrs = 0
  Asum = np.full_like(LON, 0.)
  count_ice = np.full_like(LON,0).astype(int)
  for YR in range(YRS,YRE+1):
    print(f"Processing {YR}/{MM:02d}")
    flice   = f"{fname0}-{YR}{MM:02d}-{fsfx}.nc"
    dflice  = os.path.join(pthice,flice)
    
    assert os.path.isfile(dflice), f"Does not exist: {dflice}"

    with xarray.open_dataset(dflice) as dsn:
      AA = dsn[varnm].data.squeeze()
      units = dsn[varnm].attrs.get('units', None)
      if units == 'cm':
        cff_m =0.01      # cm --> m
      elif units == 'm' or units == 'kg m-3':
        cff_m = 1.
      else:
        raise Exception(f"Unrecognized units {units}")
  
    AA = np.where(AA > 1.e30, np.nan, AA) * cff_m  # cm --> m  
    BDmask = np.isnan(AA)  # keep nan mask
    nyrs += 1
 
    if not mnthly_clim:  
      HSint = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)

      # Fill the N.Pole hole:
      HSint = mrmom.fill_npole(HSint, hlon, hlat, HH, Rpole=2.)
      A2d = np.where(HH>=0, np.nan, HSint)
      A2d = np.expand_dims(A2d, axis=0) 

      fliceout = f"{field_name}_CryoSat_arcticAWI_{YR}{MM:02d}_{jdm}x{idm}.nc"
      pthnsidc = os.path.join(pthdata,'CryoSat_AWI_arctic_ithkn','interp_monthly')
      dfliceout = os.path.join(pthnsidc,fliceout)
      time_out = np.array([np.datetime64(f"{YR:04d}-{MM:02d}-01", "ns")])
     
      write_nc(A2d, time_out, field_name, dfliceout, varnm)

    else:
      #AA[HH >= 0] = np.nan
      #Asum += AA  # simple averaging will mix no ice and ice in the grid cell, not good
      ice_grid = np.isfinite(AA) & (AA > 0.0)
      Asum[ice_grid] += AA[ice_grid]
      count_ice[ice_grid] += 1

  if mnthly_clim:
    #Asum /= nyrs
    Asum = np.divide(Asum, count_ice, out=np.zeros_like(Asum), where = count_ice > 0)
    # Mask Land and N. Pole hole:
    Asum[BDmask] = np.nan

    HSint = msisrlx.interp2Dfld(Asum, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
    # Fill the N.Pole hole:
    HSint = mrmom.fill_npole(HSint, hlon, hlat, HH, Rpole=2.5)
    A2d = np.where(HH>=0, np.nan, HSint)
    A2d = np.expand_dims(A2d, axis=0)

    A3d[imo,:,:] = A2d

    if save_tmp:
      print(f"Saving temporary --> {tmp_file}")
      np.save(tmp_file, A2d)

if mnthly_clim:
  fliceout = f"{field_name}_CryoSat_arcticAWI_mnthclim_{YRS}-{YRE}_{jdm}x{idm}.nc"  
  pthnsidc = os.path.join(pthdata,'CryoSat_AWI_arctic_ithkn','clim')
  dfliceout = os.path.join(pthnsidc,fliceout)
  time_out = np.array(mnth_keep)

  idx_keep = time_out - 1
  A3d_keep = A3d[idx_keep, :, :]
  write_nc(A3d_keep, time_out, field_name, dfliceout, varnm)


f_chck = False
if f_chck:
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  clrmp.set_under(color=[1,1,1]) 

  AP = HSint.copy()
  AP[HH >= 0] = np.nan   # land
  AP[np.isnan(AP) & (HH < 0)] = -1.  # ocean

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

  img = m.pcolormesh(xh,yh, AP, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"{varnm}, CryoSat AWI, {YR}/{MM:02d}")

  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = 'interp_CryoSat_arcticAWI_iceflds_mesh025.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)





