"""
  Interpolate NSIDC sea ice concentration 
  to MOM6/CICE6 mesh025 grid

  gmapi indices: get_gmapi_NSIDC_to_mesh025.py

  NSIDC fields from 
  https://noaadata.apps.nsidc.org/NOAA/G02202_V6/north/daily/2025/

  see scripts/./get_NRT_seaconc.sh --yrs 2025 --ms 9 --regn south
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib  
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
import mod_misc1 as mmisc
import mod_mom6 as mmom6
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

YR = 2025
MM = 1 
DD = 15
tmpf = 1
ddS = 1    # month day to start processing
ddE = 31   # end processing
fsave = 1

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", type=str, required=True)
parser.add_argument("--yr", help=f"year of NSIDC data, default={YR}", type=int)
parser.add_argument("--mm", help=f"month of NSIDCS data to interpolate, default={MM}", type=int)
parser.add_argument("--tmpf", choices=[0,1], 
                    help=f"1: Save, start from last processed field, default={tmpf}", type=int)
parser.add_argument("--fsave", help=f"Save final dataset with all days as netcdf, default={fsave}", 
                    choices=[0,1], type=int)
parser.add_argument("--ddS", help=f"Day to start interpolation, default={ddS}", type=int)
parser.add_argument("--ddE", help=f"Day to end interpolation, default={ddE}", type=int) 
args = parser.parse_args()
  
regn = args.regn if args.regn else None
YR   = args.yr if args.yr else YR
MM   = args.mm if args.mm else MM
tmpf = args.tmpf if args.tmpf is not None else tmpf
fsave = args.fsave if args.fsave is not None else fsave
ddS = args.ddS if args.ddS is not None else ddS 
ddE = args.ddE if args.ddE is not None else ddE 

save_nc = fsave == 1
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
fgmapi  = f'NSIDC_NRTice_MOM6_gmapi_{idm}x{jdm}_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)
print(f'Loading gmapi --> {dfgmapi}')

with xarray.open_dataset(dfgmapi) as dgmapi:
  IMOM = dgmapi['mom_indx'].data
  JMOM = dgmapi['mom_jndx'].data
  INDX = dgmapi['gmapi_i'].data
  JNDX = dgmapi['gmapi_j'].data

def read_NSIDC(YR,MM,DD,regn,pthnsidc,varnm):
  if regn == 'south':
    if YR >= 2025:
      flnsidc = f"sic_pss25_{YR}{MM:02d}{DD:02d}_am2_v06r00.nc"  
    else:
      flnsidc = f"sic_pss25_{YR}{MM:02d}{DD:02d}_F17_v06r00.nc"  
  else:
    if YR >= 2025:
      flnsidc = f"sic_psn25_{YR}{MM:02d}{DD:02d}_am2_v06r00.nc"  
    else:
      flnsidc = f"sic_psn25_{YR}{MM:02d}{DD:02d}_F17_v06r00.nc"  

  with xarray.open_dataset(os.path.join(pthnsidc,flnsidc)) as ds_nsidc:
    A = ds_nsidc[varnm].data.squeeze()

  return A

# Get lon/lat for NSIDC data
pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")
Xnsidc = read_NSIDC(YR,MM,1,regn,pthnsidc,'x')
Ynsidc = read_NSIDC(YR,MM,1,regn,pthnsidc,'y')

XX, YY = np.meshgrid(Xnsidc, Ynsidc, indexing='xy')
# Determine ellipsoid parameters from NSIDC information
# Note that Radius of ellipsoid WGS84 is typically referred to major semi-axis (equatorial radius)
if regn == 'south':
  slat0 = -70.  # latitude of 0 distortion, standard lat. 
  lon0  = -90.   # orientation of the 0 longitude wth to X axis on polar grid, not NSIDC is fliped upside down
  flat_inv = 298.279411123064
  ax_maj = 6378273.
  flat = 1./flat_inv  # flattening
  eccentr = np.sqrt(2*flat - flat**2)
  R_polar = ax_maj*(1.-flat)   # polar radius or semi-minor axis


  LON, LAT = mmisc.convert_polarXY_lonlat(XX,YY, North=False, E=eccentr, RE=ax_maj, SLAT=slat0, LON0_dir=lon0)
  assert np.max(LAT) < 0., f"For southern hemisphere latitudes should be < 0"
  LON = -LON   # nor sure why but this makes sign of the longitudes right

elif regn == 'north':
  with xarray.open_dataset(dfgmapi) as dgmapi: 
    LON = dgmapi['longit'].data
    LAT = dgmapi['latit'].data

# Save temporary fields 
if save_tmp:
  tmp_dir = os.path.join(pthdata,'NRT_NOAA_NSIDC_seaconc','tmp')
  os.makedirs(tmp_dir, exist_ok=True)

icc = 0
ndays = mtime.month_days(MM,YR)
ddE = np.min([ddE, ndays])
A3d = np.zeros((ndays,jdm,idm))
time_days = np.arange(ddS,ddE+1)
print(f"Start processing data for time range: {ddS} - {ddE}")
print(f"Saving netcdf at the end: {save_nc}")

for mday in range(ddS,ddE+1):
  print(f"Processing {YR}/{MM}/{mday} ...")

  if save_tmp:
    tmp_file = os.path.join(tmp_dir, f"tmp_NSIDCinterp_{YR}{MM:02d}{mday:02d}_{regn}.npy")
    # Skip if already processed
    #tmp_file_npy = f"{tmp_file}.npy"
    if os.path.exists(tmp_file):
      print(f"Skipping {YR}/{MM}/{mday}: already computed <--- {tmp_file}")
      CIint = np.load(tmp_file)
      A3d[mday-1,:,:] = CIint
      continue
      

  AA = read_NSIDC(YR, MM, mday, regn, pthnsidc, 'cdr_seaice_conc')
  CIint = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
  CIint = np.where(HH>=0, np.nan, CIint)
  A3d[mday-1,:,:] = CIint

  if save_tmp:
    print(f"Saving tmp file --> {tmp_file}")
    np.save(tmp_file, CIint)

if not save_nc:
  print(f"Final netcdf is not saved, save_nc={save_nc}")
 
else:
  darr_cice = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
                     coords={"time": time_days,\
                             "jdim": np.arange(jdm),\
                             "idim": np.arange(idm)})
  dset = xarray.Dataset({"ice_conc": darr_cice})
  dset['ice_conc'].attrs['long_name']='ice partial area'
  # Add global attributes:
  dset.attrs['title']       = 'NRT NSIDC v6 sea ice conc daily interpolated onto mesh025 grid'
  dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
  dset.attrs['source']      = 'interp_NSIDC_mesh025.py'
  dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
  dset.attrs['region']      = regn
  dset.attrs['Grid_idm_jdm'] = f'{idm}x{jdm}'

  fliceout = f'NSIDC_iconc_interp_mesh025_{jdm}x{idm}_{YR}{MM:02d}_{regn}.nc'
  dfliceout = os.path.join(pthnsidc,fliceout)
  print(f'Dumping interpolated ice conc --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')

