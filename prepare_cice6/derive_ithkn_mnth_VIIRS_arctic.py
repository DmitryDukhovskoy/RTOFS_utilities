"""
  Derive monthly clim  VIIRS ice thickness on mesh025
  To fill gaps in summer months

  see: interp_VIIRS_arctic_ithkn_mesh025.py 

  ice thickness Arctic
  daily fields, 750-m resolution

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
import mod_regmom as mrmom 
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

regn  = 'north'
YRS   = 2021
YRE   = 2025
MM    = 9

parser = argparse.ArgumentParser()
parser.add_argument("--mm", help=f"Month of data, default={MM}", type=int)
args = parser.parse_args()
  
MM    = args.mm    if args.mm   else MM

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


pthviirs = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/VIIRS_ithkn_highres/interp_daily'
for YR in range(YRS,YRE+1):
  ndays = mtime.month_days(MM, YR)
  for DD in range(1,ndays+1):
    dnmb = mtime.datenum(YR,MM,DD)
    flviirs = f"ithkn_VIIRS_arctic_{YR}{MM:02d}{DD:02d}_1080x1440.nc"



# Mask out not needed latitudes:
LMsk = np.where(hlat < 50, 0, LMsk)

# Get gmapi 4 NSIDC grid points for interpolation
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump  = os.path.join(pthdata,"gmapi_NSIDC")
fgmapi  = f'VIIRS_ithkn_MOM6_gmapi_{idm}x{jdm}_{regn}.nc'
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


for dnmb in DAYS_LIST:
  YR, MM, DD = mtime.datevec(dnmb)[:3] 
  print(f"Interpolating VIIRS daily ithkn {YR}/{MM}/{DD}")
  A3d = np.zeros((jdm,idm))

  # Original data:
  pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
  pthice  = os.path.join(pthdata, 'VIIRS_ithkn_highres')
  jday1 = jday - 3
  jday2 = jday
  flice = f'VXSACW_B{YR}{jday1:03d}_B{YR}{jday2:03d}_H4_NP06_edgemask_IceThickness.nc' 
  varnm = 'IceThickness'

  if read_url:
    url = f"https://www.star.nesdis.noaa.gov/thredds/dodsC/IceThickVIIRSnppSectorFourDayNP06/{YR}/"
    flviirs = f"VXSACW_B{YR}{jday-3:03d}_B{YR}{jday:03d}_H4_NP06_edgemask_IceThickness.nc"
    dflice = os.path.join(url, flviirs)
  else:
    dflice  = os.path.join(pthice,flice)
    assert os.path.isfile(dflice), f"Does not exist: {dflice}"

  print(f"Processing {YR}/{MM:02d}/{DD:02d} <-- {dflice}")


  with xarray.open_dataset(dflice) as dsn:
    AA = dsn[varnm].data.squeeze()
    units = dsn[varnm].attrs.get('units', None)
    if units == 'cm':
      cff_m =0.01      # cm --> m
    elif units == 'm' or units == 'meter':
      cff_m = 1.
    else:
      raise Exception(f"Unrecognized units {units}")

  AA = np.where(AA > 1.e30, np.nan, AA) * cff_m  # cm --> m  

  HSint = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat, info_step=5000)

  # Fill the N.Pole hole: - no need for VIIRS
  #HSint = mrmom.fill_npole(HSint, hlon, hlat, HH, Rpole=2.)
  A2d = np.where(HH>=0, np.nan, HSint)
  A2d = np.expand_dims(A2d, axis=0)

  fliceout = f"ithkn_VIIRS_arctic_{YR}{MM:02d}{DD:02d}_{jdm}x{idm}.nc"
  pthnsidc = os.path.join(pthdata,'VIIRS_ithkn_highres','interp_daily')
  dfliceout = os.path.join(pthnsidc,fliceout)
  os.makedirs(pthnsidc, exist_ok=True)

  time_out = np.array([np.datetime64(f"{YR:04d}-{MM:02d}-{DD:02d}", "ns")])

  darr_cice = xarray.DataArray(A2d, dims=("time","jdim","idim"),\
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
  dset.attrs['title']       = f'Arctic ice thickness from VIIRS hig-res 4-day composite L3 data' 
  dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
  dset.attrs['source']      = 'interp_VIIRS_arctic_ithkn_mesh025.py'
  dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
  dset.attrs['region']      = 'north'

  print(f'Dumping interpolated ithkn --> {dfliceout}\n')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')



f_chck = False
if f_chck:
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  clrmp.set_under(color=[1,1,1]) 

  AP = A2d.squeeze()
  AP[HH >= 0] = np.nan   # land
  AP[np.isnan(AP) & (HH < 0)] = -1.  # ocean

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.1, 0.15, 0.8, 0.8])
      
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l', ax=ax1)
  xh, yh = m(hlon, hlat)
  #xL, yL = m(LON, LAT)

  m.drawparallels(np.arange(60, 90, 5), labels=[1,0,0,0])
  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,1])
  m.drawcoastlines()

  img = m.pcolormesh(xh,yh, AP, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"{varnm}, VIIRS 4-day composite, {YR}/{MM:02d}/{DD:02d}")

  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = 'interp_VIIRS_arctic_ithkn_mesh025.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)





