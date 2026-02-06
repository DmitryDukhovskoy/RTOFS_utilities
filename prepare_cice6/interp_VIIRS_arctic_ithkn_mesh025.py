"""
  Interpolate VIIRS ice thickness to mesh025

  ice thickness Arctic
  daily fields, 750-m resolution

  1 or several days
  using OpeNDaP or from downloaded files (not recommended - too big)

https://coastwatch.noaa.gov/cwn/products/viirs-sea-ice-concentration-ice-thickness-ice-surface-temperature.html

Use NSIDC Polar Stereogr Proj:
semi_major_axis: 6378137.0
inverse_flattening: 298.257223563

  Due to very large VIIRS grid, recommended to run several serial jobs saving temporary files with
  gmapi indices, using get_gmapi_subVIIRS_arctic_to_mesh025.py  e.g.
  run get_gmapi_subVIIRS_arctic_to_mesh025.py --iS0 0 --iE0 200
  run get_gmapi_subVIIRS_arctic_to_mesh025.py --iS0 201 --iE0 400
  ...
  
  then use this script to combine all saved pieces and creating the final netcdf

  get_gmapi_VIIRS_arctic_to_mesh025.py

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
YR    = 2021
ndays = 1
dskip = 1

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help=f"Year of data, default={YR}", type=int)
parser.add_argument("--mm", help=f"Month of data", type=int)
parser.add_argument("--dd", help=f"Month day of data", type=int)
parser.add_argument("--jday", help=f"Year day, instead of MM/DD", type=int)
parser.add_argument("--ndays", help=f"How many days to process, default=1", type=int)
parser.add_argument("--dskip", help=f"N days to skip when processing, default={dskip}", type=int)
parser.add_argument("--furl", help="Read data from VIIRS website (=1) or use downloaded (=0)",
                    choices=[0,1], required=True, type=int)
args = parser.parse_args()
  
YR    = args.yr    if args.yr   else YR
MM    = args.mm    if args.mm   else None
DD    = args.dd    if args.dd   else None
jday  = args.jday  if args.jday else None
ndays = args.ndays if args.ndays else ndays
dskip = args.dskip if args.dskip else dskip
furl  = args.furl  if args.furl is not None else None

if jday is None and (MM is None or DD is None):
  raise ValueError(f"Either jday or MM and DD have to be specified")

# Set initial dates:
if jday is None:
  jday = int(mtime.date2jday([YR,MM,DD]))
  dnmb0 = mtime.jday2dnmb(YR, jday)
else:
  dnmb0 = mtime.jday2dnmb(YR, jday)

YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]

# List of days to process:
DAYS_LIST = []
for idd in range(0,ndays,dskip):
  dnmb = dnmb0 + idd
  DAYS_LIST.append(dnmb)
  YR, MM, DD = mtime.datevec(dnmb)[:3]
  print(f"Interpolation will be done:     {YR}/{MM}/{DD}")
read_url = bool(furl)

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





