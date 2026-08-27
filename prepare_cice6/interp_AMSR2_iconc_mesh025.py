"""
  Interpolate AMSR2 sea ice concentration 
  to MOM6/CICE6 mesh025 grid

  get_gmapi_AMSR2_osisaf_to_mesh025

  AMSR2 L4 OSI SAF EUMETSAT sea ice concentration 
  on native grid
  daily fields

  AMSR2 fields from 
  Downloaded from:
  https://data.marine.copernicus.eu/product/SEAICE_ARC_PHY_AUTO_L4_MYNRT_011_024/services
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib  
import xarray as xr
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
import mod_colormaps as mclrmps
import mod_regmom as mrmom 
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--regn",
      help="hemisphere: north or south",
      choices=['north', 'south'],
      required=True)
parser.add_argument("--sdate", help="Start date of AMSR2 interpolation, YYYYMMDD", type=int, required=True)
parser.add_argument("--edate", help="Etart date of AMSR2 interpolation, YYYYMMDD", type=int, required=True)
parser.add_argument("--tmpf",  
                    help=f"1: Save, start from last processed field, default=1", 
                    choices=[0,1],
                    default=1,
                    type=int)
parser.add_argument("--fsave", help=f"Save final dataset with all days as netcdf, default=1", 
                    choices=[0,1], 
                    default=1,
                    type=int)
args = parser.parse_args()
  
regn = args.regn if args.regn else None
sdate = args.sdate
edate = args.edate
tmpf = args.tmpf if args.tmpf is not None else tmpf
fsave = args.fsave if args.fsave is not None else fsave

save_nc = fsave == 1
save_tmp = tmpf == 1

# Dates:
dnmbS = int(mtime.rdate2datenum(sdate))
YRs, MMs, DDs = mtime.datevec(dnmbS)[:3]
dnmbE = int(mtime.rdate2datenum(edate))
YRe, MMe, DDe = mtime.datevec(dnmbE)[:3]

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

with xr.open_dataset(dftopo_mom) as dstopo:
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

# Get gmapi 4 AMSR2 grid points for interpolation
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump = os.path.join(pthdata,'gmapi_mapping2mesh025')
fgmapi = f'AMSR2_to_mesh025_gmapi_{idm}x{jdm}_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)
print(f'Loading gmapi --> {dfgmapi}')

with xr.open_dataset(dfgmapi) as dgmapi:
  IMOM = dgmapi['mom_indx'].data
  JMOM = dgmapi['mom_jndx'].data
  INDX = dgmapi['gmapi_i'].data
  JNDX = dgmapi['gmapi_j'].data

def read_AMSR2(pthdata, regn, dnmb, varnm):
  YR, MM, DD = mtime.datevec(dnmb)[:3]
  pthamsr = os.path.join(pthdata,f"AMSR2_OSISAF_L4/{regn}/{YR}")
  flamsr = f"osisaf_{regn}_nrt_amsr2_L4_{MM:02d}{YR}.nc"
  dfl = os.path.join(pthamsr,flamsr)
  assert os.path.isfile(dfl), f"Missing data: {dfl}"

  if varnm == 'ice_conc':
    with xr.open_dataset(dfl) as ds_amsr:
      AA = ds_amsr[varnm].isel(time=DD-1).values.squeeze() * 0.01  # % to fractions
  elif varnm == 'longitude' or varnm == 'latitude':
    with xr.open_dataset(dfl) as ds_amsr:
      AA = ds_amsr[varnm].values
  else:
    raise ValueError(f"Unknown AMSR2 variable: {varnm}")

  return AA

def write_nc(dfliceout, time_dnmb, A3d):
  yr1, mm1, dd1 = mtime.datevec(time_dnmb[0])[:3]
  nrecs, jdim, idim = A3d.shape
  assert nrecs == len(time_dnmb), f"Mismatch in time axis of A3d ({nrecs}) vs time_dnmb = {len(time_dnmb)}"

  # Days wrt to the 1st day of the month:
  time_days = time_dnmb - mtime.datenum([yr1, mm1, 1])
  darr_cice = xr.DataArray(A3d, dims=("time","jdim","idim"),\
                     coords={"time": time_days,\
                             "jdim": np.arange(jdm),\
                             "idim": np.arange(idm)})

  dset = xr.Dataset({"ice_conc": darr_cice})
  dset['ice_conc'].attrs['long_name'] = 'ice partial area'
  dset["time"].attrs = {
       "long_name": f"days since {yr1}/{mm1:02d}/{dd1:02d}"
  } 

  # Add global attributes:
  dset.attrs['title']       = 'AMSR2 L4 OSI SAF EUMETSAT sea ice concentration daily interpolated onto mesh025 grid'
  dset.attrs['institution'] = 'NOAA NWS OMD'
  dset.attrs['source']      = 'interp_AMSR2_iconc_mesh025.py'
  dset.attrs['region']      = regn
  dset.attrs['Grid_idm_jdm'] = f'{idm}x{jdm}'

  print(f'Dumping interpolated ice conc --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')


# Get lon/lat for AMSR2 data
lon1d = read_AMSR2(pthdata, regn, dnmbS, 'longitude')
lat1d = read_AMSR2(pthdata, regn, dnmbS, 'latitude')

LON, LAT = np.meshgrid(lon1d, lat1d)
JDIM, IDIM = LON.shape

# Save temporary fields 
if save_tmp:
  tmp_dir = os.path.join(pthdata, f'AMSR2_OSISAF_L4/{regn}','tmp')
  os.makedirs(tmp_dir, exist_ok=True)


# Save by months
MMwrk = 0
YRwrk = 0
irec = 0 
A3d  = []
time_dnmb = []
for dnmb0 in range(dnmbS, dnmbE+1):
  YR, MM, DD = mtime.datevec(dnmb0)[:3]

  if not MM == MMwrk or not YR == YRwrk:
    # Save previous month fields:
    if irec > 0:
      # Convert lists --> np array AAin etc
      # Save processes fields    

      if not save_nc:
        print(f"Final netcdf is not saved, save_nc={save_nc}")
       
      else:
        fliceout = f'AMSR2_iconc_interp_mesh025_{jdm}x{idm}_{YRwrk}{MMwrk:02d}_{regn}.nc'
        pthamsr = os.path.join(pthdata, f'AMSR2_OSISAF_L4/{regn}/{YRwrk}')
        os.makedirs(pthamsr, exist_ok=True)
        dfliceout = os.path.join(pthamsr, fliceout)

        A3d = np.asarray(A3d)
        time_dnmb = np.asarray(time_dnmb)

        write_nc(dfliceout, time_dnmb, A3d)

    # Prepare fields for processing new month: 
    MMwrk = MM
    YRwrk = YR
    irec = 0
    A3d = []
    time_dnmb = []

  # Load saved tmp fields if exist:
  
  if save_tmp:
    tmp_file = os.path.join(tmp_dir, f'AMSR2_iconc_mesh025_{YR}{MM:02d}{DD:02d}_{regn}.npy')
    if os.path.isfile(tmp_file):
      print(f"Loading tmp file: {tmp_file}")
      CIint = np.load(tmp_file)
      A3d.append(CIint)
      time_dnmb.append(dnmb0)
      irec += 1
      continue
    
  AA = read_AMSR2(pthdata, regn, dnmb0, 'ice_conc') 

  CIint = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
  # Fill the N.Pole hole:
  CIint = mrmom.fill_npole(CIint, hlon, hlat, HH, Rpole=1.)
  CIint = np.where(HH>=0, np.nan, CIint)
  A3d.append(CIint)
  time_dnmb.append(dnmb0)
  irec += 1

  if save_tmp:
    tmp_file = os.path.join(tmp_dir, f'AMSR2_iconc_mesh025_{YR}{MM:02d}{DD:02d}_{regn}.npy')
    print(f"Saving tmp file --> {tmp_file}")
    np.save(tmp_file, CIint)

# Save final month:
if irec > 0:
	if not save_nc:
		print(f'Final netcdf is not saved, save_nc={save_nc}')
	else:
		fliceout = f'AMSR2_iconc_interp_mesh025_{jdm}x{idm}_{YRwrk}{MMwrk:02d}_{regn}.nc'
		pthamsr = os.path.join(pthdata, f'AMSR2_OSISAF_L4/{regn}/{YRwrk}')
		os.makedirs(pthamsr, exist_ok=True)
		dfliceout = os.path.join(pthamsr, fliceout)

		A3d = np.asarray(A3d)
		time_dnmb = np.asarray(time_dnmb)

		write_nc(dfliceout, time_dnmb, A3d)

f_chck = True
if f_chck:
  clrmp = mclrmps.colormap_conc()
  rmin = 0.
  rmax = 1.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  clrmp.set_under(color=[1,1,1]) 

  AP = CIint.copy()
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
  ax1.set_title(f"iconc, AMSR2, {YR}/{MM:02d}/{DD:02d}")

  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = 'interp_AMSR2_iconc_mesh025'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)



