"""
  Prepare GDAS atm surface temperature for 
  creating predictor fields for ML ithkn emulators
  - Calculate daily mean SAT (T2m)
  - subset North and South polar regions
  - Interpolate onto mesh025 grid

  For this script to work, 6-hr GDAS atm output should be downloaded 

  For deriving several daily means, with fetching GDAS from HPSS, prior to interpolation:
  use sbatch derive_dailyT2m_Ndays.sh --sdate 20250701 --edate 20250731 --dt 3
  This will save daily felds at 3 day intervals for 2025/07

  ML models developed on PPAN

  GFSv17 status with HPSS / WCOSS directories:
  https://docs.google.com/spreadsheets/d/1N3isKTVmE4ITdiULDLP5lK1RoZOzkNlHFN-NFHwrH6o/edit?gid=492588212#gid=492588212

  Fetch GDAS fields from HPSS:
  atm GDAS: scripts/DATA_HPSS/get_GDASatm.sh
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray as xr
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
import mod_sis2_relax as msisrlx
import mod_regmom as mrmom 

init_date = 20250701
init_hr = 0

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", 
                    choices=['north','south','global'],
                    type=str, required=True)
parser.add_argument("--rdate",  
             help="GDAS date(s) YYYYMMDD to calc daily average and interpolate onto mesh025", 
             type=int,
             nargs="+", 
             required=True)
parser.add_argument("--hr", help="use GDAS hour(s) for daily avrg., list, default=[0,6,12,18]",
                    type=int,
                    nargs="+")
parser.add_argument("--tmpf",  
                    help=f"1: Save, start from last processed daily field, default=1", 
                    choices=[0,1],
                    default=1,
                    type=int)
parser.add_argument("--fsave", help=f"Save final dataset with all days as netcdf, default=1", 
                    choices=[0,1], 
                    default=1,
                    type=int)

args     = parser.parse_args()
regn     = args.regn if args.regn else None
RDATES   = args.rdate
HRS      = args.hr if args.hr is not None else [0,6,12,18]
save_nc  = args.fsave == 1
save_tmp = args.tmpf == 1

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

# Get MOM6 grid:
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


def read_gdas(dflgdas, k2c=True, flip_north=True):
  assert os.path.isfile(dflgdas), f"Not found: {dflgdas}"
  with xr.open_dataset(dflgdas) as ds:
    A2d = ds["tmp2m"].isel(time=0).values

  if k2c:
    A2d = A2d - 273.15   # K --> Celsius

  if flip_north:
    # Flip array to have N. at the top:
    A2d = np.flipud(A2d)

  return A2d


def get_gdas_coord(dflgdas, flip_north=True):
  assert os.path.isfile(dflgdas), f"Not found: {dflgdas}"
  with xr.open_dataset(dflgdas) as ds:
    LON = ds["lon"].values
    LAT = ds["lat"].values

  if flip_north:
    # Flip array to have N. at the top:
    LON = np.flipud(LON)
    LAT = np.flipud(LAT)

  return LON, LAT


# Get gmapi 4 GDAS grid points for interpolation
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump = os.path.join(pthdata,'gmapi_mapping2mesh025')
fgmapi = f'GDASatm_reggrid_to_mesh025_gmapi_1440x1080_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)
print(f'Loading gmapi --> {dfgmapi}')

with xr.open_dataset(dfgmapi) as dgmapi:
  IMOM = dgmapi['mom_indx'].data
  JMOM = dgmapi['mom_jndx'].data
  INDX = dgmapi['gmapi_i'].data
  JNDX = dgmapi['gmapi_j'].data

def T2m_daily_mean(pthgdas, HRS):
  """
    Derive daily mean T2m from GDAS
    using HRS hours
  """
  T2m_mean = None
  for ik, hrz in enumerate(HRS):
    flname = f"gdas.t{hrz:02d}z.sfc.f000.nc"
    print(f"Reading {flname}")
    dflgdas = os.path.join(pthgdas, flname)
    T2m = read_gdas(dflgdas)

    if T2m_mean is None:
      T2m_mean = T2m.copy()
    else:
      T2m_mean += T2m

  T2m_mean /= len(HRS)

  return T2m_mean

def write_nc(dfliceout, time_dnmb, A3d):
  yr1, mm1, dd1 = mtime.datevec(time_dnmb[0])[:3]
  nrecs, jdim, idim = A3d.shape
  assert nrecs == len(time_dnmb), f"Mismatch in time axis of A3d ({nrecs}) vs time_dnmb = {len(time_dnmb)}"

  # Days wrt to the 1st day of the month:
  time_days = time_dnmb - mtime.datenum([yr1, mm1, 1])
  darr_t2m = xr.DataArray(A3d, dims=("time","jdim","idim"),\
                     coords={"time": time_days,\
                             "jdim": np.arange(jdm),\
                             "idim": np.arange(idm)})

  dset = xr.Dataset({"temp_2m": darr_t2m})
  dset['temp_2m'].attrs['long_name'] = 'Air temperature C at 2m'
  dset["time"].attrs = {
       "long_name": f"days since {yr1}/{mm1:02d}/{dd1:02d}"
  } 

  # Add global attributes:
  dset.attrs['title']       = 'GDAS air temperature 2m, daily mean, interpolated onto mesh025 grid'
  dset.attrs['source']      = 'calc6hr_interp_GDAS_T2m_daily_mesh025.py'
  dset.attrs['region']      = regn
                    
  print(f'Dumping interpolated daily T2m --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')
  

# Save temporary fields 
if save_tmp:
  tmp_dir = os.path.join(pthdata, f'GDAS_T2m/{regn}','tmp')
  os.makedirs(tmp_dir, exist_ok=True)


A3d = []
time_dnmb = []
for iday, gdas_date in enumerate(RDATES):
  dnmb0 = mtime.rdate2datenum(gdas_date)
  YR, MM, DD = mtime.datevec(dnmb0)[:3]
  if save_tmp:
    tmp_file = os.path.join(tmp_dir, f'GDAS_T2m_mesh025_{YR}{MM:02d}{DD:02d}_{regn}.npy')

  # GDAS data
  pthgdas = f"/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/GFSv17/gdas.{gdas_date}"

  if iday == 0:
    hrz = HRS[0]
    flname = f"gdas.t{hrz:02d}z.sfc.f000.nc"
    print(f"Reading {flname}")
    dflgdas = os.path.join(pthgdas, flname)
    LON, LAT = get_gdas_coord(dflgdas)

  # Load saved tmp fields if exist:
  if save_tmp:
    if os.path.isfile(tmp_file):
      print(f"Loading tmp file: {tmp_file}")
      Tint = np.load(tmp_file)
      A3d.append(Tint)
      time_dnmb.append(dnmb0)
      continue

  # Derive daily mean:
  T2m_day = T2m_daily_mean(pthgdas, HRS)

  # Interpolate to mesh025
  Tint = msisrlx.interp2Dfld(T2m_day, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)

  # Fill the N.Pole hole:
  Tint = mrmom.fill_npole(Tint, hlon, hlat, HH, Rpole=1.)
  Tint = np.where(HH>=0, np.nan, Tint)  # no interpolation has been done over land
  A3d.append(Tint)
  time_dnmb.append(dnmb0)

  if save_tmp:
    print(f"Saving tmp file --> {tmp_file}")
    np.save(tmp_file, Tint)

# Save final:
A3d = np.asarray(A3d)
if iday > 0:
  if not save_nc:
    print(f'Final netcdf is not saved, save_nc={save_nc}')
  else:
    if len(RDATES) == 1: 
      dstmp = f"{RDATES[0]}"
    elif len(RDATES) > 1:
      dstmp = f"{RDATES[0]}_{RDATES[-1]}"

    dnmb0 = mtime.rdate2datenum(RDATES[0])
    YR, MM, DD = mtime.datevec(dnmb0)[:3]

    fliceout = f'GDAS_T2m_interp_mesh025_{jdm}x{idm}_{dstmp}_{regn}.nc'
    pthgdas = os.path.join(pthdata, f'GDAS_T2m/{regn}/{YR}')
    os.makedirs(pthgdas, exist_ok=True)
    dfliceout = os.path.join(pthgdas, fliceout)

    time_dnmb = np.asarray(time_dnmb)

    write_nc(dfliceout, time_dnmb, A3d)


f_chck = True
if f_chck:
  day_chck = 1
  ichck = day_chck - 1 

  AP = A3d[ichck,:,:].squeeze() 

  clrmp = mclrmps.colormap_temperature_coldwarm()
  rmin = -10.
  rmax = 20.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  #AP[HH >= 0] = np.nan   # land
  #AP[np.isnan(AP) & (HH < 0)] = -1.  # ocean

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.1, 0.15, 0.8, 0.8])
      
  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
    parallels = np.arange(-80,-10,10.)
    meridians = np.arange(-360,359.,45.)
  elif regn == 'north':
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
    parallels = np.arange(50, 90, 5)
    meridians = np.arange(-360, 359., 45.)

  xh, yh = m(hlon, hlat)

  m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
  m.drawcoastlines()

  img = m.pcolormesh(xh,yh, AP, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"Air T2m, GDAS on mesh025, {YR}/{MM:02d}/{DD:02d}")

  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='both')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = 'interp_GDAS_T2m_daily_mesh025.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)


