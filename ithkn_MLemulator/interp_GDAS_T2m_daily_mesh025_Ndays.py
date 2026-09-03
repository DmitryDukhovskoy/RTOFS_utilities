"""
  Interpolate GDAS atm surface temperature (T2m) onto mesh025
  by days, save individual days
  Save as numpy binary npy / npz file

  For this script to work, 6-hr GDAS atm output should be downloaded 
  and averaged by days, saved as numpy binary

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

init_hr = 0

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", 
                    choices=['north','south'],
                    type=str, required=True)
parser.add_argument("--sdate",  
             help="GDAS start date YYYYMMDD to interpolate onto mesh025", 
             type=int,
             required=True)
parser.add_argument("--edate",  
             help="GDAS end date  YYYYMMDD to interpolate onto mesh025, default=sdate (1day)", 
             type=int,
             )
parser.add_argument("--savenc", help="=1: save netcdf daily file, =0: save npy/npz (default)",
           choices=[0,1],
           default=0,
           type=int)
parser.add_argument("--flipN", 
           help="=1: Use flipped South-North GDAS grid to have North at the top, =0 - keep original GDAS grid",
           default=1,
           choices=[0,1],
           type=int)

args     = parser.parse_args()
regn     = args.regn if args.regn else None
sdate    = args.sdate
edate    = args.edate if args.edate is not None else args.sdate
save_nc  = args.savenc == 1
flip_north = args.flipN == 1

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

fyaml = 'paths_ML.yaml'
with open(fyaml) as ff:
  pths_ml = safe_load(ff)

pthdata = pths_ml["GDAS"]["pthdata"]  # root dir for processed data
pthdaily = pths_ml["GDAS"]["pthdaily"]  # daily GDAS fields

def derive_time(dnmbS, dnmbE, ptht2m):
  DNMB = []
  time_stmp = []
  for dnmb in range(dnmbS, dnmbE+1):
    YR, MM, DD = mtime.datevec(dnmb)[:3]
    flnm = f"GDAS_flipN_T2m_{YR}{MM:02d}{DD:02d}.npy"
    dflnm = os.path.join(ptht2m, flnm)
    if os.path.isfile(dflnm):
      DNMB.append(dnmb)

  return np.asarray(DNMB)

# Dates to process:
dnmbS = mtime.rdate2datenum(sdate)
dnmbE = mtime.rdate2datenum(edate)
DNMB = derive_time(dnmbS, dnmbE, pthdaily)


# Get MOM6 grid:
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
     
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
  
with xr.open_dataset(dftopo_mom) as dstopo:
  depth = dstopo['depth'].values.squeeze()
  
# Convert all positive values -> land (100) and ocean (<0):
HH = np.where(depth < 1.e-20, 100., -depth)
 
jdm, idm = HH.shape
LMsk = np.where(HH<0, 1, 0)
    
# Mask out not needed latitudes:
if regn == 'south':
  LMsk = np.where(hlat > -55, 0, LMsk)
elif regn == 'north':
  LMsk = np.where(hlat < 55, 0, LMsk)


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


# Get gmapi 4 GDAS grid points for interpolation
pthdata = pths_ml["GDAS"]["pthdata"]
pthdump = os.path.join(pthdata,'gmapi_mapping2mesh025')
fgmapi = f'GDASatm_reggrid_to_mesh025_gmapi_1440x1080_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)
print(f'Loading gmapi --> {dfgmapi}')

with xr.open_dataset(dfgmapi) as dgmapi:
  LON = dgmapi['longit'].values      # GDAS reg. grid
  LAT = dgmapi['latit'].values       # GDAS reg. grid
  IMOM = dgmapi['mom_indx'].data
  JMOM = dgmapi['mom_jndx'].data
  INDX = dgmapi['gmapi_i'].data
  JNDX = dgmapi['gmapi_j'].data


def write_nc(dfliceout, time_dnmb, A3d):
  yr1, mm1, dd1 = mtime.datevec(time_dnmb[0])[:3]

  nrecs, jdim, idim = A3d.shape
  assert nrecs == len(time_dnmb), f"Mismatch in time axis of A3d ({nrecs}) vs time_dnmb = {len(time_dnmb)}"

  # Days wrt to the 1st day of the month:
  time_days = time_dnmb - mtime.datenum([yr1, mm1, 1])

  darr_t2m = xr.DataArray(A3d, dims=("time","jdim","idim"),
                     coords={"time": time_days,
                             "jdim": np.arange(jdim),
                             "idim": np.arange(idim)})

  dset = xr.Dataset({"temp_2m": darr_t2m})
  dset['temp_2m'].attrs['long_name'] = 'Air temperature C at 2m'
  dset["time"].attrs = {
       "long_name": f"days since {yr1}/{mm1:02d}/{dd1:02d}"
  } 

  # Add global attributes:
  dset.attrs['title']       = 'GDAS air temperature 2m, daily mean, interpolated onto mesh025 grid'
  dset.attrs['source']      = 'interp_GDAS_T2m_daily_mesh025_Ndays.py'
  dset.attrs['region']      = regn
                    
  print(f'Dumping interpolated daily T2m --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')
  

for iday, dnmb0 in enumerate(DNMB):
  A3d = []
  time_dnmb = []

  YR, MM, DD = mtime.datevec(dnmb0)[:3]

  # Read GDAS daily mean:
  if flip_north:
    dflt2m = os.path.join(pthdaily, f"GDAS_flipN_T2m_{YR}{MM:02d}{DD:02d}.npy")
  else:
    dflt2m = os.path.join(pthdaily, f"GDAS_notflipN_T2m_{YR}{MM:02d}{DD:02d}.npy")

  assert os.path.isfile(dflt2m), f"File missing: {dflt2m}"

  print(f"Loading T2m: {dflt2m}")
  T2m_day = np.load(dflt2m)

  # Interpolate to mesh025
  Tint = msisrlx.interp2Dfld(T2m_day, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)

  # Fill the N.Pole hole:
  Tint = mrmom.fill_npole(Tint, hlon, hlat, HH, Rpole=1.)
  Tint = np.where(HH>=0, np.nan, Tint)  # no interpolation has been done over land

  A3d.append(Tint)
  time_dnmb.append(dnmb0)

  A3d = np.asarray(A3d)
  time_dnmb = np.asarray(time_dnmb)

  pthintrp = pths_ml["GDAS"]["pthit2m"].format(regn=regn, YR=YR)
  if save_nc:
    flt2m_out = f"GDAS_t2m_interp_mesh025_1080x1440_{YR}{MM:02d}{DD:02d}_{regn}.nc"
    dfliceout = os.path.join(pthintrp, flt2m_out)
    write_nc(dfliceout, time_dnmb, A3d)
  else:
    # Save numpy binary
    flt2m_out = f"GDAS_t2m_interp_mesh025_1080x1440_{YR}{MM:02d}{DD:02d}_{regn}.npz"
    dfliceout = os.path.join(pthintrp, flt2m_out)
    print(f"Saving T2m file --> {dfliceout}")
    np.savez(dfliceout, A3d=A3d, TM=time_dnmb)
        
f_chck = True
if f_chck:
  day_chck = 1
  ichck = day_chck - 1 

  AP = A3d[ichck,:,:].squeeze() 

  clrmp = mclrmps.colormap_temperature_coldwarm()
  clrmp = mclrmps.colormap_difference_negpos()
  rmin = -20.
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

  btx = 'interp_GDAS_T2m_daily_mesh025_Ndays.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)


