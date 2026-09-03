"""
  Derive dynamic predictor: the number of freeze degree days
  Following Zubov's relation: h2 + 50h = 8 IFDD, 
  IFDD = sum of (Tfrz - Tair), when Tair < Tfrz

  Check whether sqrt(IFDD) is used as a predictor
  Use GDAS T2m daily fields saved every N days
  Interpolated IFDD onto mesh025

  Need mapping indices gmapi to map ERA5 --> GLORYS grid
  Save as numpy binary npy / npz file or netcdf

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
from scipy.interpolate import interp1d

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
parser.add_argument("--dxy", help=f"Min dist (km) between data points (~corr.scale), to skip close i,j points",
                    type=int, required=True)
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--rdate", help="Date of ithkn prediction YYYMMDD", type=int, required=True)
parser.add_argument("--savenc", 
           help="=1: save netcdf full grid, =0: save npy/npz only grid points within LMsk (default)",
           choices=[0,1],
           default=0,
           type=int)
parser.add_argument("--flipN", 
           help="=1: Use flipped South-North GDAS grid to have North at the top, =0 - keep original GDAS grid",
           default=1,
           choices=[0,1],
           type=int)
parser.add_argument("--nintgr", help="Days of IFDD accumulation, default=90",
           default=90,
           type=int)
args  = parser.parse_args()

dxy        = args.dxy
regn       = args.regn
rdate      = args.rdate
intgr_time = args.nintgr
save_nc    = args.savenc == 1
flip_north = args.flipN == 1

Tfrz = -1.85    # ocea freezing T

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

dnmbR = mtime.rdate2datenum(rdate)
YRr, MMr, DDr = mtime.datevec(dnmbR)[:3]

pthdata = pths_ml["GDAS"]["pthdata"]  # root dir for processed data
ptht2m = pths_ml["GDAS"]["pthit2m"].format(regn=regn, YR=YRr)
#regn_name = pths_ml[regn]["name"]
regn_lat0 = pths_ml[regn]["lat0"]  # bounding lat 


def gdas_t2m_file(pths_ml, flip_north, dnmb):
  YR, MM, DD = mtime.datevec(dnmb)[:3]

  # Read GDAS daily mean on GDAS grid:
  pthdaily = pths_ml["GDAS"]["pthdaily"].format(YR=YR)  # daily T2m, GDAS grid
  if flip_north:
    dflt2m = os.path.join(pthdaily, f"GDAS_flipN_T2m_{YR}{MM:02d}{DD:02d}.npy")
  else:
    dflt2m = os.path.join(pthdaily, f"GDAS_notflipN_T2m_{YR}{MM:02d}{DD:02d}.npy")

  return dflt2m

def find_dates(dnmbS, dnmbE, pths_ml, flip_north):
  """
    Find existing T2m daily fields interpoalted onto mesh025
    for the given time range

    For time interpoaltion:
      1st existing time record has to be before or = dnmbS 
      Last record = dnmbE or after that
  """
  DNMB = []

  for dnmb in range(dnmbS, dnmbE+1):
    YR, MM, DD = mtime.datevec(dnmb)[:3]
    dflnm = gdas_t2m_file(pths_ml, flip_north, dnmb)

    found = False
    if os.path.isfile(dflnm):
      DNMB.append(dnmb)
      found = True

    elif dnmb == dnmbS and not found:
      # Find previous closest file for 1st date:
      for ddp in range(dnmb-1, dnmb-10, -1):
        dflnm = gdas_t2m_file(pths_ml, flip_north, ddp)
        if os.path.isfile(dflnm):
          DNMB.append(ddp)
          found = True
          break

      if not found:
        print(f"Last Search: {dflnm}")
        raise RuntimeError(f"Could not find starting record for {mtime.datestr(dnmbS)}")

    elif dnmb == dnmbE and not found:
      for ddp in range(dnmb+1, dnmb+10):
        dflnm = gdas_t2m_file(pths_ml, flip_north, ddp)
        if os.path.isfile(dflnm):
          DNMB.append(ddp)
          found = True
          break

      if not found:
        print(f"Last Search: {dflnm}")
        raise RuntimeError(f"Could not find starting record for {mtime.datestr(dnmbE)}")
         
  return np.asarray(DNMB)

# Dates to process:
dnmbS = dnmbR - intgr_time + 1
dnmbE = dnmbR
DNMB = np.arange(dnmbS, dnmbE+1)

DNMB_data = find_dates(dnmbS, dnmbE, pths_ml, flip_north)

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
  LMsk = np.where(hlat > regn_lat0, 0, LMsk)
elif regn == 'north':
  LMsk = np.where(hlat < regn_lat0, 0, LMsk)

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


def write_nc(dfliceout, dnmbR, A2d, Tfrz, intgr_time):
  yr1, mm1, dd1 = mtime.datevec(dnmbR)[:3]
  time_dnmb = np.asarray([dnmbR])
  A3d = np.expand_dims(A2d, axis=0)

  nrecs, jdim, idim = A3d.shape
  assert nrecs == len(time_dnmb), f"Mismatch in time axis of A3d ({nrecs}) vs time_dnmb = {len(time_dnmb)}"

  # Days wrt to the 1st day of the month:
  time_days = time_dnmb - dnmbR

  darr_t2m = xr.DataArray(A3d, dims=("time","jdim","idim"),
                     coords={"time": time_days,
                             "jdim": np.arange(jdim),
                             "idim": np.arange(idim)})

  dset = xr.Dataset({"IFDD": darr_t2m})
  dset['IFDD'].attrs['long_name'] = 'Integrated Freeze Degree Days, Zubov'
  dset["time"].attrs = {
       "long_name": f"days since {yr1}/{mm1:02d}/{dd1:02d}"
  } 

  # Add global attributes:
  dset.attrs['title'] = f'Intgerated Freeze Degree Days wrt Tfrz={Tfrz:.3f} Ndays={intgr_time} on mesh025'
  dset.attrs['source'] = 'derive_intgrFDD_mesh025.py'
                    
  print(f'Dumping interpolated IFDD --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')
 

def read_gdas_t2m(pths_ml, flip_north, dnmb0):
  # Read GDAS daily T2m field for day dnmb0
  dflt2m = gdas_t2m_file(pths_ml, flip_north, dnmb0)
  assert os.path.isfile(dflt2m), f"File missing: {dflt2m}"
  print(f"Loading T2m: {dflt2m}")
  T2m_day = np.load(dflt2m)

  return T2m_day 

IFDD = None

# Use Zubov's IFDD
for irec0, dnmb0 in enumerate(DNMB):
  YR, MM, DD = mtime.datevec(dnmb0)[:3]

  print(f"{irec0}: Processing GDAS T2m {YR}/{MM:02d}/{DD:02d}")

  #Find bracket days for time interpolation to get T2m for dnmb0
  ii0 = np.searchsorted(DNMB_data, dnmb0, side='right') - 1

  assert 0 <= ii0 < len(DNMB_data), f"{dnmb0} is outside GDAS data range"
  
  if ii0 == len(DNMB_data)-1:
    assert dnmb0 == DNMB_data[ii0], f"Expected last dnmb {dnmb0} matches DNMB_data"

  if dnmb0 == DNMB_data[ii0]:
    # Date matches DNMB_data, no interpolation
    T2m0 = read_gdas_t2m(pths_ml, flip_north, dnmb0)

  else:
    # Interpolate in time
    dnmb1 = DNMB_data[ii0]
    dnmb2 = DNMB_data[ii0+1]
    T2m1  = read_gdas_t2m(pths_ml, flip_north, dnmb1)
    T2m2  = read_gdas_t2m(pths_ml, flip_north, dnmb2)

    # Interpolate lineraly in time
    wght  = (dnmb0 - dnmb1) / (dnmb2 - dnmb1)
    T2m0 = T2m1 + wght * (T2m2 - T2m1)

  #  Use Zubov definition of accumulated freeze degree days
  #  sum(Tfrz-T), when SAT < Tfrz
  # = 0 otherwise
  freeze_degree = np.maximum(Tfrz - T2m0, 0.0)

  if IFDD is None:
    IFDD = freeze_degree.copy()
  else:
    IFDD += freeze_degree


# Interpolate to mesh025
IFDDi = msisrlx.interp2Dfld(IFDD, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)

# Fill the N.Pole hole:
#IFDDi = mrmom.fill_npole(Tint, hlon, hlat, HH, Rpole=1.)
IFDDi = np.where(HH>=0, np.nan, IFDDi)  # no interpolation has been done over land

# Save only active grid points:
JG, IG = np.where(LMsk == 1)

assert np.all(~np.isnan(IFDDi[JG,IG])), "Unexpected nan values in IFDDi at active grid points, check LMsk"

pthprd = pths_ml["PRED"]["pthprd"]
flfrz0 = pths_ml["PRED"]["flfrzdgr"].format(Ndays=intgr_time, dxy=dxy, rdate=rdate, regn=regn)
if save_nc:
  flfrz = f"{flfrz0}.nc"
  dfliceout = os.path.join(pthprd, flfrz)
  write_nc(dfliceout, dnmbR, IFDDi, Tfrz, intgr_time)
else:
  # Save numpy binary
  flfrz = f"{flfrz0}.npz"
  dfliceout = os.path.join(pthprd, flfrz)
  print(f"Saving interpoalted IFDD[JG,IG] file --> {dfliceout}")
  np.savez(dfliceout, IFDD=IFDDi[JG,IG], JG=JG, IG=IG, jdim=jdm, idim=idm)
      
f_chck = True
if f_chck:
  day_chck = 1
  ichck = day_chck - 1 

  clrmp = mclrmps.colormap_temperature_coldwarm()
  rmin = 0.
  rmax = 1000.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  AP = IFDDi.copy()
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
  ax1.set_title(f"IFDD from T2m GDAS on mesh025, {YR}/{MM:02d}/{DD:02d}")

  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='both')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = 'derive_intgrFDD_mesh025.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)


