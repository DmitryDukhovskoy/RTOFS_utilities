"""
  Derive dynamic predictor: the number of freeze degree days
  Following Zubov's relation: h2 + 50h = 8 IFDD, 
  IFDD = sum of (Tfrz - Tair), when Tair < Tfrz

  Heat degree days
  in analogy to Zubov's IFDD, but opposite
  which may be important for melt season

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

  ML models developed / trained on PPAN

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


parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--rdate", help="Date of ithkn prediction YYYMMDD", type=int, required=True)
parser.add_argument(
  "--field",
  help="Output field: ifdd (Freeze Degr. Days) or ihdd (Heat Degr. Days), or both",
  choices=['ifdd','ihdd','all'],
  required=True
)
parser.add_argument(
  "--savenc", 
  help="=1: save netcdf full grid, =0: save npy/npz only grid points within LMsk (default)",
  choices=[0,1],
  default=0,
  type=int
)
parser.add_argument(
  "--flipN", 
  help="=1: Use flipped South-North GDAS grid to have North at the top, =0 - keep original GDAS grid",
  default=1,
  choices=[0,1],
  type=int
)
args  = parser.parse_args()

regn       = args.regn
rdate      = args.rdate
fld_name   = args.field.lower()
save_nc    = args.savenc == 1
flip_north = args.flipN == 1

if fld_name == 'all':
  save_ifdd = True
  save_ihdd = True
else:
  save_ifdd = fld_name == 'ifdd'
  save_ihdd = fld_name == 'ihdd'

dnmbR = mtime.rdate2datenum(rdate)
YRr, MMr, DDr = mtime.datevec(dnmbR)[:3]

fyaml = "config_ML.yaml"
with open(fyaml) as ff:
  config_ml = safe_load(ff)

pthdata   = config_ml["GDAS"]["pthdata"]  # root dir for processed data
regn_lat0 = config_ml["regn"][regn]["lat_bnd"]
dxy       = config_ml["params"]["dxy"]             # ice length scale
Tfrz      = config_ml["params"]["Tfrz"]
intgr_frz = config_ml["params"]["intgr_frz"]

def gdas_t2m_file(config_ml, flip_north, dnmb):
  YR, MM, DD = mtime.datevec(dnmb)[:3]

  # Read GDAS daily mean on GDAS grid:
  pthdaily = config_ml["GDAS"]["pthdaily"].format(YR=YR)  # daily T2m, GDAS grid
  if flip_north:
    dflt2m = os.path.join(pthdaily, f"GDAS_flipN_T2m_{YR}{MM:02d}{DD:02d}.npy")
  else:
    dflt2m = os.path.join(pthdaily, f"GDAS_notflipN_T2m_{YR}{MM:02d}{DD:02d}.npy")

  return dflt2m

def find_dates(dnmbS, dnmbE, config_ml, flip_north):
  """
    Find existing T2m daily fields interpoalted onto mesh025
    for the given time range

    For time interpoaltion:
      1st existing time record has to be before or = dnmbS 
      Last record = dnmbE or after that
  """
  DNMB = []

  for dnmb in range(dnmbS, dnmbE+1):
    dflnm = gdas_t2m_file(config_ml, flip_north, dnmb)

    found = False
    if os.path.isfile(dflnm):
      DNMB.append(dnmb)
      found = True

    elif dnmb == dnmbS and not found:
      # Find previous closest file for 1st date:
      for ddp in range(dnmb-1, dnmb-10, -1):
        dflnm = gdas_t2m_file(config_ml, flip_north, ddp)
        if os.path.isfile(dflnm):
          DNMB.append(ddp)
          found = True
          break

      if not found:
        print(f"Last Search: {dflnm}")
        raise RuntimeError(f"Could not find starting record for {mtime.datestr(dnmbS)}")

    elif dnmb == dnmbE and not found:
      for ddp in range(dnmb+1, dnmb+10):
        dflnm = gdas_t2m_file(config_ml, flip_north, ddp)
        if os.path.isfile(dflnm):
          DNMB.append(ddp)
          found = True
          break

      if not found:
        print(f"Last Search: {dflnm}")
        raise RuntimeError(f"Could not find starting record for {mtime.datestr(dnmbE)}")
         
  return np.asarray(DNMB)

# Dates to process:
dnmbS = int(dnmbR - intgr_frz + 1)
dnmbE = int(dnmbR)
DNMB = np.arange(dnmbS, dnmbE+1)

DNMB_data = find_dates(dnmbS, dnmbE, config_ml, flip_north)
assert DNMB_data[0] <= DNMB[0], "Missing GDAS record before integration period"
assert DNMB_data[-1] >= DNMB[-1], "Missing GDAS record after integration period"
assert np.all(np.diff(DNMB_data) > 0), "DNMB_data must be strictly increasing"

# Get MOM6 grid:
pthgrid = config_ml["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
     
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
  
with xr.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

# Only for positive bottom depths:
assert np.min(HH) >= 0, "Bottom depth expected > 0"  
HH = np.where(HH < 1.e-20, np.nan, -HH)
HH = np.where(np.isnan(HH), 1., HH)
  
jdm, idm = HH.shape
LMsk = np.where(HH<0, 1, 0)
    
# Mask out not needed latitudes:
if regn == 'south':
  LMsk = np.where(hlat > regn_lat0, 0, LMsk)
elif regn == 'north':
  LMsk = np.where(hlat < regn_lat0, 0, LMsk)

# Get gmapi 4 GDAS grid points for interpolation
pthdata = config_ml["GDAS"]["pthdata"]
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


def debug_ifdd(ipp, Tsurf, hlat, hlon, SATprv, Tfrz, fld_pnts):
  """
  Debugging Intgr Freeze Dgr Days
  """
  T2m_ij = np.asarray(Tsurf)
  LATI = hlat[JG,IG]
  LONI = hlon[JG,IG]
  xpp = LONI[ipp]
  ypp = LATI[ipp]
  t2m_prv = SATprv[ipp,:]
  ndays_frz = np.count_nonzero(t2m_prv < Tfrz)
  ihdd = np.nansum(Tfrz - t2m_prv[t2m_prv < Tfrz])
  print(f"  Check pnt: x={xpp:.2f}W, y={ypp:.2f}N, N integr. days: {intgr_frz}")
  print(f"  N days T < {Tfrz:.2f}: {ndays_frz}, min/max T: {np.min(t2m_prv):.2f} / {np.max(t2m_prv):.2f}"
        f" IntgrFreeze: {ihdd:.2f}")
  print(f"  ALL: Min / max T2m: {np.min(T2m_ij):.2f} / {np.max(T2m_ij):.2f}")
  print(f"  ALL: min/max IFDD: {np.min(fld_pnts):.1f} / {np.max(fld_pnts):.1f}")


def debug_ihdd(ipp, Tsurf, hlat, hlon, SATprv, Tfrz, fld_pnts):
  """
  Debugging Intgr Heat Dgr Days
  """
  T2m_ij = np.asarray(Tsurf)
  LATI = hlat[JG,IG]
  LONI = hlon[JG,IG]
  xpp = LONI[ipp]
  ypp = LATI[ipp]
  t2m_prv = SATprv[ipp,:]
  ndays_pos = np.count_nonzero(t2m_prv > Tfrz)
  ihdd = np.nansum(t2m_prv[t2m_prv > Tfrz] - Tfrz)
  print(f"  Check pnt: x={xpp:.2f}W, y={ypp:.2f}N, N integr. days: {intgr_frz}")
  print(f"  N days T > {Tfrz:.2f}: {ndays_pos}, min/max T: {np.min(t2m_prv):.2f} / {np.max(t2m_prv):.2f}"
        f" IntgrHeat: {ihdd:.2f}")
  print(f"  ALL: Min / max T2m: {np.min(T2m_ij):.2f} / {np.max(T2m_ij):.2f}")
  print(f"  ALL: min/max IHDD: {np.min(fld_pnts):.1f} / {np.max(fld_pnts):.1f}")


def write_nc(dfliceout, dnmbR, A2d, Tfrz, intgr_frz, fld_save):
  yr1, mm1, dd1 = mtime.datevec(dnmbR)[:3]
  time_dnmb = np.asarray([dnmbR])
  A3d = np.expand_dims(A2d, axis=0)

  nrecs, jdim, idim = A3d.shape
  assert nrecs == len(time_dnmb), f"Mismatch in time axis of A3d ({nrecs}) vs time_dnmb = {len(time_dnmb)}"

  FLD = fld_save.upper()
  if fld_save == 'ifdd':
    fld_long = 'Intgerated Freeze Degree Days'
  elif fld_save == 'ihdd':
    fld_long = 'Intgerated Heat Degree Days'
   
  # Days wrt to the 1st day of the month:
  time_days = time_dnmb - dnmbR

  darr_t2m = xr.DataArray(A3d, dims=("time","jdim","idim"),
                     coords={"time": time_days,
                             "jdim": np.arange(jdim),
                             "idim": np.arange(idim)})

  dset = xr.Dataset({FLD: darr_t2m})
  dset[FLD].attrs['long_name'] = fld_long
  dset["time"].attrs = {
       "long_name": f"days since {yr1}/{mm1:02d}/{dd1:02d}"
  } 

  # Add global attributes:
  dset.attrs['title'] = f'{fld_long} wrt Tfrz={Tfrz:.3f} Ndays={intgr_frz} on mesh025'
  dset.attrs['source'] = 'derive_IHDD_IFDD_mesh025.py'
                    
  print(f'Dumping interpolated {FLD} --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')
 

def read_gdas_t2m(config_ml, flip_north, dnmb0):
  # Read GDAS daily T2m field for day dnmb0
  dflt2m = gdas_t2m_file(config_ml, flip_north, dnmb0)
  assert os.path.isfile(dflt2m), f"File missing: {dflt2m}"
  print(f"Loading T2m: {dflt2m}")
  T2m_day = np.load(dflt2m)

  return T2m_day 

IFDD = None
IHDD = None

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
    T2m0 = read_gdas_t2m(config_ml, flip_north, dnmb0)

  else:
    # Interpolate in time
    dnmb1 = DNMB_data[ii0]
    dnmb2 = DNMB_data[ii0+1]
    T2m1  = read_gdas_t2m(config_ml, flip_north, dnmb1)
    T2m2  = read_gdas_t2m(config_ml, flip_north, dnmb2)

    # Interpolate lineraly in time
    wght  = (dnmb0 - dnmb1) / (dnmb2 - dnmb1)
    T2m0 = T2m1 + wght * (T2m2 - T2m1)

  #  Use Zubov definition of accumulated freeze degree days
  #  sum(Tfrz-T), when SAT < Tfrz
  # = 0 otherwise
  # Similarly, integr. heat degree days
  freeze_degree = np.maximum(Tfrz - T2m0, 0.0)
  heat_degree = np.maximum(T2m0 - Tfrz, 0.0)

  if IFDD is None:
    IFDD = freeze_degree.copy()
  else:
    IFDD += freeze_degree

  if IHDD is None:
    IHDD = heat_degree.copy()
  else:
    IHDD += heat_degree

# Save only active grid points:
JG, IG = np.where(LMsk == 1)

if save_ifdd: 
  # Interpolate to mesh025
  IFDDi = msisrlx.interp2Dfld(IFDD, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)

  # Fill the N.Pole hole:
  #IFDDi = mrmom.fill_npole(Tint, hlon, hlat, HH, Rpole=1.)
  IFDDi = np.where(HH>=0, np.nan, IFDDi)  # no interpolation has been done over land

  #assert np.all(~np.isnan(IFDDi[JG,IG])), "Unexpected nan values in IFDDi at active grid points, check LMsk"
  assert np.all(np.isfinite(IFDDi[JG, IG])), "Non-finite IFDD values at active grid points"

  pthprd = config_ml["PRED"]["pthprd"]
  flfrz0 = config_ml["PRED"]["flifdd"].format(Ndays=intgr_frz, dxy=dxy, rdate=rdate, regn=regn)
  if save_nc:
    flfrz = f"{flfrz0}.nc"
    dfliceout = os.path.join(pthprd, flfrz)
    write_nc(dfliceout, dnmbR, IFDDi, Tfrz, intgr_frz, "ifdd")
  else:
    # Save numpy binary
    flfrz = f"{flfrz0}.npz"
    dfliceout = os.path.join(pthprd, flfrz)
    print(f"Saving interpoalted IFDD[JG,IG] file --> {dfliceout}")
    np.savez(dfliceout, IFDD=IFDDi[JG,IG], JG=JG, IG=IG, jdim=jdm, idim=idm)

if save_ihdd:      
  # Interpolate to mesh025
  IHDDi = msisrlx.interp2Dfld(IHDD, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)

  IHDDi = np.where(HH>=0, np.nan, IHDDi)  # no interpolation has been done over land

  #assert np.all(~np.isnan(IHDDi[JG,IG])), "Unexpected nan values in IHDDi at active grid points, check LMsk"
  assert np.all(np.isfinite(IHDDi[JG, IG])), "Non-finite IHDDi values at active grid points"

  pthprd = config_ml["PRED"]["pthprd"]
  flfrz0 = config_ml["PRED"]["flihdd"].format(Ndays=intgr_frz, dxy=dxy, rdate=rdate, regn=regn)
  if save_nc:
    flfrz = f"{flfrz0}.nc"
    dfliceout = os.path.join(pthprd, flfrz)
    write_nc(dfliceout, dnmbR, IHDDi, Tfrz, intgr_frz, "ihdd")
  else:
    # Save numpy binary
    flfrz = f"{flfrz0}.npz"
    dfliceout = os.path.join(pthprd, flfrz)
    print(f"Saving interpoalted IHDD[JG,IG] file --> {dfliceout}")
    np.savez(dfliceout, IHDD=IHDDi[JG,IG], JG=JG, IG=IG, jdim=jdm, idim=idm)


f_chck = True
if f_chck:
  day_chck = 1
  ichck = day_chck - 1 
  #fld_chck = 'ihdd'
  fld_chck = 'ifdd'

  clrmp = mclrmps.colormap_temperature_coldwarm()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  if fld_chck == 'ihdd':
    AP = IHDDi.copy()
    rmin = 0.
    rmax = 300.
  elif fld_chck == 'ifdd':
    AP = IFDDi.copy()
    rmin = 0.
    rmax = 1000.
  
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
  ax1.set_title(f"{fld_chck} from T2m GDAS on mesh025, {YR}/{MM:02d}/{DD:02d}")

  ax3 = fig1.add_axes([0.1, 0.1, 0.8, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='both')
  ticks = np.linspace(rmin, rmax, 11)
  clb.set_ticks(ticks)
  clb.ax.set_xticklabels([f"{t:.0f}" for t in ticks])
  clb.ax.tick_params(direction='in', length=12)

  btx = 'derive_IHDD_IFDD_mesh025.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)


