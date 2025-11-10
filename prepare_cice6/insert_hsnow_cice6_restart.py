"""
  Modify snow depth on sea ice in the CICE6 restart file
  by direct insertion of snow depth fields

  Assumed that all non-nan non-zero values are inserted 
  to the grid values where aice > 0

  Snow is distributed across the thikn. categories proportional 
  to the aice (ice partial area)

  Here, snow depth climatology (1998-2007) from NASA SSM/I gridded fields
  are used

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
import argparse

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
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

rest_date = 20250103
rest_hr   = 0
hunits    = 'cm'

parser = argparse.ArgumentParser()
parser.add_argument("--rdate", help=f"restart date input file, default={rest_date}", type=int)
parser.add_argument("--rhr", help=f"input file, restart hour = 0, ..., 23, default={rest_hr}", type=int)
parser.add_argument("--rdate_out", help="output file, restart date if different from input", type=int)
parser.add_argument("--rhr_out", help="output file, restart hour if date is different from input", type=int)
args = parser.parse_args()

rest_date     = args.rdate if args.rdate else rest_date
rest_hr       = args.rhr if args.rhr else rest_hr
rest_date_out = args.rdate_out if args.rdate_out else rest_date
rest_hr_out   = args.rhr_out if args.rhr_out else rest_hr

change_rest_time = (rest_date != rest_date_out) or (rest_hr != rest_hr_out)

# Get date numbers:
# Input restart file
dnmbR = mtime.rdate2datenum(rest_date*100+rest_hr)  # restart day nmb
yrR,mmR,ddR,hrR = mtime.datevec(dnmbR, round_hrs=True)[:4]
nsecR = hrR*3600

# Dates of the output fields in the new restart:
dnmbN = mtime.rdate2datenum(rest_date_out*100+rest_hr_out)
yrN, mmN, ddN, hrN = mtime.datevec(dnmbN, round_hrs=True)[:4]
nsecN = hrN*3600
 
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


def extract_suffix(fname):
  parts = fname.split('.')
  # must end with .nc, so suffix is the part before that
  if len(parts) > 2 and parts[-1] == 'nc':
    suffix = parts[-2]
    if not suffix.isdigit():
      return suffix
  return None

pthrest = os.path.join(pths_ufs[node_nm]["MOM6"]["pthrest"],'new')
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]

# CICE parameters:
puny      = 1.e-11
c0        = 0.0
c1        = 1.0
c2        = 2.0
p5        = 0.5
Lsub      = 2.835e6    # latent heat sublimation fw (J/kg)
Lvap      = 2.501e6    # latent heat vaporization fw (J/kg)
Lfresh    = Lsub - Lvap # latent heat of melting of fresh ice (J/kg)
cp_ice    = 2106.       # specific heat of fresh ice (J/ kg/K)
rhos      = 330.        # density of snow (kg/m3)
hs_min    = 1.e-4       # min snow thickness for computing Tsno (m)
nsal      = 0.407
msal      = 0.573
min_salin = 0.1      # threshold for brine pocket treatment
saltmax   = 3.2        # max S at ice base
hg        = 1.e20    # bad values, land mask, etc.
nslyr     = 1   # snow layers
Tmin      = -100.   # minimum snow T

# Snow depth climatology, Interpolated fields mesh025:
pthsnow = os.path.join(pthdata,'snow_nasa','monthly_clim')
flhsn = 'SSMI_hsnow_mnthclim_1998_2007_mesh025_1440x1080_south.nc'
dflhsn = os.path.join(pthsnow,flhsn)
print(f"Reading interpolated hsnow {dflhsn}")
with xarray.open_dataset(dflhsn) as ds_snow:
  HSi = ds_snow['snow_depth'].isel(time=mmN-1).data.squeeze()
  LON = ds_snow['lon'].data
  LAT = ds_snow['lat'].data
  units = ds_snow['snow_depth'].attrs.get('units', None)
  if units is not None:
    print(f"'snow_depth' units: {units}")
    hunits = units
  else:
    print("No 'units' attribute found for 'snow_depth', use default: {hunits}")

units_m = hunits == 'm'

# Restart from a GFS17 rt13  run:
#flrst_in = f"cice_model.res.{yrR}{mmR:02d}{ddR:02d}.{nsecR:06d}.nc"
# Restart with inserted iconc from NSIDC NRT:
flrst_in = f"cice_model.res.{yrR}{mmR:02d}{ddR:02d}.{hrR:02d}.iconc.nc"
dflrst_in = os.path.join(pthrest, flrst_in)
print(f"Reading restart: {dflrst_in}")
ds_in = xarray.open_dataset(dflrst_in)
ds_out = ds_in.copy(deep=True)
ds_in.close()


# Insert snow:
assert nslyr == 1, f"Code needs to be modified for nslyr>1, nslyr={nslyr}"
aicen = ds_out['aicen'].data  # partial area by cats
vsnon = ds_out['vsnon'].data  # snow vol per m2 of ice area
qsnon = ds_out['qsno001'].data  # snow enthalpy by cats for 1 snow layer
vicen = ds_out['vicen'].data   # ice vol per unit area of grid cell m3/m2
ncat, jdim, idim = vsnon.shape

# Aggregated ice partial area:
aice = np.sum(aicen, axis=0).squeeze()

# Select points to insert:
Jins, Iins = np.where((HSi > puny) & (~np.isnan(HSi)) & (aice > puny))
Xins = LON[Jins,Iins]
Yins = LAT[Jins,Iins]
npnts = len(Jins)

print(f"Found {npnts} points for insertion, min/max lat={np.min(Yins):.1f}/{np.max(Yins):.1f}"
       f" lon={np.min(Xins):.1f}/{np.max(Xins):.1f}")

# Note qsnon, qice < 0 !
vsnon_new = vsnon.astype(ds_out['vsnon'].dtype).copy()
qsnon_new = qsnon.astype(ds_out['qsno001'].dtype).copy()

dvol_sum = 0.
print("Snow depth insertion ...")
for ipp in range(npnts):
  if ipp%10000 == 0:
    print(f"   {ipp/npnts*100.:.2f}% done ...")
  j0 = Jins[ipp]
  i0 = Iins[ipp]

  # Note hsnow = vsn / aice for aice > 0
  # for cat n: vsn(n) = hsnow(n) * aice(n) 
  ai  = aice[j0,i0]           # aggreageted ice partial area 
  ain = aicen[:,j0,i0]        # partial areas by cats
  vin = vicen[:,j0,i0]
  vsn = vsnon[:,j0,i0]        # snow volume per unit grid-cell area m2
  if units_m:
    hsnow = HSi[j0,i0]        # m of snow over sea ice
  else:
    hsnow = HSi[j0,i0]*0.01   # m of snow over sea ice 

  # Distribute new snow depth evenly by cats in snow vol m3/m2:
  #ain = np.where(ain<1.e-20, 1.e-20, ain)  
  vsn_new = hsnow * ain        # m3/m2 per category
  if hsnow <= hs_min:
    vsn_new = vsn_new * 0.

  # Update snow enthalpy: J/m3  
  # see icepack_therm_vertical.F90 in icepack
  #
  # snow enthalpy should be: qsn_min <= qsn <= qsn_max
  # In theory, qsn_max = -rhos_Lfresh (latent heat of metling at 0C)
  # Make it a little lower to keep snow from melting right away
  #hsn_new = vsn_new / ain
  qsn = qsnon[:,j0,i0]        # enthalpy, J/kg < 0
  qsn_min = -rhos * Lfresh + (Tmin + 0.01) * cp_ice * rhos  # enth. of the coldest possible snow
  qsn_max = -rhos * Lfresh - 0.01 * cp_ice * rhos  # a little colder than 0C snow
  qT0 = -Lfresh*rhos      # enth. of pure snow at 0C

  # Update new enthalpy of new snow:
  # Clip to min/max enthalpy, set to 0 where no snow:
  qsn_new = np.clip(qsn, qsn_min, qsn_max)
  qsn_new = np.where(ain <= puny, 0., qsn_new)      # no ice
  qsn_new = np.where(vsn_new <= 0, 0., qsn_new)     # no snow

  vtot_init = np.nansum(vsn)
  vtot_new  = np.nansum(vsn_new)
  #print(f"tot vsnon change = {vtot_new-vtot_init}") 

  dvol_sum = dvol_sum + (vtot_new-vtot_init)
  vsnon_new[:,j0,i0] = vsn_new
  qsnon_new[:,j0,i0] = qsn_new

  diff = np.nansum(vsn_new - vsnon[:, j0, i0])
  diff2 = np.nansum(vsn_new -vsn)
  diff3 = np.nansum(vsnon[:,j0,i0] - vsnon_new[:,j0,i0])
  if diff == 0 and abs(diff2) > 0:
    print(f"No change at {j0},{i0}, expected diff={diff2}")
  #else:
  #  print(f"diff={diff}, diff2={diff2}")  

  #assert abs(diff) > 0, f"no change at {j0},{i0}, expected diff={diff2}"

  if diff3 == 0 and abs(diff2) > 0:
    print(f"No change in the arrays at {j0},{i0}, expected diff={diff2}")

# Checking:
print(f"dvol_sum = {dvol_sum}")
total_vsnon_init = np.nansum(vsnon)
total_vsnon_new  = np.nansum(vsnon_new)
print(f"Tot snow volume change (m3/m2): {total_vsnon_new - total_vsnon_init}")

# Update data set:
#ds_out = ds_out.assign(vsnon=vsnon_new, qsno001=qsnon_new)
ds_out['vsnon'].values[:] = vsnon_new
ds_out['qsno001'].values[:] = qsnon_new


# Sanity checking:
assert "vsnon" in ds_out and "qsno001" in ds_out, "Missing updated snow fields vsnon and qsno001"
assert ds_out["vsnon"].shape == vsnon_new.shape, "Check shape of vsnon "
assert ds_out["qsno001"].shape == qsnon_new.shape, "Check shape of qsnon "

# Attributes:
from datetime import datetime
istep1_val = ds_out.attrs.get('istep1', None)
ds_out.attrs.update({
    "title": "CICE6 restart with inserted hsnow from SSM/I NASA gridded fields for S. Ocean",
    "source": "insert_hsnow_cice6_restart.py",
    "contact": "dmitry.dukhovskoy@noaa.gov",
    "istep1": np.int32(istep1_val) if istep1_val is not None else np.int32(0), 
    "myear": np.int32(yrN),
    "mmonth": np.int32(mmN),
    "mday": np.int32(ddN),
    "msec": np.int32(nsecN),
    "history": f"Modified {datetime.now().isoformat()}",
})

# Save:flrst_in
sfx = extract_suffix(flrst_in)
if sfx is None:
  #flrst_out = f"cice_model.res.{yrN}{mmN:02d}{ddN:02d}.{nsecN:06d}.newhsnow.nc"
  flrst_out = f"cice_model.res.{yrN}{mmN:02d}{ddN:02d}.{hrN:02d}.hsnow.nc"
else:
  flrst_out = f"cice_model.res.{yrN}{mmN:02d}{ddN:02d}.{hrN:02d}.{sfx}.snow.nc"
dflrst_out = os.path.join(pthrest,flrst_out)
print(f"Saving CICE restart --> {dflrst_out}")
ds_out.to_netcdf(dflrst_out, encoding={var: {'_FillValue': None} for var in ds_out.data_vars}, format='NETCDF3_64BIT')
ds_out.close()




