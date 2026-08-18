"""
  Interpolate / remap AMSR2 sea ice conc fields
  from mesh025 grid --> GLORYS for ML predictions

  AMSR2 fields processed on gaea:
  interp from AMSR2 grid --> mesh025:
  interp_AMSR2_iconc_mesh025.py

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray as xr
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, cm
from yaml import safe_load
import argparse
#from pathlib import Path


#ROOT = Path(__file__).resolve().parent

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
import mod_glorys as mglr
import mod_colormaps as mclrmps


parser = argparse.ArgumentParser()
parser.add_argument("--regn",
      help="hemisphere: north or south",
      choices=['north', 'south'],
      required=True)
parser.add_argument("--rdate", help="Date of AMSR2 interpolation, YYYYMMDD", type=int, required=True)
parser.add_argument("--tmpf", help="=1: Field saved in tmp numpy file, =0: in monthly netcdf",
                    choices=[0,1],
                    required=True,
                    type=int)
parser.add_argument("--fsave", help=f"Save final dataset with all days as netcdf, default=1",
                    choices=[0,1],
                    default=1,
                    type=int)
parser.add_argument("--pcheck", help="=1: Plot to check interpolation, =0: no",
                    choices=[0,1],
                    default=0,
                    type=int)
args = parser.parse_args()

regn    = args.regn if args.regn else None
rdate   = args.rdate
use_tmp = args.tmpf == 1
fsave   = args.fsave
plot_check = args.pcheck == 1

dnmb = mtime.rdate2datenum(rdate)
YR, MM, DD = mtime.datevec(dnmb)[:3]

regions = {
    "north": ("Arctic", 65.0),
    "south": ("Antarctic", -60.0),
}
regn_name, lat0 = regions[regn]


fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

# Load gmapi:
pthgmapi = config_predictor["linregr"]["pthgmapi"]
dfgmapi  = os.path.join(pthgmapi, "gmapi_closenghb_mesh025_to_GLORYS_north.nc")
with xr.open_dataset(dfgmapi) as ds:
  LONG  = ds["glorys_longit"].values
  LATG  = ds["glorys_latit"].values
  IGLR  = ds["glorys_indx"].values
  JGLR  = ds["glorys_jndx"].values
  IM025 = ds["mesh025_indx"].values
  JM025 = ds["mesh025_jndx"].values

LONG = (LONG + 360) % 360


# Find interp points on GLORYS grid:
# GLORYS grid:
# Read GLORYS grid:
pthithkn = config_predictor["linregr"]["pthithkn"]
pthice = os.path.join(pthithkn,f"{YR}")

# Find file:
dflglr = mglr.find_file(rdate, pthice)
assert dflglr is not None, f"GLORYS file not found for {sdate} in {pthice}"

with xr.open_dataset(dflglr) as dsice:
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values

# 0 <= lon < 360
LON = (LON + 360) % 360
hlon, hlat = np.meshgrid(LON, LAT)

# Land mask:
LMsk = None
pthssh = os.path.join(config_predictor["linregr"]["pthssh"],f"{YR}")
dflssh = mglr.find_file(rdate, pthssh)
with xr.open_dataset(dflssh) as dszos:
  SSH = dszos['zos'].isel(time=0).values.squeeze()

LMsk = np.where(np.isfinite(SSH),1,0)

# Points where ice thickness is predicted:
# Define domain:
if regn == 'north':
  DOMAIN = (hlat > lat0) & (LMsk == 1)
elif regn == 'south':
  DOMAIN = (hlat < lat0) & (LMsk == 1)

JG, IG = np.where(DOMAIN)


# Create a look-up table (dictonary) for matching 
# every glorys indices JG,IG ---> mesh025 JM,IM
gmapi = {(jg,ig):(jm,im)
         for jg,ig,jm,im in zip(JGLR, IGLR, JM025, IM025)}

# Read mesh025 grid
pthdata = '/archive/Dmitry.Dukhovskoy/data'
pthice    = os.path.join(pthdata, 'ithkn_clim_combined')
fliceout  = 'ithkn_mnthclim_cryo_avhrr_ices_1440x1080_north.nc'
dfliceout = os.path.join(pthice,fliceout)

with xr.open_dataset(dfliceout) as dsice:
  LONM025 = dsice['lon'].data
  LATM025 = dsice['lat'].data


pthamsr = '/work/Dmitry.Dukhovskoy/data/AMSR2_iconc_interp'
if use_tmp:
  flnm = f"AMSR2_iconc_mesh025_{YR}{MM:02d}{DD:02d}_{regn}.npy"
  dflnm = os.path.join(pthamsr, flnm)
  print(f"Loading tmp file: {dflnm}")
  IC_m25 = np.load(dflnm)

  
# Find mesh025 --> GLORYS i, j pairs:
# Alternative to distance-approach, build a lookup table:
# every mesh025 JM,IM ---> glorys indices JG,IG
print("Finding GLORYS I, J to match mesh025 I,J")
JM = np.empty(len(JG), dtype=int)
IM = np.empty(len(IG), dtype=int)

for k, (jj, ii) in enumerate(zip(JG, IG)):
  if k > 0 and k % 10000 == 0:
    prc = k/len(JG)*100.
    print(f"  {prc:.2f}% processed")
  key = (int(jj), int(ii))
  if key not in gmapi:
    # mesh025 indices may be outside GLORYS subset region for ML 
    #raise ValueError(f"Missing gmapi entry for GLORYS index {key}")
    continue
  JM[k], IM[k] = gmapi[key]

ICgl = np.where(LMsk == 0, np.nan, 0)
ICgl[JG,IG] = IC_m25[JM,IM]

def write_nc(dfliceout, time_dnmb, A2d):
  yr1, mm1, dd1 = mtime.datevec(dnmb)[:3]
  nrecs = 1
  jdim, idim = A2d.shape
  A3d = np.expand_dims(A2d, axis=0) 
 
  # Days wrt to the Jan 1st:
  time_days = np.array([dnmb - mtime.datenum([yr1, 1, 1])])
  darr_cice = xr.DataArray(A3d, dims=("time","jdim","idim"),\
                     coords={"time": time_days,\
                             "jdim": np.arange(jdim),\
                             "idim": np.arange(idim)})
    
  dset = xr.Dataset({"ice_conc": darr_cice})
  dset['ice_conc'].attrs['long_name'] = 'ice partial area'
  dset["time"].attrs = {
       "long_name": f"days since {yr1}/01/01"
  } 
  
  # Add global attributes:
  dset.attrs['title']       = 'AMSR2 L4 OSI SAF EUMETSAT sea ice concentration daily interpolated onto mesh025 grid and then onto GLORYS grid'
  dset.attrs['institution'] = 'NOAA NWS OMD'
  dset.attrs['source']      = 'interp_iconc_mesh025_to_GLORYS.py'
  dset.attrs['region']      = regn

  print(f'Dumping interpolated ice conc --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')

if fsave:
  dfliceout = os.path.join(pthamsr, f"AMSR2_iconc_GLORYSgrid_{rdate}_north.nc")
  write_nc(dfliceout, dnmb, ICgl)

if plot_check:
  clrmp = mclrmps.colormap_conc()
  rmin = 0.
  rmax = 1.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  clrmp.set_under(color=[1,1,1])

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.1, 0.13, 0.8, 0.8])


  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l', ax=ax1)
    # Subset region
    JJ = np.where(hlat[:, 0] <= -50)[0]

  elif regn == 'north':
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l', ax=ax1)
    # Subset region
    JJ = np.where(hlat[:, 0] >= 50)[0]

  hlat_s = hlat[JJ, :]
  hlon_s = hlon[JJ, :]
  AP_s   = ICgl[JJ, :]

  xh, yh = m(hlon_s, hlat_s)

  if regn == 'north':
    m.drawparallels(np.arange(60, 90, 10), labels=[0,0,0,0])
  elif regn == 'south':
    m.drawparallels(np.arange(-80, -50, 10), labels=[0,0,0,0])

  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,0])
  m.drawcoastlines()

  img = ax1.pcolormesh(xh, yh, AP_s, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"AMSR2 OSI SAF iconc interp to GLORYS grid \n{YR}/{MM:02d}/{DD:02d}")

  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'interp_iconc_mesh025_to_GLORYS.py'
  bottom_text(btx)





