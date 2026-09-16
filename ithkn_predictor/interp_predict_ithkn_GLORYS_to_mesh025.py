"""
  Interpolate / remap 
  ML predicted ithkn from GLORYS to mesh025 grid
  see: predict_ML_ithkn_Ndays.py

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
import mod_icepredict as micepr
import mod_mom6 as mmom6


parser = argparse.ArgumentParser()
parser.add_argument("--regn",
      help="hemisphere: north or south",
      choices=['north', 'south'],
      required=True)
parser.add_argument("--rdate", help="Date of ithkn interpolation, YYYYMMDD", type=int, required=True)
parser.add_argument("--model", help="Model name",
                    choices=['clim','ols1','ols2','rf1','rf2','rf3','gbr1','gbr2','gbr3'],
                    required=True, type=str)
parser.add_argument("--iconc", help="Ice conc field used as a predictor",
                    choices=['glorys','amsr2','nsidc'],
                    type=str,
                    default="glorys")
parser.add_argument("--fsave", help=f"Save final dataset with all days as netcdf, default=1",
                    choices=[0,1],
                    default=1,
                    type=int)
parser.add_argument("--pcheck", help="=1: Plot to check interpolation, =0: no",
                    choices=[0,1],
                    default=0,
                    type=int)
args = parser.parse_args()

model   = args.model
regn    = args.regn if args.regn else None
rdate   = args.rdate
iconc_fld = args.iconc
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
dfgmapi  = os.path.join(pthgmapi, "gmapi_closenghb_GLORYS_to_UFSmesh025_north.nc")
with xr.open_dataset(dfgmapi) as ds:
  LONG  = ds["glorys_longit"].values
  LATG  = ds["glorys_latit"].values
  IGLR  = ds["glorys_indx"].values
  JGLR  = ds["glorys_jndx"].values
  IM025 = ds["mesh025_indx"].values
  JM025 = ds["mesh025_jndx"].values

LONG = (LONG + 360) % 360

# Read predicted ithkn:
# Load prediction and grid points:
# Training linregr params:
MODEL_NAMES = micepr.models_info()
model_name = MODEL_NAMES[model]

pthfcst = os.path.join(config_predictor["linregr"]["pthfcst"],f"{model_name}")
flfcst = f"{model_name}_ithkn_fcast_{rdate}.npz"
if not iconc_fld == "glorys":
  flfcst = f"{model_name}_ithkn_fcast_{iconc_fld}_{rdate}.npz"

#flfcst = f"{model_name}_{YS}_{YE}_ithkn_fcast_{rdate}.npz"
dflfcst = os.path.join(pthfcst, flfcst)
print(f"Loading fcst {dflfcst}")
data_fcst = np.load(dflfcst)
Ithkn = data_fcst['Yfcst']
JGF    = data_fcst['JG']
IGF    = data_fcst['IG']

# Get GLORYS grid
# Find file:
pthice = os.path.join(config_predictor["linregr"]["pthithkn"], f"{YR}")
dflice = mglr.find_file(rdate, pthice)
assert dflice is not None, f"GLORYS file not found for {rdate} in {pthice}"

with xr.open_dataset(dflice) as dsice:
  A2d = dsice['sithick'].isel(time=0).data.squeeze()
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values

hlon, hlat = np.meshgrid(LON, LAT)

# Replace glorys with predicted ithkn
AP = A2d * np.nan
AP[JGF,IGF] = Ithkn


# Create a look-up table (dictonary) for matching 
# every glorys indices JG,IG ---> mesh025 JM,IM
gmapi = {(jm,im):(jg,ig)
         for jm,im,jg,ig in zip(JM025, IM025, JGLR, IGLR)}

# Read mesh025 grid
#pthdata = '/archive/Dmitry.Dukhovskoy/data'
#pthice    = os.path.join(pthdata, 'ithkn_clim_combined')
#fliceout  = 'ithkn_mnthclim_cryo_avhrr_ices_1440x1080_north.nc'
#dfliceout = os.path.join(pthice,fliceout)

pthgrid = '/work/Dmitry.Dukhovskoy/GFSv17/mesh025_topo_grid'
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

LON_m25, LAT_m25 = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xr.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape
LMsk = np.where(HH<0, 1, 0)

# Ice thickness, mesh025
ITm25 = np.where(LMsk == 0, np.nan, 0)

DOMAIN = LMsk == 1
if regn == 'north':
  DOMAIN &= LAT_m25 > lat0
elif regn == 'south':
  DOMAIN &= LAT_m25 < lat0

JM, IM = np.where(DOMAIN  )

# Find GLORYS ----> mesh025 i, j pairs:
# Alternative to distance-approach, build a lookup table:
# every mesh025 JM,IM ---> glorys indices JG,IG
print("Finding mesh025 I, J to match GLORYS I,J")
JG = np.empty(len(JM), dtype=int)
IG = np.empty(len(IM), dtype=int)

for k, (jj, ii) in enumerate(zip(JM,IM)):
  if k > 0 and k % 10000 == 0:
    prc = k/len(JG)*100.
    print(f"  {prc:.2f}% processed")
  key = (int(jj), int(ii))
  if key not in gmapi:
    # mesh025 indices may be outside GLORYS subset region for ML 
    raise ValueError(f"Missing gmapi entry for GLORYS index {key}")
    continue
  JG[k], IG[k] = gmapi[key]

ITm25[JM,IM] = AP[JG,IG]

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
    
  dset = xr.Dataset({"ice_thkn": darr_cice})
  dset['ice_thkn'].attrs['long_name'] = 'ice thickness'
  dset["time"].attrs = {
       "long_name": f"days since {yr1}/01/01",
       "units" : "meters",
  } 
  
  # Add global attributes:
  dset.attrs['title']       = f'ML predicted ithkn, ML={model_name}, predictor iconc={iconc_fld}, '+\
                               'interpolated GLORYS to mesh025 grid'
  dset.attrs['institution'] = 'NOAA NWS OMD'
  dset.attrs['source']      = 'interp_predict_ithkn_GLORYS_to_mesh025.py'
  dset.attrs['region']      = regn

  print(f'Dumping interpolated ice thickness prediction --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')

if fsave:
  flice_out = f"ML_{model}_ithkn_iconc_{iconc_fld}_mesh025_{rdate}_{regn}.nc"
  dfliceout = os.path.join(pthfcst, flice_out)
  write_nc(dfliceout, dnmb, ITm25)

if plot_check:
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  clrmp.set_under(color=[1,1,1])

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.1, 0.13, 0.8, 0.8])


  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l', ax=ax1)
    # Subset region
    JJ = np.where(LAT_m25[:, 0] <= -50)[0]

  elif regn == 'north':
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l', ax=ax1)
    # Subset region
    JJ = np.where(LAT_m25[:, 0] >= 50)[0]

  hlat_s = LAT_m25[JJ, :]
  hlon_s = LON_m25[JJ, :]
  AP_s   = ITm25[JJ, :]

  xh, yh = m(hlon_s, hlat_s)

  if regn == 'north':
    m.drawparallels(np.arange(60, 90, 10), labels=[0,0,0,0])
  elif regn == 'south':
    m.drawparallels(np.arange(-80, -50, 10), labels=[0,0,0,0])

  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,0])
  m.drawcoastlines()

  img = ax1.pcolormesh(xh, yh, AP_s, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"ithkn {model} inp iconc={iconc_fld} inrtp to mesh025\n{YR}/{MM:02d}/{DD:02d}")

  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'interp_predict_ithkn_GLORYS_to_mesh025.py'
  bottom_text(btx)





