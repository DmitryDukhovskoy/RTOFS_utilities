"""
  Interpolate / remap 
  GLORYS daily ice thickness to UFS mesh025
  saved by months

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
parser.add_argument("--sdate", help="Start Date of ithkn interpolation, YYYYMMDD", type=int, required=True)
parser.add_argument("--edate", help="End Date of ithkn interpolation, YYYYMMDD", type=int, required=True)
parser.add_argument("--pcheck", help="=1: Plot to check interpolation, =0: no",
                    choices=[0,1],
                    default=0,
                    type=int)
args = parser.parse_args()

sdate   = args.sdate
edate   = args.edate
regn    = args.regn
plot_check = args.pcheck == 1

fsave = True

# Dates:
dnmbS = int(mtime.rdate2datenum(sdate))
YRs, MMs, DDs = mtime.datevec(dnmbS)[:3]
dnmbE = int(mtime.rdate2datenum(edate))
YRe, MMe, DDe = mtime.datevec(dnmbE)[:3]

regions = {
    "north": ("Arctic", 60.0),
    "south": ("Antarctic", -55.0),
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


# Get GLORYS grid
# Find file:
pthice = os.path.join(config_predictor["linregr"]["pthithkn"], f"{YRs}")
dflice = mglr.find_file(sdate, pthice)
assert dflice is not None, f"GLORYS file not found for {sdate} in {pthice}"

with xr.open_dataset(dflice) as dsice:
  LON = dsice['longitude'].values
  LAT = dsice['latitude'].values

glon, glat = np.meshgrid(LON, LAT)


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


def write_nc(dfliceout, time_dnmb, A3d):
  yr1, mm1, dd1 = mtime.datevec(time_dnmb[0])[:3]
  nrecs, jdim, idim = A3d.shape
 
  # Days wrt to the Jan 1st:
  time_days = time_dnmb - mtime.datenum([yr1, 1, 1])
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
  dset.attrs['title']       = 'GLORYS field interpolated GLORYS to mesh025 grid, closest neighbour'
  dset.attrs['institution'] = 'NOAA NWS OMD'
  dset.attrs['source']      = 'interp_GLORYSithkn_to_mesh025_month.py'
  dset.attrs['region']      = regn

  print(f'Dumping GLORYS interpolated ice thickness  --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')

def read_glorys(rdate):
  dnmbR = mtime.rdate2datenum(rdate*100)  # restart day nmb
  YR, MM, DD = mtime.datevec(dnmbR)[:3]
  pthice = f"/uda/Global_Ocean_Physics_Reanalysis/global/daily/sithick/{YR}"
  dflice = mglr.find_file(rdate, pthice)
  if not os.path.isfile(dflice):
    raise RuntimeError(f"File not found: {dflice}")

  with xr.open_dataset(dflice) as dsice:
    A2d = dsice['sithick'].isel(time=0).data.squeeze()
    #LON = dsice['longitude'].values
    #LAT = dsice['latitude'].values

  A2d[np.isnan(A2d)] = 0.

  return A2d


# Save by months
MMwrk = 0
YRwrk = 0
irec = 0 
A3d  = []
time_dnmb = []
for dnmb0 in range(dnmbS, dnmbE+1):
  YR, MM, DD = mtime.datevec(dnmb0)[:3]

  if not MM == MMwrk:
    # Save previous month fields:
    if irec > 0:
      # Convert lists --> np array AAin etc
      # Save processed fields    
      if not fsave:
        print(f"Final netcdf is not saved, save_nc={save_nc}")
      else:
        flice_out = f"GLORYS_ithkn_interp_mesh025_{YRwrk}{MMwrk:02d}_{regn}.nc"
        pthintrp = '/work/Dmitry.Dukhovskoy/data/GLORYS_ithkn_interp_UFSmesh025' 
        dfliceout = os.path.join(pthintrp, flice_out)

        A3d = np.asarray(A3d)
        time_dnmb = np.asarray(time_dnmb)

        write_nc(dfliceout, time_dnmb, A3d)

    # Prepare fields for processing new month: 
    MMwrk = MM
    YRwrk = YR
    irec = 0
    A3d = []
    time_dnmb = []

  rdate = YR*10000 + MM*100 + DD
  ITglr = read_glorys(rdate)
  ITm25 = np.where(LMsk == 0, np.nan, 0)
  ITm25[JM,IM] = ITglr[JG,IG]

  irec += 1
  time_dnmb.append(dnmb0)
  A3d.append(ITm25)

if fsave:
  A3d = np.asarray(A3d)
  time_dnmb = np.asarray(time_dnmb)
  flice_out = f"GLORYS_ithkn_interp_mesh025_{YR}{MM:02d}_{regn}.nc"
  pthintrp = '/work/Dmitry.Dukhovskoy/data/GLORYS_ithkn_interp_UFSmesh025'
  dfliceout = os.path.join(pthintrp, flice_out)
  write_nc(dfliceout, time_dnmb, A3d)

if plot_check:
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  clrmp.set_under(color=[1,1,1])

  iday = 10
  A2d = A3d[iday-1,:,:].squeeze()

  print(f"Plotting day={iday}")
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

  glat_s = LAT_m25[JJ, :]
  glon_s = LON_m25[JJ, :]
  AP_s   = A2d[JJ, :]

  xh, yh = m(glon_s, glat_s)

  if regn == 'north':
    m.drawparallels(np.arange(60, 90, 10), labels=[0,0,0,0])
  elif regn == 'south':
    m.drawparallels(np.arange(-80, -50, 10), labels=[0,0,0,0])

  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,0])
  m.drawcoastlines()

  img = ax1.pcolormesh(xh, yh, AP_s, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f"ithkn GLORYS inrtp to mesh025 {YR}/{MM:02d}/{iday:02d}")

  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'interp_GLORYSithkn_to_mesh025_month.py'
  bottom_text(btx)





