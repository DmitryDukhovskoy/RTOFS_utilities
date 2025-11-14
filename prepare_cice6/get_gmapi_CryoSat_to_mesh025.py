"""
  Gridded estimates of Antarctic sea ice physical properties derived from 
  CryoSat-2 Baseline-D SAR and SARIn data spanning July 2010 through August 2021. 
  Data are generated using the CryoSat-2 Waveform-Fitting method for Antarctic sea ice (CS2WFA).

  Fons, S., Kurtz, N., & Bagnardi, M. (2022). 
  Antarctic Sea Ice Thickness Estimates from CryoSat-2: 2010-2021 (0.1.1) [Data set]. 
  Zenodo. https://doi.org/10.5281/zenodo.7327711

  Grid seems to be similar to 
  NASA SSMI & AMSR snow thickness on sea ice in Antarctic
  and to NSIDC (Polar Sterographic)

  Derive gmapi indices for bi-linear interpolation

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
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6

YR = 2020
MM = 1
regn = 'south'

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--yr", help=f"year of CryoSat data, default={YR}", type=int)
parser.add_argument("--mm", help=f"month of CryoSat data, default={MM}", type=int)
args = parser.parse_args()

regn = args.regn if args.regn else regn
YR   = args.yr if args.yr else YR
MM   = args.mm if args.mm else MM

syst_info = os.uname()
machine = syst_info.nodename

if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
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

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

# Lon/lats have been derived in the code
# that computed monthly clim for NASA hsnow
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthice  = os.path.join(pthdata, 'CryoSat2_antarctic_ice_snow_thkn')
flice   = f"CS2WFA_25km_{YR}{MM:02d}.nc"
dflice  = os.path.join(pthice,flice)
print(f"Loading {dflice}")
# Note the NetCDF files have subgroups inside
# Open the subgroup 'sea_ice_thickness' to read ice thikn.
# other variables are in the root group:
with xarray.open_dataset(dflice) as dshi:
  LON = dshi['lon'].data.squeeze()
  LAT = dshi['lat'].data.squeeze()

LON = np.where(LON > 180., LON-360., LON)

# Similar grid, snow data
pthsnow = os.path.join(pthdata, 'snow_nasa')
flhs = 'AMSR_Antarctic_hsnow_month_clim_1998_2007.nc'
dflhs = os.path.join(pthsnow,flhs)
print(f"Loading {dflhs}")
with xarray.open_dataset(dflhs) as dshs:
  XX = dshs['xpolar'].data.squeeze()
  YY = dshs['ypolar'].data.squeeze()
  LONS = dshs['lon'].data.squeeze()
  LATS = dshs['lat'].data.squeeze()

import mod_regmom as mrmom

# Find lat bounds of CryoSat data:
# MOM6 points should be inside the CryoSat domain 
# to be able to find 4-vertice of the bounding box for interpolation
lat_min = np.min(LAT)
lat_max = np.max(LAT)

if regn == 'south':
  lat_max = -50.
else:
  lat_min = 50.

row_min = np.min(hlat, axis=1)
row_max = np.max(hlat, axis=1)
jS = np.argmax(row_min >= lat_min)
jE = len(row_max) - np.argmax(row_max[::-1] <= lat_max) - 1
#jE = np.argmin(row_max <= lat_max) - 1

jdm, idm = hlon.shape
icc = -1
INDX = None
JNDX = None
IMOM = []
JMOM = []
for ii in range(idm):
  if ii%50 == 0:
    print(f' icc={icc} {ii/idm*100:.2f}% done ...')
  for jj in range(jS,jE):
    if HH[jj,ii] >= 0:
      continue
    x0 = hlon[jj,ii]
    y0 = hlat[jj,ii]
    if y0 < lat_min or y0 > lat_max:
      continue

    icc += 1
    ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=1.)
    if len(ixx)==0 or len(jxx)==0:
     continue
    ixx = np.expand_dims(ixx, axis=0)
    jxx = np.expand_dims(jxx, axis=0)

    if icc == 0:
      INDX = ixx.copy()
      JNDX = jxx.copy()
    else:
      INDX = np.append(INDX, ixx, axis=0)
      JNDX = np.append(JNDX, jxx, axis=0)

    IMOM.append(ii)
    JMOM.append(jj)

IMOM = np.array(IMOM)
JMOM = np.array(JMOM)

npnts = len(IMOM)
darr_imom = xarray.DataArray(IMOM, dims=("npoints"),\
                   coords={"npoints": np.arange(npnts)})
darr_jmom = xarray.DataArray(JMOM, dims=("npoints"),\
                   coords={"npoints": np.arange(npnts)})
darr_indx = xarray.DataArray(INDX, dims=("npoints","nvert"),\
                   coords={"npoints": np.arange(npnts),\
                           "nvert": np.arange(4)})
darr_jndx = xarray.DataArray(JNDX, dims=("npoints","nvert"),\
                   coords={"npoints": np.arange(npnts),\
                           "nvert": np.arange(4)})
dset = xarray.Dataset({"mom_indx": darr_imom, \
                       "mom_jndx": darr_jmom, \
                       "gmapi_i": darr_indx,\
                       "gmapi_j": darr_jndx})

dset['mom_indx'].attrs['long_name'] = 'MOM6 grid I indices corresponding gmapi'
dset['mom_jndx'].attrs['long_name'] = 'MOM6 grid J indices corresponding gmapi'
dset['gmapi_i'].attrs['long_name'] = 'I indices NSIDC grid for interpolation'
dset['gmapi_j'].attrs['long_name'] = 'J indices NSIDC grid for interpolation'

# Global attributes
dset.attrs['title']       = 'Grid mapping between Polar sterographic south and MOM6 grids'
dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
dset.attrs['source']      = 'get_gmapi_CryoSat_to_mesh025.py'
dset.attrs['history']     = 'CryoSat snow and ice thickness cell mean, polar stereogr., monthly'
dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
dset.attrs['region']      = regn

pthdump = os.path.join(pthdata,'gmapi_NSIDC')
fgmapi  = f'CryoSat_MOM6_gmapi_{idm}x{jdm}_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)

print(f'Saving gmapi --> {dfgmapi}')
dset.to_netcdf(dfgmapi, format='NETCDF4', engine='netcdf4')


f_chck = False
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:
  # Note the NetCDF files have subgroups inside
  # Open the subgroup 'sea_ice_thickness' to read ice thikn.
  # other variables are in the root group:
  with xarray.open_dataset(dflice, group='sea_ice_thickness') as dshi:
    C2d = dshi['sea_ice_thickness'].data.squeeze()

  #clrmp = mclrmps.colormap_uv()
  #rmin = -1.
  #rmax = 1.

  #clrmp = mclrmps.colormap_conc()
  #rmin = 0.
  #rmax = 1.

  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  LON1 = LON.copy()
  LON2 = LON.copy()
  LON1 = np.where(LON1 < -175, np.nan, LON1)
  LON2 = np.where(LON2 > 172, np.nan, LON2)
  LON3 = np.where(LON < 0, LON+360., LON)
  LON3 = np.where(LON3 > 350., np.nan, LON3)
  lon_cntr1 = [x for x in range(-170,0,10)]  # grey -180:0
  lon_cntr2 = [x for x in range(10,178,10)]  # blue: 0 to 180 E
  lat_cntr = [x for x in range(-85,-20,5)]

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.12, 0.12, 0.8, 0.8])

  # Plot on grid:
  img = ax1.pcolormesh(C2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  # Check longitudes:
  cs = ax1.contour(LON1, lon_cntr1, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
  ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
  cs2 = ax1.contour(LON2, lon_cntr2, linestyles='solid', colors=[(0.,0.5,0.9)], linewidths=1)
  ax1.clabel(cs2, inline=True, fontsize=10, fmt="%.1f")
  cs3 = ax1.contour(LON1,[0], linestyles='solid', colors=[(1,0.,0.)], linewidths=2)
  ax1.clabel(cs3, inline=True, fontsize=12, fmt="%.1f")
  ax1.contour(LON3,[180], linestyles='solid', colors=[(0.6,0.,0.9)], linewidths=1)
  ax1.contour(LAT, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
  ax1.contour(LAT,[-75], linestyles='solid', colors=[(1,0.,0.)], linewidths=1)  
  ax1.axis('scaled')
  ax1.invert_yaxis() 
  ax1.set_ylabel('Inverted j index')
  ax1.set_ylabel('j index')
  ax1.set_xlabel('i index')

  ax1.set_title(f'CryoSat ice thickness {YR}/{MM:02d}')

  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  if rmin < 0:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  else:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')

  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'get_gmapi_CryoSat_to_mesh025.py'
  bottom_text(btx) 





