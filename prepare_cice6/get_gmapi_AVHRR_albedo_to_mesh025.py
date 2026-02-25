"""
  Derive mapping indices gmapi to interpolate AVHRR high-res data --> mesh025 grid

  NOAA Climate Data Record (CDR) of AVHRR Polar Pathfinder Extended (APP-X) Cryosphere, Version 2

  NOAA Climate Data Record (CDR) of the eXtended AVHRR Polar Pathfinder (APP-X) 
  cryosphere contains 19 geophysical variables over the Arctic and Antarctic for the period 1982 - present. 

  https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C01580

Cite as: Key, Jeffrey; Wang, Xuanji; Liu, Yinghui; and NOAA CDR Program (2019). NOAA Climate Data Record of AVHRR Polar Pathfinder Extended (APP-X), Version 2. [indicate subset used]. NOAA National Centers for Environmental Information. doi:10.25921/AE96-0E57 [access date].

"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray
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

import mod_time as mtime
from mod_utils_fig import bottom_text
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_regmom as mrmom

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"Region to process",
                    choices=['north','south'], required=True, type=str)
args = parser.parse_args()

regn = args.regn if args.regn is not None else None

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

jdm, idm = hlon.shape


# Get grid data using OpenDAP
#Access dataset through OPeNDAP using the DAP2 protocol.
YR=2025
MM=1
DD=3
DDe = DD + 8
dnmb = mtime.datenum([YR,MM,DD])
jday = int(mtime.date2jday([YR,MM,DD]))
#url = f"https://www.star.nesdis.noaa.gov/thredds/dodsC/IceThickAVHRRnppSectorFourDayNP06/{YR}/"
if regn == 'north':
  url = f"https://www.ncei.noaa.gov/thredds/dodsC/avhrr-polar-pathfinder-ext-files/nhem/{YR}/"
  flname = f"Polar-APP-X_v02r00_Nhem_0400_d{YR}{MM:02d}{DD:02d}_c{YR}{MM:02d}{DDe:02d}.nc"
else:
  url = f"https://www.ncei.noaa.gov/thredds/dodsC/avhrr-polar-pathfinder-ext-files/shem/{YR}/"
  flname = f"Polar-APP-X_v02r00_Shem_0200_d{YR}{MM:02d}{DD:02d}_c{YR}{MM:02d}{DDe:02d}.nc"

dflinp = os.path.join(url, flname)

print(f"Reading {dflinp}")
with xarray.open_dataset(dflinp) as ds_ices:
  LON = ds_ices['longitude'].values
  LAT = ds_ices['latitude'].values

# Convert to -180, 180:
LON = (LON + 180.) % 360. - 180.


# Find lat bounds of AVHRR data:
# MOM6 points should be inside the AVHRR domain 
# to be able to find 4-vertice of the bounding box for interpolation
if regn == 'north':
  lat_min = np.min(LAT)
  lat_min = max([50.,lat_min])
  lat_max = 90. # override np.max(LAT) for Polar stereogr. projection, the code should work
                  # for correctly finding 4 points around hlat>np.max(LAT) but only
                  # for Polar stereogr. projection by grabbing points over the N. Pole

if regn == 'south':
lat_min = np.min(LAT)
lat_max = np.max(LAT)
lat_max = max([lat_max, -50])


ignore_north_lim = lat_max >= 90.  # make this true to avoid a discont. line across the Arctic Ocean
if ignore_north_lim:
  print(f"WARN: Indices north of northernmost lat={np.max(LAT):.4f} will be searched")

row_min = np.min(hlat, axis=1)
row_max = np.max(hlat, axis=1)
jS0 = np.argmax(row_min >= lat_min)
jE = len(row_max) - np.argmax(row_max[::-1] <= lat_max) - 1
#jE = np.argmin(row_max <= lat_max) - 1

def check_cntr(icc, ii, jj, iS0, jS0, iE0, jE0, ichk=200):
  if icc % ichk == 0:
    Npnts = (iE0 - iS0) * (jE0 - jS0)
    nIs_row = iE0 - iS0
    nJs_col = jE0 - jS0
    # Assuming loop is by i-indx first
    ncols_done = ii - iS0
    npnts_done = ncols_done * nJs_col + (jj - jS0 + 1)
    print(f' icc={icc} ii={ii} jj={jj} {(float(npnts_done)/Npnts*100.):.2f}% done ...')


INDX = None
JNDX = None
INDX_list = []
JNDX_list = []
IMOM = []
JMOM = []
icc = -1
iS = 0 
jS = jS0
iE = idm-1
Npnts = (iE+1) * ((jE+1) - jS)

print(f"Running gmapi i-indx={iS}:{iE}, j-indx={jS}:{jE}")

# Read saved tmp file if exist, and update indx arrays
#if use_tmp:
#  IMOM = combine_tmp('IMOM')
#  JMOM = combine_tmp('JMOM')
#  INDX = combine_tmp('INDX')
#  JNDX = combine_tmp('JNDX')
#  INDX_list = list(INDX)
#  JNDX_list = list(JNDX)
#  IMOM = list(IMOM)
#  JMOM = list(JMOM)
#  assert len(IMOM) == len(JMOM), f"IMOM and JMOM different lengths"
#  assert INDX.shape[0] == JNDX.shape[0], f"INDX and JNDX different lengths"

# Set of processed points for checking:
#pnts_done = set(zip(JMOM,IMOM))

Npnts = (iE + 1) * (jE + 1 - jS)
#print(f"TMP updated: Running gmapi i-indx={iS}:{iE}, j-indx={jS}:{jE}")

for ii in range(iS, iE+1):
  if ii >= idm:
    continue
  for jj in range(jS, jE+1):
    if jj >= jdm:
      continue

    #if (jj,ii) in pnts_done:
    #  icc += 1 
    #  check_cntr(icc, ii, jj, iS, jS, iE, jE, ichk=5000)
    #  continue

    if HH[jj,ii] >= 0:
      continue
    #ii, jj = mutil.find_indx_lonlat(358.047,89.816, hlon,hlat)
    x0 = hlon[jj,ii]
    y0 = hlat[jj,ii]
    if y0 < lat_min or y0 > lat_max:
      continue

    ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=1., ignore_north_lim=ignore_north_lim)
    ixx = np.asarray(ixx).reshape(-1)
    jxx = np.asarray(jxx).reshape(-1)

    check_cntr(icc, ii, jj, iS, jS, iE, jE, ichk=1000)
    if ixx.size == 0 or jxx.size == 0:
      continue

    icc += 1

    INDX_list.append(ixx.copy())
    JNDX_list.append(jxx.copy())
    IMOM.append(ii)
    JMOM.append(jj)

print("Main loop finished")

IMOM = np.array(IMOM)
JMOM = np.array(JMOM)
INDX = np.asarray(INDX_list)   # shape (N, 4)
JNDX = np.asarray(JNDX_list)

jdim, idim = LON.shape

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
darr_lon = xarray.DataArray(LON, dims=("jdim","idim"),\
                  coords={"jdim": np.arange(jdim),\
                          "idim": np.arange(idim)})
darr_lat = xarray.DataArray(LAT, dims=("jdim","idim"),\
                  coords={"jdim": np.arange(jdim),\
                          "idim": np.arange(idim)})

dset = xarray.Dataset({"mom_indx": darr_imom, \
                       "mom_jndx": darr_jmom, \
                       "gmapi_i": darr_indx,\
                       "gmapi_j": darr_jndx,\
                       "longit":  darr_lon,\
                       "latit":   darr_lat})

dset['mom_indx'].attrs['long_name'] = 'MOM6 grid I indices corresponding gmapi'
dset['mom_jndx'].attrs['long_name'] = 'MOM6 grid J indices corresponding gmapi'
dset['gmapi_i'].attrs['long_name'] = 'I indices AVHRR grid for interpolation'
dset['gmapi_j'].attrs['long_name'] = 'J indices AVHRR grid for interpolation'
dset['longit'].attrs['long_name']  = 'Longitudes AVHRR polar grid'
dset['latit'].attrs['long_name']   = 'Latitudes AVHRR polar grid'

# Global attributes
dset.attrs['title']       = 'Grid mapping between AVHRR on NSIDC EASE 25km Polar Stereographic to MOM6 mesh025 grids'
dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
dset.attrs['source']      = 'get_gmapi_AVHRR_albedo_to_mesh025.py'
dset.attrs['history']     = 'Arctic sea ice thickness from AVHRR 750m L3 product, NOAA/NESDIS/STAR/SOCD'
dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
dset.attrs['region']      = regn

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump = os.path.join(pthdata,'gmapi_mapping2mesh025')
fgmapi  = f'AVHRR_albedo_MOM6_gmapi_{idm}x{jdm}_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)

print(f'Saving gmapi --> {dfgmapi}')
dset.to_netcdf(dfgmapi, format='NETCDF4', engine='netcdf4')


f_chck = False
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.12, 0.12, 0.8, 0.8])

  # Plot on grid:
  #img = ax1.pcolormesh(C2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  # Check longitudes:

  JD, ID = zip(*pnts_done)

  ax1.cla()
  ax1.contour(HH, [0], linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
  ax1.plot(ID, JD, '.', markersize=1)


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

  btx = 'get_gmapi_AVHRR_albedo_to_mesh025.py'
  bottom_text(btx) 



