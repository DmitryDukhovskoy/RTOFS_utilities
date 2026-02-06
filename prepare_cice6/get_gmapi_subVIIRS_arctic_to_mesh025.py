"""
  Derive mapping indices gmapi to interpolate VIIRS high-res data --> mesh025 grid

  This script runs for a subset of all points to speed up the process,
  saves those in tmp files (npy) to be combined later 
  Run it for iS - iE, specifying the start i indx and N i-pnts to process

  VIIRS ice thickness
  on native grid   
  daily and 4-day composites 
                    
  ice and snow depth, Arctic

https://coastwatch.noaa.gov/cwn/products/viirs-sea-ice-concentration-ice-thickness-ice-surface-temperature.html

Use NSIDC Polar Stereogr Proj:
semi_major_axis: 6378137.0
inverse_flattening: 298.257223563
straight_vertical_longitude_from_pole: -45.0
latitude_of_projection_origin: 90.0
standard_parallel: 70.0
false_easting: 0.0
false_northing: 0.0
srid: EPSG:3413
proj4text: +proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
crs_wkt: PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3413"]]

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

import mod_time as mtime
from mod_utils_fig import bottom_text
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

regn = 'north'
iS0 = 0
iE0 = None
parser = argparse.ArgumentParser()
parser.add_argument("--iS0", help=f"For running N jobs, specify i start index from 0,...,max(idim), otherwise={iS0}",
                            required=True, type=int)
parser.add_argument("--iE0", help=f"If running N jobs, specify i end index to process: > iS0 < idm",
                             required=True, type=int)
args = parser.parse_args()

iS0  = args.iS0 if args.iS0 is not None else None
iE0  = args.iE0 if args.iE0 is not None else None

print("====  WARNING: Saving temporary files with gmapi only ===")

def check_cntr(icc, ii, jj, iS0, jS0, iE0, jE0, ichk=200):
  if icc % ichk == 0:
    Npnts = (iE0 - iS0) * (jE0 - jS0)
    nIs_row = iE0 - iS0
    nJs_col = jE0 - jS0
    # Assuming loop is by i-indx first
    ncols_done = ii - iS0
    npnts_done = ncols_done * nJs_col + (jj - jS0 + 1)
    print(f' icc={icc} ii={ii} jj={jj} {(float(npnts_done)/Npnts*100.):.2f}% done ...')

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
YR=2021
MM=9
DD=15
dnmb = mtime.datenum([YR,MM,DD])
jday = int(mtime.date2jday([YR,MM,DD]))
#url = f"https://www.star.nesdis.noaa.gov/thredds/dodsC/CoastWatch/VIIRS/npp/IceThick/FourDaySector/NP06"
url = f"https://www.star.nesdis.noaa.gov/thredds/dodsC/IceThickVIIRSnppSectorFourDayNP06/{YR}/"
flviirs = f"VXSACW_B{YR}{jday-3:03d}_B{YR}{jday:03d}_H4_NP06_edgemask_IceThickness.nc"
dflviirs = os.path.join(url, flviirs)

print(f"Reading {dflviirs}")
with xarray.open_dataset(dflviirs) as ds_ices:
  LON = ds_ices['longitude'].values
  LAT = ds_ices['latitude'].values

# Convert to -180, 180:
LON = (LON + 180.) % 360. - 180.


import mod_regmom as mrmom

# Find lat bounds of CryoSat data:
# MOM6 points should be inside the CryoSat domain 
# to be able to find 4-vertice of the bounding box for interpolation
lat_min = np.min(LAT)
lat_min = max([50.,lat_min])
lat_max = 90. # override np.max(LAT) for Polar stereogr. projection, the code should work
                # for correctly finding 4 points around hlat>np.max(LAT) but only
                # for Polar stereogr. projection by grabbing points over the N. Pole

ignore_north_lim = lat_max >= 90.  # make this true to avoid a discont. line across the Arctic Ocean
if ignore_north_lim:
  print(f"WARN: Indices north of northernmost lat={np.max(LAT):.4f} will be searched")

row_min = np.min(hlat, axis=1)
row_max = np.max(hlat, axis=1)
jS0 = np.argmax(row_min >= lat_min)
jE = len(row_max) - np.argmax(row_max[::-1] <= lat_max) - 1
#jE = np.argmin(row_max <= lat_max) - 1

INDX = None
JNDX = None
IMOM = []
JMOM = []
INDX_list = []
JNDX_list = []
init = False
icc = -1
iS = iS0
jS = jS0
iE = iE0
iE = min(iE, idm)
Npnts = (iE - iS) * (jE - jS)

print(f"Running gmapi i-indx={iS}:{iE}, j-indx={jS}:{jE}")

# TMP files:
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump = os.path.join(pthdata,'gmapi_NSIDC','tmp')
os.makedirs(pthdump, exist_ok=True)

fltmp_im = f"VIIRS_gmapi_tmp_IMOM_iS{iS:04d}_iE{iE:04d}.npy"
fltmp_jm = f"VIIRS_gmapi_tmp_JMOM_iS{iS:04d}_iE{iE:04d}.npy"
fltmp_ix = f"VIIRS_gmapi_tmp_INDX_iS{iS:04d}_iE{iE:04d}.npy"
fltmp_jx = f"VIIRS_gmapi_tmp_JNDX_iS{iS:04d}_iE{iE:04d}.npy"
dflim = os.path.join(pthdump, fltmp_im)
dfljm = os.path.join(pthdump, fltmp_jm)
dflix = os.path.join(pthdump, fltmp_ix)
dfljx = os.path.join(pthdump, fltmp_jx)

# Load saved tmp file if exist, and update counter
if not (os.path.exists(dflim) and os.path.exists(dfljm) and
      os.path.exists(dflix) and os.path.exists(dfljx)):    
  print(f"Tmp files *_iS{iS:04d}_iE{iE:04d}.npy do not exist, start from iS={iS}")
else:
  IMOM = np.load(dflim)
  JMOM = np.load(dfljm)
  INDX = np.load(dflix)
  JNDX = np.load(dfljx)

  IMOM = IMOM.tolist()
  JMOM = JMOM.tolist()
  INDX_list = list(INDX)
  JNDX_list = list(JNDX)

  init = True
  icc = len(IMOM)-1
  iS = IMOM[-1]
  jS = JMOM[-1]+1 
  print(f"Starting from TMP: Updating iS={iS} and jS={jS}")

assert iS >= iS0, f"Check iS={iS} starting iS0={iS0}"    

print(f"TMP updated: Running gmapi i-indx={iS}:{iE}, j-indx={jS}:{jE}")

for ii in range(iS, iE+1):
  if ii == idm:
    continue
  pass1 = (init and ii == iS)
  if pass1:
    njj = jS - jS0
    nii = ii - iS0 + 1
    print(f'Starting from last record:  icc={icc}')
    check_cntr(icc, ii, jS, iS0, jS0, iE, jE, ichk=icc)

    init = False
    if jS > jE:
      continue   
  else:
    jS = jS0

  for jj in range(jS, jE+1):
    if jj == jdm:
      continue

    if HH[jj,ii] >= 0:
      continue
    #ii, jj = mutil.find_indx_lonlat(358.047,89.816, hlon,hlat)
    x0 = hlon[jj,ii]
    y0 = hlat[jj,ii]
    if y0 < lat_min or y0 > lat_max:
      continue

    ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=0.1, ignore_north_lim=ignore_north_lim)
    ixx = np.asarray(ixx).reshape(-1)
    jxx = np.asarray(jxx).reshape(-1)
    if ixx.size == 0 or jxx.size == 0:
      continue

    icc += 1

    INDX_list.append(ixx.copy())
    JNDX_list.append(jxx.copy())
    IMOM.append(ii)
    JMOM.append(jj)

    check_cntr(icc, ii, jj, iS0, jS0, iE, jE, ichk=200)
    if icc > 0 and icc % 200 == 0:
      print(f"Dumping tmp files IMOM, JMOM, INDX, JNDX --> {pthdump} npy files, iS0={iS0} iE0={iE0}")
      INDX = np.asarray(INDX_list)   # shape (N, 4)
      JNDX = np.asarray(JNDX_list)
      np.save(dflim, np.array(IMOM))
      np.save(dfljm, np.array(JMOM))
      np.save(dflix, INDX)
      np.save(dfljx, JNDX)

# Dump at the end:
print(f'Finished iS0={iS0} : iE0={iE0}')
print(f"Dumping IMOM --> {dflim}")
print(f"Dumping JMOM --> {dfljm}")
print(f"Dumping INDX --> {dflix}")
print(f"Dumping JNDX --> {dfljx}")

INDX = np.asarray(INDX_list)   # shape (N, 4)
JNDX = np.asarray(JNDX_list)
np.save(dflim, np.array(IMOM))
np.save(dfljm, np.array(JMOM))
np.save(dflix, INDX)
np.save(dfljx, JNDX)


f_chck = False
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:
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
  #img = ax1.pcolormesh(C2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  # Check longitudes:
  ax1.contour(HH, [0], linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)


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





