"""
  Monthly snow climatology from EWG snow data
  Data prepared from S. Union field observations
  Warren climatology (?)

  The Arctic Meteorology and Climate Atlas

  From EWG 
  https://nsidc.org/data/search#keywords=Arctic+snow+climatology/sortKeys=score,,desc/facetFilters=%257B%257D/pageNumber=1/itemsPerPage=25

  Grid info - see Documentation directory / technical_documentation.pdf

All of the gridded fields on the Atlas are in EASE-Grid format. EASE-Grid is a set of
equal-area projections and grids developed at the National Snow and Ice Data Center
(EWG) to be a tool for users of global-scale gridded data. The Atlas gridded fields are
in an azimuthal equal-area projection centered on the North Pole. The grid cell size is
250 km. There are a total of 529 cells in a 23- by 23- cell array, with cells numbered 0 to
22. Cell 11,11 is centered on the North Pole.

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
import pandas as pd

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

regn = 'north'

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

# EWG Atlas grid
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthgrid = os.path.join(pthdata,'Warren_snow_clim_EWG_atlas/DATA/GRIDDED_FIELDS/EASE_INFO')
dflgrid = os.path.join(pthgrid,'N65-250km.latlon.dat')
pthdump  = os.path.join(pthdata,"gmapi_NSIDC")

# Read the ASCII file
dasc = pd.read_csv(
       dflgrid,
       sep=r"\s+",  
       comment="#",
       engine="python",
)

#print(dasc)
lat = dasc["Lat"].to_numpy()
lon = dasc["Lon"].to_numpy()
JDX = dasc["Row"].to_numpy(dtype=int)
IDX = dasc["Col"].to_numpy(dtype=int)

jdim = np.max(JDX) + 1
idim = np.max(IDX) + 1
LAT = np.zeros((jdim, idim)) * np.nan
LON = np.zeros((jdim, idim)) * np.nan
assert JDX.min() == 0 and IDX.min() == 0

LAT[JDX,IDX] = lat
LON[JDX,IDX] = lon
assert LAT.shape == (jdim, idim)

#pthnsidc = os.path.join(pthdata,f"NRT_NOAA_EWG_seaconc/{YR}")


import mod_misc1 as mmisc
import mod_regmom as mrmom

# Find lat bounds of EWG data:
# MOM6 points should be inside the EWG domain 
# to be able to find 4-vertice of the bounding box for interpolation
lat_min = 60.
lat_max = np.max(LAT)
ignore_north_lim = lat_max >= 90.
if ignore_north_lim:
  print(f"WARN: Indices north of northernmost lat={np.max(LAT):.4f} will be searched")

row_min = np.min(hlat, axis=1)  # min lat in each row
row_max = np.max(hlat, axis=1)  # max lat in each row
jS = np.argmax(row_min >= lat_min)
jE = len(row_max) - np.argmax(row_max[::-1] <= lat_max) - 1 # note reverse indexing for [::-1]
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
  for jj in range(jS,jE+1):
    if HH[jj,ii] >= 0:
      continue
    x0 = hlon[jj,ii]
    y0 = hlat[jj,ii]
    if y0 < lat_min or y0 > lat_max:
      continue

    icc += 1
    ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=3., ignore_north_lim=ignore_north_lim)
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

# EWG coord dimensions:
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
dset['gmapi_i'].attrs['long_name'] = 'I indices EWG grid for interpolation'
dset['gmapi_j'].attrs['long_name'] = 'J indices EWG grid for interpolation'
dset['longit'].attrs['long_name']  = 'Longitudes derived from EWG polar grid'
dset['latit'].attrs['long_name']   = 'Latitudes derived from EWG polar grid'

# Global attributes
dset.attrs['title']       = 'Grid mapping btw EWG Atlas polar sterogr. EASE-Grid and MOM6 mesh025 grid'
dset.attrs['institution'] = 'NOAA NWS NCEP MDC'
dset.attrs['source']      = 'get_gmapi_EWG_snow_Arctic_to_mesh025.py'
dset.attrs['history']     = 'Environmental Working Group Arctic Meteorology and Climate Atlas, Version 1'
dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
dset.attrs['region']      = regn
dset.attrs['Grid_idm_jdm'] = f'{idm}x{jdm}'

fgmapi  = f'EWG_Atlas_MOM6_gmapi_{idm}x{jdm}_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)

print(f'Saving gmapi --> {dfgmapi}')
dset.to_netcdf(dfgmapi, format='NETCDF4', engine='netcdf4')

f_chck = False
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:
  clrmp = mclrmps.colormap_conc()
  rmin = 0.
  rmax = 1.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  LON1 = LON.copy()
  LON2 = LON.copy()
  LON1 = np.where(LON1 < -175, np.nan, LON1)
  LON2 = np.where(LON2 > 172, np.nan, LON2)
  LON3 = np.where(LON < 0, LON+360., LON)
  LON3 = np.where(LON3 > 350., np.nan, LON3)
  lon_cntr1 = [x for x in range(-170,0,10)]  # grey -180:0
  lon_cntr2 = [x for x in range(10,178,10)]  # blue: 0 to 180 E
  if regn == 'south':
    lat_cntr = [x for x in range(-85,-20,5)]
  else:
    lat_cntr = [x for x in range(50,89,5)]

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.12, 0.12, 0.8, 0.8])

  if f_xy: 
    ax1.pcolormesh(XX,YY,C2d, cmap=clrmp, vmin=rmin, vmax=rmax)
    #ax1.invert_yaxis()
    # plot on XX,YY:
    cs = ax1.contour(XX,YY,LON1, lon_cntr1, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LON2, lon_cntr2, linestyles='solid', colors=[(0.,0.5,0.9)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LAT, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LON1,[0], linestyles='solid', colors=[(1,0.,0.)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LAT,[-75], linestyles='solid', colors=[(1,0.,0.)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LON3,[180], linestyles='solid', colors=[(0.6,0.,0.9)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    ax1.axis('scaled')
  else:
    # Plot on grid:
    ax1.pcolormesh(C2d, cmap=clrmp, vmin=rmin, vmax=rmax)
    cs = ax1.contour(LON1, lon_cntr1, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs2 = ax1.contour(LON2, lon_cntr2, linestyles='solid', colors=[(0.,0.5,0.9)], linewidths=1)
    ax1.clabel(cs2, inline=True, fontsize=10, fmt="%.1f")
    cs3 = ax1.contour(LON1,[0], linestyles='solid', colors=[(1,0.,0.)], linewidths=2)
    ax1.clabel(cs3, inline=True, fontsize=12, fmt="%.1f")
    ax1.contour(LON3,[180], linestyles='solid', colors=[(0.6,0.,0.9)], linewidths=1)
    cs = ax1.contour(LAT, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    ax1.axis('scaled')
    ax1.invert_yaxis() 
    ax1.set_ylabel('Inverted j index')
    ax1.set_xlabel('i index')

  ax1.set_title('Derived lon/lat from EWG polar sterogr. projection')
  btx = 'get_gmapi_EWG_snow_Arctic_to_mesh025.py'
  bottom_text(btx) 


