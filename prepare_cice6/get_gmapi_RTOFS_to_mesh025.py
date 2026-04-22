"""
  Find RTOFS indices for bilinear interpolation onto mesh025 global grid
  Polar regions only (bounded by 50 N/S)


"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import time
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

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_read_hycom as mhycom
import mod_regmom as mrmom

init_date = 20250704  #
init_hr = 0
regn = 'global'  
fhr = 0

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
#parser.add_argument("--init", help=f"init date", choices=[20250103, 20250704], required=True, type=int)
parser.add_argument("--init", help=f"init date", choices=[20251231, 20250704], default=init_date, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--stmp", 
        help="Save temporary files at stmp=percent processed pnts, continue from saved (>0 - 100=yes, 0=no)",
        required=True, type=int)
args = parser.parse_args()

regn      = args.regn if args.regn else regn
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
stmp      = args.stmp 
use_tmp   = stmp > 0

assert stmp > 1 and stmp < 100, 'keep stmp > 1% and <100%'

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


# Read RTOFS - CICE4
# Init date:
dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
yrI, mmI, ddI, hrI = mtime.datevec(dnmbI)[:4]


pthice = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/cice4_fcast/{init_date}"
if fhr == 0:
  flice = f"rtofs_glo.t{init_hr:02d}z.n00.cice_inst.nc"
else:
  flice = f"rtofs_glo.t{init_hr:02d}z.f{fhr:02d}.cice_inst.nc"
dflice = os.path.join(pthice, flice)

# Get grid
with xarray.open_dataset(dflice) as dcice:
  LON  = dcice["TLON"].data
  LAT  = dcice["TLAT"].data

JDIM, IDIM = LON.shape
JDIM = JDIM + 1   # ocean grid has + 1 row

# Read RTOFS topo:
# Note that RTOFS grid has +1 row at the top compared to CICE6
pthtopo = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/topo_grid/'
ftopo  = 'depth_GLBb0.08_09m11'
HH0 = mhycom.read_topo(pthtopo, ftopo, IDIM, JDIM)
HH0 = HH0[:-1,:]     # discard the extra row


# Find lat bounds of RTOFS data:
# MOM6 points should be inside the RTOFS domain 
# to be able to find 4-vertice of the bounding box for interpolation
lat_min = np.min(LAT)
lat_max = np.max(LAT)
lat_sh = -50
lat_nh = 50
if regn == 'south':
  lat_max = lat_sh
elif regn == 'north':
  lat_min = lat_nh

if regn == 'north' or regn == 'global':
  lat_max = 90. # override np.max(LAT) for Polar stereogr. projection, the code should work
                # for correctly finding 4 points around hlat>np.max(LAT) but only
                # for Polar stereogr. projection by grabbing points over the N. Pole 

ignore_north_lim = lat_max >= 90.
if ignore_north_lim:
  print(f"WARN: Indices north of northernmost lat={np.max(LAT):.4f} will be searched")

row_min = np.min(hlat, axis=1)  # min lat in each row
row_max = np.max(hlat, axis=1)  # max lat in each row
jS = np.argmax(row_min >= lat_min)
jE = len(row_max) - np.argmax(row_max[::-1] <= lat_max) - 1 # note reverse indexing for [::-1]
#jE = np.argmin(row_max <= lat_max) - 1
jdm, idm = hlon.shape

# Output dir, temp file:
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump = os.path.join(pthdata,'gmapi_mapping2mesh025')
pthtmp  = os.path.join(pthdata,'gmapi_mapping2mesh025','tmp')
fltmp   = f'RTOFS008_mesh025_gmapi_{regn}_tmp.npz'
flgmapi = f'RTOFS008_mesh025_gmapi_{idm}x{jdm}_{regn}.nc'
dfltmp  = os.path.join(pthtmp, fltmp)
dflgmapi = os.path.join(pthdump, flgmapi)

INDX = None
JNDX = None
INDX_list = []
JNDX_list = []
IMOM = []
JMOM = []

# if saving tmp file, check if file exists:
if use_tmp:
  if os.path.isfile(dfltmp):
    print(f"Loading temp. file {dfltmp}")
    data = np.load(dfltmp)
    IMOM = data["imom"].tolist()
    JMOM = data["jmom"].tolist()
    INDX_list.append(data["indx"])
    JNDX_list.append(data["jndx"])
  else:
    print(f"No saved processed pnts, Temp. file does not exist {dfltmp}") 

 
JDIM, IDIM = hlon.shape
IJDIM = IDIM*JDIM
if regn == 'global':
  PMsk = ((hlat <= lat_sh) | (hlat >= lat_nh)) & (HH <= 0)
elif regn == 'north':
  PMsk = (hlat >= lat_nh) & (HH <= 0)
elif regn == 'south':
  PMsk = (hlat <= lat_sh) & (HH <= 0)

Npnts_mom = np.sum(PMsk)
step_dump = int(stmp/100. * Npnts_mom)
if use_tmp > 0:
  print(f"Temporary files dumped every icc = {step_dump} step")
  
pnts_saved = set(zip(IMOM, JMOM))
icc = -1
time0 = time.perf_counter()
for ii in range(idm):
  for jj in range(jS,jE+1):
    if HH[jj,ii] >= 0:
      continue
    x0 = hlon[jj,ii]
    y0 = hlat[jj,ii]

    # Skip grid points outside the region:
    if regn == 'global':
      if lat_nh > y0 > lat_sh:
        continue
    else:
      if y0 < lat_min or y0 > lat_max:
        continue

    icc += 1
    if icc%5000 == 0:
      prc_done = icc / Npnts_mom * 100.
      time1 = time.perf_counter()
      dltT = (time1 - time0)/60.
      print(f'     icc={icc} {prc_done:.2f}% done, time={dltT:.3}min ...')
      time0 = time.perf_counter()

    if use_tmp and (ii,jj) in pnts_saved:
      continue

    #tt0 = time.perf_counter()
    #print(f"icc={icc}, ii={ii}, jj={jj}, x0={x0}, y0={y0}")
    ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=.8, 
                                       ignore_north_lim=ignore_north_lim, 
                                       wrap_long=True)
    #tt1 = time.perf_counter()
    #print(f" dT = {tt1-tt0} sec, estimated tot time={(tt1-tt0)*Npnts_mom/3600} hrs")

    if len(ixx)==0 or len(jxx)==0:
     continue
    ixx = np.expand_dims(ixx, axis=0)
    jxx = np.expand_dims(jxx, axis=0)

    INDX_list.append(ixx)
    JNDX_list.append(jxx)
    IMOM.append(ii)
    JMOM.append(jj)

    if use_tmp and icc > 0 and icc % step_dump == 0:
      # Saving tmp file:
      print(f"  icc={icc} Dumping tmp file --> {dfltmp}")
      np.savez(
               dfltmp, 
               imom=np.array(IMOM), 
               jmom=np.array(JMOM), 
               indx=np.concatenate(INDX_list, axis=0),
               jndx=np.concatenate(JNDX_list, axis=0)
               ) 

IMOM = np.array(IMOM)
JMOM = np.array(JMOM)
INDX = np.concatenate(INDX_list, axis=0)
JNDX = np.concatenate(JNDX_list, axis=0)

# RTOFS coord dimensions:
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
dset['gmapi_i'].attrs['long_name'] = 'I indices RTOFS grid for interpolation'
dset['gmapi_j'].attrs['long_name'] = 'J indices RTOFS grid for interpolation'
dset['longit'].attrs['long_name']  = 'Longitudes derived from RTOFS polar grid'
dset['latit'].attrs['long_name']   = 'Latitudes derived from RTOFS polar grid'

# Global attributes
dset.attrs['title']       = f'Grid mapping between RTOFS CICE4 and MOM6 mesh025 grid'
dset.attrs['institution'] = 'NOAA NWS NCEP EMC'
dset.attrs['source']      = 'get_gmapi_RTOFS_to_mesh025.py'
dset.attrs['history']     = f'RTOFS CICE4 grid  {flice}'
dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
dset.attrs['region']      = regn

fgmapi  = f'RTOFS_NRTice_MOM6_gmapi_{idm}x{jdm}_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)

print(f'Saving gmapi --> {dfgmapi}')
dset.to_netcdf(dfgmapi, format='NETCDF4', engine='netcdf4')

f_chck = False
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:
  # pnt x0/y0: 74.125/-50.697 outside or near the boundary: i/j=0/766
  # WARN: pnt x0/y0: 74.107/68.533 outside or near the S/N boundary: i/j=4017/3296, skipping
  #x0 = 74.125
  #y0 = -50.697
  ii0, jj0 = mhycom.find_indx_lonlat(74.107, 68.533, hlon, hlat)
  x0 = hlon[jj0,ii0]
  y0 = hlat[jj0,ii0]
  ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=.3,
                                       ignore_north_lim=ignore_north_lim,
                                       wrap_long=True)

  ii0 = 293
  jj0 = 1070
  x0 = hlon[jj0,ii0]
  y0 = hlat[jj0,ii0]

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.12, 0.12, 0.8, 0.8])

  # Lon/lat space:
  X = LON[jxx,ixx]
  Y = LAT[jxx,ixx]
  ax1.plot(LON, LAT, '.', color=[0.8,0.8,0.8])   # RTOFS grid
  ax1.plot(x0, y0, 'r*')     # pnt of interest
  ax1.plot(X, Y, '.-')
  ax1.plot([X[0], X[-1]], [Y[0], Y[-1]], '-')

  # Reset:
  #ax1.autoscale()
  #ax1.relim()

  dx=0.15
  dy=0.15
  ax1.set_xlim([x0-dx, x0+dx])
  ax1.set_ylim([y0-dy, y0+dy])

  # Show lon/ lat
  plt.cla()
  ax1.contour(HH0, [0], colors=[(0.8,0.8,0.8)])
  ax1.contour(LON, [x0], colors=[(0.,0.8,1)])
  ax1.contour(LAT, [y0], colors=[(0.,1,0.5)])
  ih0, jh0 = mhycom.find_indx_lonlat(x0, y0, LON, LAT)  # RTOFS indices
  ax1.plot(ih0, jh0, 'r*')     # pnt of interest





