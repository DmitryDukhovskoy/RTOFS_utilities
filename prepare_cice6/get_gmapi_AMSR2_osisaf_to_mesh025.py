"""
  Plotting indces from AMSR2 grid --> mesh0.25 grid
  4 vertices for bilinear interpolation

  AMSR2 L4 OSI SAF EUMETSAT sea ice concentration 
  on native grid
  daily fields

    string :title = "Level-3 Sea Ice Concentration Analysis (AMSR2) from OSI SAF EUMETSAT" ;
    string :institution = "EUMETSAT OSI SAF" ;
    string :history = "Created 2025-05-13 02:01:07" ;
    string :contact = "osisaf-manager@met.no" ;
    string :references = "Product User Manual and Algorithm Theoretical Basis Document available at https://osi-saf.eumetsat.int/documentation/products-documentation" ;
    string :subset :source = "ARCO data downloaded from the Marine Data Store using the MyOcean Data Portal" ;

Downloaded from:
https://data.marine.copernicus.eu/product/SEAICE_ARC_PHY_AUTO_L4_MYNRT_011_024/services

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray as xr
import time
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

init_hr = 0
regn = 'north'  
fhr = 0
chsize = 360   # for saving multiple tmp files by chunks 0:143, 144:287, etc 

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south", required=True, type=str)
parser.add_argument("--rdate", help="AMSR2 data date to read grid, YYYYMMDD", required=True, type=int)
parser.add_argument("--stmp", 
    help="Save temp files at stmp-perc processed pnts & continue from saved (1 - 100perc, 100-only final save)",
    default=10,
    type=int)
parser.add_argument("--chsize", 
    help=f"Chunk size (diviser of i-dim): numb of i indices 1,2,4,8,10, .., idm, default={chsize}", 
    default=chsize,
    type=int)
parser.add_argument("--kchunk", 
        help="Not required for final, for tmp chunk number: 1,2,... (idm/chsize)", 
        type=int)
parser.add_argument("--final", help="Finish by combining all tmp files 1=yes, 0=no", 
                    default=0, choices=[0,1], type=int)
args = parser.parse_args()

regn    = args.regn if args.regn is not None else regn
rdate   = args.rdate
stmp    = args.stmp
use_tmp = stmp > 0
chsize  = args.chsize 
kchunk  = args.kchunk if args.kchunk is not None else None
final   = args.final == 1

assert stmp == 0 or (1 < stmp < 100), 'keep stmp=0 or 1<stmp<100'


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
jdm, idm = hlon.shape

with xr.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

# Check chunks:
assert idm % chsize == 0, f"chsize={chsize} is not an integer divisor of {idm}"

if not final:
  assert 1 <= kchunk <= idm // chsize, \
    f"valid kchunk range for chsize={chsize} is [1, {idm // chsize}], " \
    f"given kchunk={kchunk}"
  iS = (kchunk-1) * chsize
  iE = kchunk * chsize
else:
  iS = 0
  iE = idm


dnmb0 = mtime.rdate2datenum(rdate)
YR, MM, DD = mtime.datevec(dnmb0)[:3]

# Read AMSR2 grid:
pthamsr = f'/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/AMSR2_OSISAF_L4/{regn}/{YR}'
flamsr = f"osisaf_{regn}_nrt_amsr2_L4_{MM:02d}{YR}.nc"
dfl = os.path.join(pthamsr,flamsr)
assert os.path.isfile(dfl), f"Missing data: {dfl}"

with xr.open_dataset(dfl) as ds_amsr:
  AA = ds_amsr['ice_conc'].isel(time=DD-1).values.squeeze() * 0.01  # % to fractions
  lon1d = ds_amsr['longitude'].values
  lat1d = ds_amsr['latitude'].values

LON, LAT = np.meshgrid(lon1d, lat1d)
JDIM, IDIM = LON.shape

# Find lat bounds of AMSR2 data:
# MOM6 points should be inside the AMSR2 domain 
# to be able to find 4-vertice of the bounding box for interpolation
lat_min = np.min(LAT)
lat_max = np.max(LAT)
if regn == 'south':
  lat_sh = min(lat_max, -50)
  lat_max = lat_sh  
elif regn == 'north':
  lat_nh = max(50, lat_min)
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
jE = jdm

# Output dir, temp file:
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump = os.path.join(pthdata,'gmapi_mapping2mesh025')
pthtmp  = os.path.join(pthdata,'gmapi_mapping2mesh025','tmp')
flgmapi = f'AMSR2_grid0p1_to_mesh025_gmapi_{idm}x{jdm}_{regn}.nc'
dflgmapi = os.path.join(pthdump, flgmapi)

if final:
  print(f"   Final file will be combined from tmp files and saved")
else:
  fltmp   = f'AMSR2_to_mesh025_gmapi_iS{iS:04d}_iE{iE:04d}_{kchunk:03d}.npz'
  dfltmp  = os.path.join(pthtmp, fltmp)
  print(f"   Final file will not be saved,  tmp file={fltmp}")

def read_tmp_file(dfltmp): 
  """
    Update gmapi indices from tmp files
  """
  INDX_list = []
  JNDX_list = []
  IMOM = []
  JMOM = []
  if os.path.isfile(dfltmp):
    print(f"Loading temp. file {dfltmp}")
    data = np.load(dfltmp)
    IMOM = data["imom"].tolist()
    JMOM = data["jmom"].tolist()
    INDX_list.append(data["indx"])
    JNDX_list.append(data["jndx"])
  else:
    print(f"No saved processed pnts, Temporary file does not exist {dfltmp}") 

  return IMOM, JMOM, INDX_list, JNDX_list

def save_tmp(dfltmp, IMOM, JMOM, INDX_list, JNDX_list):
    np.savez(
        dfltmp,
        imom=np.array(IMOM),
        jmom=np.array(JMOM),
        indx=np.concatenate(INDX_list, axis=0),
        jndx=np.concatenate(JNDX_list, axis=0)
    )



def combine_tmp_files(pthtmp):
  tmp_files = sorted([
      f for f in os.listdir(pthtmp)
      if f.startswith("AMSR2_to_mesh025_gmapi_") and f.endswith(".npz")
  ])

  if len(tmp_files) == 0:
    raise RuntimeError(f"No tmp files found in {pthtmp}")

  IMOM_all = []
  JMOM_all = []
  INDX_all = []
  JNDX_all = []

  for fl in tmp_files:
    print(f"Reading {fl}")
    dfl = os.path.join(pthtmp, fl)
    d = np.load(dfl, allow_pickle=True)

    imom = np.asarray(d["imom"])
    jmom = np.asarray(d["jmom"])
    indx = np.asarray(d["indx"])
    jndx = np.asarray(d["jndx"])

    # skip empty files
    if imom.size == 0 or indx.size == 0:
      continue

    # fix shape issues
    if indx.ndim == 1:
      indx = indx.reshape(1, -1)
    if jndx.ndim == 1:
      jndx = jndx.reshape(1, -1)

    IMOM_all.append(imom)
    JMOM_all.append(jmom)
    INDX_all.append(indx)
    JNDX_all.append(jndx)

  # final merge
  IMOM = np.concatenate(IMOM_all)
  JMOM = np.concatenate(JMOM_all)
  #INDX = np.vstack(INDX_all)
  #JNDX = np.vstack(JNDX_all)

  return IMOM, JMOM, INDX_all, JNDX_all


def save_netcdf(IMOM, JMOM, INDX_list, JNDX_list, LON, LAT, dfgmapi):
  """
    Final NetCDF writer
  """
  IMOM = np.array(IMOM)
  JMOM = np.array(JMOM)
  npnts = len(IMOM)
  if isinstance(INDX_list, np.ndarray):
      INDX = INDX_list
      JNDX = JNDX_list

  elif isinstance(INDX_list, list):
      INDX = np.vstack([np.asarray(x) for x in INDX_list])
      JNDX = np.vstack([np.asarray(x) for x in JNDX_list])

  else:
      raise TypeError(f"Unrecognized type: {type(INDX_list)}")

  INDX = np.asarray(INDX)
  JNDX = np.asarray(JNDX)

  assert INDX.shape == (npnts, 4), f"INDX shape wrong: {INDX.shape}"
  assert JNDX.shape == (npnts, 4), f"JNDX shape wrong: {JNDX.shape}"

  # AMSR2 coord dimensions:
  jdim, idim = LON.shape

  darr_imom = xr.DataArray(IMOM, dims=("npoints"),\
                     coords={"npoints": np.arange(npnts)})
  darr_jmom = xr.DataArray(JMOM, dims=("npoints"),\
                     coords={"npoints": np.arange(npnts)})
  darr_indx = xr.DataArray(INDX, dims=("npoints","nvert"),\
                     coords={"npoints": np.arange(npnts),\
                             "nvert": np.arange(4)})
  darr_jndx = xr.DataArray(JNDX, dims=("npoints","nvert"),\
                     coords={"npoints": np.arange(npnts),\
                             "nvert": np.arange(4)})
  darr_lon = xr.DataArray(LON, dims=("jdim","idim"),\
                    coords={"jdim": np.arange(jdim),\
                            "idim": np.arange(idim)})
  darr_lat = xr.DataArray(LAT, dims=("jdim","idim"),\
                    coords={"jdim": np.arange(jdim),\
                            "idim": np.arange(idim)})

  dset = xr.Dataset({"mom_indx": darr_imom, \
                         "mom_jndx": darr_jmom, \
                         "gmapi_i": darr_indx,\
                         "gmapi_j": darr_jndx,\
                         "longit":  darr_lon,\
                         "latit":   darr_lat})

  dset['mom_indx'].attrs['long_name'] = 'MOM6 grid I indices corresponding gmapi'
  dset['mom_jndx'].attrs['long_name'] = 'MOM6 grid J indices corresponding gmapi'
  dset['gmapi_i'].attrs['long_name'] = 'I indices AMSR2 grid for interpolation'
  dset['gmapi_j'].attrs['long_name'] = 'J indices AMSR2 grid for interpolation'
  dset['longit'].attrs['long_name']  = 'Longitudes derived from AMSR2 polar grid'
  dset['latit'].attrs['long_name']   = 'Latitudes derived from AMSR2 polar grid'

  # Global attributes
  dset.attrs['title']       = f'Grid mapping between AMSR2 {regn} polar region onto UFS mesh025 grid'
  dset.attrs['institution'] = 'NOAA NWS OMD'
  dset.attrs['source']      = 'get_gmapi_AMSR2_osisaf_to_mesh025.py'
  dset.attrs['history']     = f'AMSR2 grid:  {dfl}'
  dset.attrs['region']      = regn

  print(f'Saving gmapi --> {dfgmapi}')
  dset.to_netcdf(dfgmapi, format='NETCDF4', engine='netcdf4')


INDX = None
JNDX = None
INDX_list = []
JNDX_list = []
IMOM = []
JMOM = []

if not final:
  IMOM, JMOM, INDX_list, JNDX_list = read_tmp_file(dfltmp)
else:
  # For final stage: combine all temporary files that exist:
  IMOM, JMOM, INDX_list, JNDX_list = combine_tmp_files(pthtmp)

JDIM, IDIM = hlon.shape
IJDIM = IDIM*JDIM

# Points to process:
PMsk = (HH < 0.)
PMsk &= (hlat >= lat_min)
if regn == 'global':
  PMsk &= (hlat <= lat_sh) | (hlat >= lat_nh)
elif regn == 'north':
  PMsk &= (hlat >= lat_nh)
elif regn == 'south':
  PMsk &= (hlat <= lat_sh)

JMSK, IMSK = np.where(PMsk)

Npnts_mom = np.sum(PMsk)
step_dump = max(100, int(stmp/100. * Npnts_mom))
if not final:
  print(f"Temporary files dumped every icc = {step_dump} steps")
  
pnts_saved = set(zip(IMOM, JMOM))

# Find not matched PMsk
#JM, IM = np.where(PMsk)
#mask_pnts = set(zip(IM,JM))
#for im, jm in mask_pnts:
#  if (im,jm) not in pnts_saved:
#    print(f"Missing jm={jm} im={im}")


print(f" Processed in tmp files: {len(pnts_saved)}, total pnts={Npnts_mom}\n")

icc = -1
time0 = time.perf_counter()
for ii in range(iS, iE):
  for jj in range(jS,jE):
    x0 = hlon[jj,ii]
    y0 = hlat[jj,ii]

    # Skip grid points outside the region:
    if not PMsk[jj, ii]:
      continue

    #print(f"icc={icc} ii={ii} jj={jj}")

    icc += 1
    if icc%2000 == 0:
      prc_done = icc / Npnts_mom * 100.
      time1 = time.perf_counter()
      dltT = (time1 - time0)/60.
      print(f'     icc={icc} {prc_done:.2f}% done, time={dltT:.3}min ...')
      time0 = time.perf_counter()

    if (ii,jj) in pnts_saved:
      continue

    #tt0 = time.perf_counter()
    #print(f"icc={icc}, ii={ii}, jj={jj}, x0={x0}, y0={y0}")
    ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=.2, 
                                       ignore_north_lim=ignore_north_lim, 
                                       wrap_long=True)
    #tt1 = time.perf_counter()
    #print(f" dT = {tt1-tt0} sec, estimated tot time={(tt1-tt0)*Npnts_mom/3600} hrs")

    if len(ixx)==0 or len(jxx)==0:
     continue
    ixx = np.expand_dims(ixx, axis=0)
    jxx = np.expand_dims(jxx, axis=0)
    assert ixx.shape == (1, 4), "unexpected ixx shape"
    assert jxx.shape == (1, 4), "unexpected jxx shape"

    INDX_list.append(ixx)
    JNDX_list.append(jxx)
    IMOM.append(ii)
    JMOM.append(jj)

    if not final and icc > 0 and icc % step_dump == 0:
      # Saving tmp file:
      print(f"  icc={icc} Dumping tmp file --> {dfltmp}")
      save_tmp(dfltmp, IMOM, JMOM, INDX_list, JNDX_list)

if not final:
  # Saving tmp file:
  print(f"Final tmp:  icc={icc} Dumping tmp file --> {dfltmp}")
  save_tmp(dfltmp, IMOM, JMOM, INDX_list, JNDX_list)

  print("Finished without saving final netcdf")

else:
  fgmapi  = f'AMSR2_to_mesh025_gmapi_{idm}x{jdm}_{regn}.nc'
  dfgmapi = os.path.join(pthdump, fgmapi)
  print(f"Final saving NetCDF --> {dfgmapi}")
  save_netcdf(IMOM, JMOM, INDX_list, JNDX_list, LON, LAT, dfgmapi) 



f_chck = False
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:
  # icc=197 ii=360 jj=1078
  ii0 = 360
  jj0 = 1078
  x0 = hlon[jj0,ii0]
  y0 = hlat[jj0,ii0]

  #ii0, jj0 = mhycom.find_indx_lonlat(219.125, -50.699, hlon, hlat)
  #x0 = hlon[jj0,ii0]
  #y0 = hlat[jj0,ii0]
  ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=.3,
                                       ignore_north_lim=ignore_north_lim,
                                       wrap_long=True)

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.12, 0.12, 0.8, 0.8])

  # Lon/lat space:
  X = LON[jxx,ixx]
  Y = LAT[jxx,ixx]
  #ax1.plot(LON, LAT, '.', color=[0.8,0.8,0.8])   # AMSR2 grid
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
  ih0, jh0 = mhycom.find_indx_lonlat(x0, y0, LON, LAT)  # AMSR2 indices
  ax1.plot(ih0, jh0, 'r*')     # pnt of interest





