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
chsize = 144   # for saving multiple tmp files by chunks 0:143, 144:287, etc 

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
#parser.add_argument("--init", help=f"init date", choices=[20250103, 20250704], required=True, type=int)
parser.add_argument("--init", help=f"init date", choices=[20251231, 20250704], default=init_date, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--stmp", 
    help="Save temp files at stmp-perc processed pnts & continue from saved (1 - 100, 100-only final save)",
    required=True, type=int)
parser.add_argument("--chsize", help="Chunk size: numb of i indices 1,..., idm", type=int)
parser.add_argument("--kchunk", help=f"If save tmp, chunk number: 1:idm/chunk size, default={chsize}", 
                    type=int)
parser.add_argument("--final", help="Finish by combining all tmp files 1=yes, 0=no", 
                    default=0, choices=[0,1], type=int)
args = parser.parse_args()

regn      = args.regn if args.regn is not None else regn
init_date = args.init if args.init is not None else init_date
init_hr   = args.ihr if args.ihr is not None else init_hr
stmp      = args.stmp
use_tmp   = stmp > 0

chsize = args.chsize if args.chsize is not None else chsize
kchunk = args.kchunk if args.kchunk is not None else None
final  = args.final == 1

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

if not final:
  assert kchunk > 0, f"kchunk should be > 0, kchunk={kchunk}"
  iS = (kchunk-1) * chsize
  iE = kchunk * chsize

  assert iE <= IDIM, print(f"  WARN:   check kchunk={kchunk} chsize={chsize} ==> iE={iE} > {IDIM}")
else:
  iS = 0
  iE = idm


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
jE = jdm

# Output dir, temp file:
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump = os.path.join(pthdata,'gmapi_mapping2mesh025')
pthtmp  = os.path.join(pthdata,'gmapi_mapping2mesh025','tmp')
flgmapi = f'RTOFS008_mesh025_gmapi_{idm}x{jdm}_{regn}.nc'
dflgmapi = os.path.join(pthdump, flgmapi)

if final:
  print(f"   Final file will be combined from tmp files and saved")
else:
  fltmp   = f'RTOFS008_mesh025_gmapi_iS{iS:04d}_iE{iE:04d}_{kchunk:03d}.npz'
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
      if f.startswith("RTOFS008_mesh025_gmapi_") and f.endswith(".npz")
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

  # RTOFS coord dimensions:
  jdim, idim = LON.shape

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
  dset.attrs['title']       = f'Grid mapping between RTOFS CICE4 polar regions and UFS mesh025 grid'
  dset.attrs['institution'] = 'NOAA NWS OMD'
  dset.attrs['source']      = 'get_gmapi_RTOFS_to_mesh025.py'
  dset.attrs['history']     = f'RTOFS CICE4 grid  {flice}'
  dset.attrs['contact']     = 'dmitry.dukhovskoy@noaa.gov'
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
# Note: PMsk is ~12 points > jS:jE, iS:iE
# because PMsk leaving points at jS-1 
PMsk = (HH < 0.)
PMsk &= (hlat >= lat_min)
if regn == 'global':
  PMsk &= (hlat <= lat_sh) | (hlat >= lat_nh)
elif regn == 'north':
  PMsk &= (hlat >= lat_nh)
elif regn == 'south':
  PMsk &= (hlat <= lat_sh)

Npnts_mom = np.sum(PMsk)
step_dump = max(100, int(stmp/100. * Npnts_mom))
if not final:
  print(f"Temporary files dumped every icc = {step_dump} step")
  
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

    icc += 1
    if icc%5000 == 0:
      prc_done = icc / Npnts_mom * 100.
      time1 = time.perf_counter()
      dltT = (time1 - time0)/60.
      print(f'     icc={icc} {prc_done:.2f}% done, time={dltT:.3}min ...')
      time0 = time.perf_counter()

    if (ii,jj) in pnts_saved:
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
  fgmapi  = f'RTOFS_CICE4_gmapi_{idm}x{jdm}_{regn}.nc'
  dfgmapi = os.path.join(pthdump, fgmapi)
  print(f"Final saving NetCDF --> {dfgmapi}")
  save_netcdf(IMOM, JMOM, INDX_list, JNDX_list, LON, LAT, dfgmapi) 



f_chck = False
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:
  # pnt x0/y0: 74.125/-50.697 outside or near the boundary: i/j=0/766
  # WARN: pnt x0/y0: 74.107/68.533 outside or near the S/N boundary: i/j=4017/3296, skipping
  #x0 = 74.125
  #y0 = -50.697
  #ii0, jj0 = mhycom.find_indx_lonlat(74.107, 68.533, hlon, hlat)
  ii0 = 293
  jj0 = 1070
  x0 = hlon[jj0,ii0]
  y0 = hlat[jj0,ii0]

  ii0, jj0 = mhycom.find_indx_lonlat(219.125, -50.699, hlon, hlat)
  x0 = hlon[jj0,ii0]
  y0 = hlat[jj0,ii0]
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
  #ax1.plot(LON, LAT, '.', color=[0.8,0.8,0.8])   # RTOFS grid
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





