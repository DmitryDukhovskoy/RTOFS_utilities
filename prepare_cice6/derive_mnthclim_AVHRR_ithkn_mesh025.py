"""
  Derive monthly clim iithkn from a range of years for specified months
  using AVHRR L4 data
  Note some dates may be missing in AVHRR varies for summer and winter seasons

  Interpolate to mesh 025
  for south / north regions


  First create tmp files for all months without final saving, e.g.
  Create, interp and save tmp for Jan - Month (3 tmp files)
  run derive_mnthclim_AVHRR_ithkn_mesh025.py --regn north --ms 1 --me 3 --ys 2015 --ye 2025 --fsave 0

   Do this for all 12 months (saving 1 or 2 months at a time)

  Final:
  run derive_mnthclim_AVHRR_ithkn_mesh025.py --regn north --ms 1 --me 12 --ys 2015 --ye 2025 --fsave 1


    NOAA Climate Data Record (CDR) of AVHRR Polar Pathfinder Extended (APP-X) Cryosphere, Version 2

  Albedo exists only when solar senith angle > 0
  For winter months, most of the polar region = NaN
  and 2 am albedo = NaN for most months except summer

  NOAA Climate Data Record (CDR) of the eXtended AVHRR Polar Pathfinder (APP-X) 
  cryosphere contains 19 geophysical variables over the Arctic and Antarctic for the period 1982 - present. 

  https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C01580


The Polar Pathfinder - Extended Climate Data Record (CDR), utilizes the Advanced Very High Resolution Radiometer (AVHRR) and Visible Infrared Imaging Radiometer Suite (VIIRS) instruments and contains several geophysical variables over the Arctic and Antarctic from 1982–present. The data products are mapped to a 25 km Equal-Area Scalable Earth (EASE) grid at two local solar times: 04:00 and 14:00 for the Arctic, and 02:00 and 14:00 for the Antarctic. 


Cite as: Key, Jeffrey; Wang, Xuanji; Liu, Yinghui; and NOAA CDR Program (2019). NOAA Climate Data Record of AVHRR Polar Pathfinder Extended (APP-X), Version 2. [indicate subset used]. NOAA National Centers for Environmental Information. doi:10.25921/AE96-0E57 [access date].


  gmapi indices:
  get_gmapi_AVHRR_albedo_to_mesh025.py

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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)  

fld_name = 'ithkn'  # albedo or ithkn
tmpf = 1  # for climatology, save temporary monthly and start from last saved

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"Region to process",
                    choices=['north','south'], required=True, type=str)
parser.add_argument("--ms",  help="Start clim month", required=True, type=int)
parser.add_argument("--me",  help="End clim month default=ms", type=int)
parser.add_argument("--ys", help="Year to start clim", choices=[x for x in range(1982,2026)],
                    required=True, type=int)
parser.add_argument("--ye", help="Year to end clim, >=ys, default=ys", type=int)
parser.add_argument("--fsave", help=f"Final Save collecting all tmp monthly files 0=no, 1=yes",
                    required=True, choices=[0,1], type=int)
args = parser.parse_args()

regn   = args.regn if args.regn is not None else None
MMS    = args.ms
MME    = args.me if args.me is not None else MMS
YRS    = args.ys
YRE    = args.ye if args.ye is not None else YRS
save_tmp = tmpf == 1
save_final = bool(args.fsave)


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

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
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
jdim, idim = HH.shape
LMsk = np.where(HH<0, 1, 0)

# Mask out not needed latitudes:
if regn == 'south':
  LMsk = np.where(hlat > -55, 0, LMsk)
else:
  LMsk = np.where(hlat < 50, 0, LMsk)

# Get gmapi 4 NSIDC grid points for interpolation
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthdump  = os.path.join(pthdata,"gmapi_mapping2mesh025")
fgmapi  = f'AVHRR_albedo_MOM6_gmapi_{idim}x{jdim}_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)
print(f'Loading gmapi --> {dfgmapi}')

dgmapi = xarray.open_dataset(dfgmapi)
IMOM = dgmapi['mom_indx'].data
JMOM = dgmapi['mom_jndx'].data
INDX = dgmapi['gmapi_i'].data
JNDX = dgmapi['gmapi_j'].data
LON  = dgmapi['longit'].values  # original gird, longit
LAT  = dgmapi['latit'].values   # original grid, latit

# Convert to -180, 180:
LON = (LON + 180.) % 360. - 180.

# Parse url file names to extract time
import re
def parse_filename(fname):
  # Extract time fomr the file name
  # expexted format: *_1400_d20251231_c20260109.nc

  pattern = r'_(\d{4})_d(\d{8})_c(\d{8})'
  m = re.search(pattern, fname)
  
  if not m:
    raise ValueError(f"Time cannot be extracted, url filename {fname}  does not match format")
  
  time_str = m.group(1)   # '1400'
  date_str = m.group(2)   # '20251231'
  creation_str = m.group(3)

  hour = int(time_str[:2])
  minute = int(time_str[2:])

  return hour, minute, date_str, creation_str

# AVHRR data from OPeNDAP server
# Files naming is irregular in terms of N days composite (?)
from urllib.request import urlopen
import xml.etree.ElementTree as ET

def list_files(year, month, day, regn="north"):
  """
    Download directory listing (XML)
    Convert XML into structured Python object using xml.etree.ElementTree
    Search for filenames matching date requested, e.g. d20251231_cYYYYMMDD <-- unknown
  """
  base = "nhem" if regn == "north" else "shem"
  cat_url = (
      f"https://www.ncei.noaa.gov/thredds/catalog/"
      f"avhrr-polar-pathfinder-ext-files/{base}/{year}/catalog.xml"
  )
 
  print(f"Reading xml catalog from {cat_url}")
  with urlopen(cat_url) as response:
    xml_data = response.read()

  root = ET.fromstring(xml_data)

  ns = {"t": "http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0"}
  target = f"d{year}{month:02d}{day:02d}"

  files = []
  for ds in root.findall(".//t:dataset", ns):
    name = ds.attrib.get("name", "")
    if target in name:
      files.append(name)

  return files

if save_final:
  pthice   = os.path.join(pthdata, 'AVHRR_albedo_ithkn','tmp')
  fliceout = f'AVHRR_{fld_name}_mnthclim_{YRS}-{YRE}_{idim}x{jdim}_{regn}.nc'
  dfliceout = os.path.join(pthice,fliceout)
  if os.path.isfile(dfliceout):
    print(f"Allready exists {dfliceout} ...")
    raise RuntimeError(f"rename/ delete existing file {dfliceout}")

# Save temporary clim fields:
if save_tmp:
  tmp_dir = os.path.join(pthdata, 'AVHRR_albedo_ithkn','tmp')
  os.makedirs(tmp_dir, exist_ok=True)

A3d = np.zeros((12,jdim,idim))

for MM in range(MMS, MME+1):
  imo = MM - 1

  # Check if tmp file exists:
  if save_tmp:
    tmp_file = os.path.join(tmp_dir, f"tmp_AVHRR_{fld_name}_{YRS}-{YRE}_{MM:02d}.npy")
    if os.path.exists(tmp_file):
      print(f"Skipping month {MM}: already computed {tmp_file}")
      A2d = np.load(tmp_file)
      A3d[imo,:,:] = A2d
      continue

  irecs = 0
  ASUM = np.full_like(LON,0)
  count_ice = np.full_like(LON,0)
  for YR in range(YRS,YRE+1):
    mdays = mtime.month_days(MM,YR)
    for DD in range(1,mdays+1):
      date_req= YR*10000 + MM*100 + DD
      dnmb_req = mtime.dateint2datenum(date_req)  # requested date numb
      YR, MM, DD = mtime.datevec(dnmb_req)[:3]

      # Supposedly, should be 2 overpasses: 4:00 and 14:00 for the Arctic Ocean local time
      # and 2:00 and 14:00 for S. Ocean - average
      url_files = list_files(YR, MM, DD, regn)

      # Find 2 matching files for different times
      nfiles = len(url_files)

      if nfiles <2:
        # Assumed 2 passes: 4:00 and 14:00 
        print(f"No url files found for {date_req}")
        continue

      fl1 = url_files[0]
      hr1, _, dstmp1, cstmp1 = parse_filename(fl1)
      fl2 = None
      for fll in url_files[1:]:
        hr2, _, dstmp2, cstmp2 = parse_filename(fll)
        if hr2 != hr1 and dstmp2 == dstmp1 and cstmp2 == cstmp1:
          fl2 = fll
          break

      print(f"Found 2 times matching {date_req}:")
      print(f"{fl1} \n{fl2}")

      if regn == 'north':
        url = f"https://www.ncei.noaa.gov/thredds/dodsC/avhrr-polar-pathfinder-ext-files/nhem/{YR}/"
      else:
        url = f"https://www.ncei.noaa.gov/thredds/dodsC/avhrr-polar-pathfinder-ext-files/shem/{YR}/"

      dflinp1 = os.path.join(url, fl1)
      dflinp2 = os.path.join(url, fl2)

      if fld_name == 'albedo':
        varnm = 'cdr_surface_albedo'
      elif fld_name == 'ithkn':
        varnm = 'cdr_sea_ice_thickness'

      print(f"Reading {dflinp1}")
      with xarray.open_dataset(dflinp1) as ds_ices:
        AA1 = ds_ices[varnm].values.squeeze()

      print(f"Reading {dflinp2}")
      with xarray.open_dataset(dflinp2) as ds_ices:
        AA2 = ds_ices[varnm].values.squeeze()

      if AA1.shape != LON.shape:
        raise RuntimeError(f"Check data / lon/lat shapes")

      # Average ignoring NaNs
      Mnan1 = np.isnan(AA1)
      Mnan2 = np.isnan(AA2)
      Mvalid = ~Mnan1 | ~Mnan2  # at least 1 good value

      AA = np.full_like(AA1, np.nan)

      # Average nonnans in AA2 and AA1 or simply copy if one is NaN:
      AA[Mvalid] = (np.nan_to_num(AA1[Mvalid]) + np.nan_to_num(AA2[Mvalid])) / (
                    (~Mnan1[Mvalid]).astype(int) + (~Mnan2[Mvalid]).astype(int) 
                    )

      # Count & average only non-zero thicknesses
      # Do not count 0 thickn. this will bias ice thickness 
      ice_grid = np.isfinite(AA) & (AA > 0.0)
      ASUM[ice_grid] += AA[ice_grid]
      count_ice[ice_grid] += 1
      irecs += 1

  # Average all records for this month:
  print(f"Averaging, nrecs={irecs}, max count_ice={np.max(count_ice)}")
  Aavrg = np.divide(ASUM, count_ice, out=np.zeros_like(ASUM), where = count_ice > 0)
  Aavrg[np.isnan(AA)] = np.nan
  
  print(f"Interpolating {fld_name} ...")
  CIint = msisrlx.interp2Dfld(Aavrg, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat, land_mask=True)
  CIint = np.where(HH>=0, np.nan, CIint)
  A2d = np.expand_dims(CIint, axis=0)
  A3d[imo,:,:] = A2d

  if save_tmp:
    print(f"Saving temporary --> {tmp_file}")
    np.save(tmp_file, A2d)


if save_final:
  A3d = A3d.astype('float32')
  JD = np.arange(jdim, dtype='int32')
  ID = np.arange(idim, dtype='int32')

  time_obs = np.array([x for x in range(1,13)])

  darr_hi = xarray.DataArray(A3d, dims=("time","jdim","idim"),
                     coords={"time": time_obs,
                             "jdim": JD,
                             "idim": ID,})
  darr_lon = xarray.DataArray(hlon, dims=("jdim","idim"),
                   coords={"jdim": JD,\
                           "idim": ID,})
  darr_lat = xarray.DataArray(hlat, dims=("jdim","idim"),
                   coords={"jdim": JD,\
                           "idim": ID,})

  if fld_name == 'ithkn':
    dset_hi = xarray.Dataset({
      "ice_thkn": darr_hi,
      "lon": darr_lon,
      "lat": darr_lat,
    })
    dset_hi['ice_thkn'].attrs.update({
      "long_name": "sea ice thickness from CDR AVHRR",
      "units": "m",
    })
  elif fld_name == 'albedo':
    dset_hi = xarray.Dataset({
      "albedo": darr_hi,
      "lon": darr_lon,
      "lat": darr_lat,
    })

    dset_hi['albedo'].attrs.update({
      "long_name": "ice / snow surface broadband albedo",
      "units": "1",
    })

  dset_hi["time"].attrs.update({
      "long_name": "time",
      "units": "Months",
  })

  dset_hi['lon'].attrs.update({
    "long_name": "Longitudes",
    "units": "degrees_east",
  })
  dset_hi['lat'].attrs.update({
    "long_name": "Latitudes",
    "units": "degrees_north",
  })

  dset_hi.attrs.update({
    "title": f"CDR AVHRR APP-X {fld_name} climatology for {YRS}-{YRE} interpolated to mesh025 grid",
    "info": "https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C01580",
    "institution": "NOAA NWS MDC",
    "source": "derive_mnthclim_AVHRR_ithkn_mesh025",
    "contact": "dmitry.dukhovskoy@noaa.gov",
    "region": regn,
  })

  print(f'Dumping interpolated {fld_name} --> {dfliceout}')
  dset_hi.to_netcdf(dfliceout,
          encoding={var: {'_FillValue': np.float32(1e30)} for var in dset_hi.data_vars},
          format='NETCDF4')


f_chck = False
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 1.
  if fld_name == 'ithkn':
    rmax = 3.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
    parallels = np.arange(-80,-10,10.)
    meridians = np.arange(-360,359.,45.)
  elif regn == 'north':
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
    parallels = np.arange(50, 90, 5)
    meridians = np.arange(-360, 359., 45.)

  xh, yh = m(hlon,hlat) # GFS coords

  lon_cntrs = [x for x in range(-180,180,45)]
  lat_cntrs = [x for x in range(-80,90,10)]

  A2d = np.squeeze(A2d)

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.08, 0.1, 0.82, 0.82])

  # Plot original field on grid:
  img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  m.drawparallels(parallels, labels=[0,0,0,0])
  m.drawmeridians(meridians, labels=[0,0,0,0])

  sttl = f'ithkn AVHRR clim MM={MM:02d}'
  ax1.set_title(sttl, fontsize=10)


  ax3 = fig1.add_axes([0.1, 0.06, 0.8, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'derive_mnthclim_AVHRR_ithkn_mesh025.py'
  bottom_text(btx, pos=[0.08,0.02], fsz=8) 



