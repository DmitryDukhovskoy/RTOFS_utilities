"""
  Interpolate N-day average albedo (varies for summer and winter seasons)
  for south / north regions

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

fld_name = 'albedo'  # albedo or ithkn

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"Region to process",
                    choices=['north','south'], required=True, type=str)
parser.add_argument("--date", help=f"YYYYMMDD of desired field", required=True, type=int) 
parser.add_argument("--fname", help=f"field to interp, default={fld_name}", type=str,
                    choices=['albedo','ithkn'])
args = parser.parse_args()

regn = args.regn if args.regn is not None else None
date_req = args.date if args.date else None
fld_name = args.fname if args.fname else fld_name

dnmb_req = mtime.dateint2datenum(date_req)  # requested date numb
YR, MM, DD = mtime.datevec(dnmb_req)[:3]

    
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


# Supposedly, should be 2 overpasses: 4:00 and 14:00 for the Arctic Ocean local time
# and 2:00 and 14:00 for S. Ocean - average
url_files = list_files(YR, MM, DD, regn)

# Find 2 matching files for different times
nfiles = len(url_files)
assert nfiles >=2, f"{nfiles} url files found for {date_req}"

fl1 = url_files[0]
hr1, _, dstmp1, cstmp1 = parse_filename(fl1)
fl2 = None
for fll in url_files[1:]:
  hr2, _, dstmp2, cstmp2 = parse_filename(fll)
  if hr2 != hr1 and dstmp2 == dstmp1 and cstmp2 == cstmp1:
    fl2 = fll
    break

if fl2 is None:
  raise Exception(f"Could not find 2nd time for {fl1}")

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

print(f"Interpolating {fld_name} ...")

CIint = msisrlx.interp2Dfld(AA, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat, land_mask=True)
CIint = np.where(HH>=0, np.nan, CIint)

A3d = np.expand_dims(CIint, axis=0)
A3d = A3d.astype('float32')
JD = np.arange(jdim, dtype='int32')
ID = np.arange(idim, dtype='int32')

dnmb_ref = mtime.datenum([1900, 1, 1])
DV = mtime.datevec(dnmb_ref)[:3]
time_obs = dnmb_req - dnmb_ref

darr_hi = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
                   coords={"time": [time_obs],\
                           "jdim": JD,\
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
    "units": f"days since {DV[0]}-{DV[1]:02d}-{DV[2]:02d}",
    "calendar": "standard",
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
  "title": f"CDR AVHRR APP-X {fld_name} interpolated to mesh025 grid",
  "info": "https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C01580",
  "institution": "NOAA NWS NCEP MDC",
  "AVHRR_file1": fl1,
  "AVHRR_file2": fl2,
  "source": "interp_AVHRR_albedo_ithkn_mesh025.py",
  "contact": "dmitry.dukhovskoy@noaa.gov",
  "region": regn,
})

pthice   = os.path.join(pthdata, f'AVHRR_albedo_ithkn')
fliceout = f'AVHRR_{fld_name}_{date_req}_mesh025_{idim}x{jdim}_{regn}.nc'
dfliceout = os.path.join(pthice,fliceout)
print(f'Dumping interpolated {fld_name} --> {dfliceout}')
dset_hi.to_netcdf(dfliceout,
        encoding={var: {'_FillValue': np.float32(1e30)} for var in dset_hi.data_vars},
        format='NETCDF4')


# Plot AA1 and AA2 --> average 
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

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.08, 0.4, 0.4, 0.4])

  # Plot original field on grid:
  img = ax1.pcolormesh(AA, cmap=clrmp, vmin=rmin, vmax=rmax)
  sttl = f'{fld_name} AVHRR\n {fl1}+\n{fl2}'
  ax1.set_title(sttl, fontsize=10)    

  ax1.contour(LAT, lat_cntrs, linestyles='solid', linewidths=1, colors=[(0.4,0.4,0.4)])
  ax1.contour(LON, lon_cntrs, linestyles='solid', linewidths=1, colors=[(0.4,0.4,0.4)])
  if regn == 'north':
    ax1.contour(LAT,[60], linestyles='solid', linewidths=1, colors=[(1,0.4,0)])
  else:
    ax1.contour(LAT,[-55], linestyles='solid', linewidths=1, colors=[(1,0.4,0)])
  ax1.invert_yaxis()
  ax1.axis('scaled') 

  # Interpolated field
  ax2 = plt.axes([0.55, 0.4, 0.4, 0.4])  
  m.drawparallels(parallels, labels=[0,0,0,0])
  m.drawmeridians(meridians, labels=[0,0,0,0])
  m.drawcoastlines()
  ax2.pcolormesh(xh, yh, CIint, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax2.set_title('Interpolated to mesh025')

  ax3 = fig1.add_axes([0.2, 0.35, 0.6, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'interp_AVHRR_albedo_ithkn_mesh025.py'
  bottom_text(btx, pos=[0.08,0.3]) 



