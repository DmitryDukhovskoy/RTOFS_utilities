"""
  Derive monthly clim of snow depth or ice thickn in the Arctic region
  using: 
  Interpolated NSIDC CryoSat snow or ice thickn. monthly fileds 2018-2021
  winter months only

  Warren (EWG Atlas) snow depth climatology - for summer months

  Both data sets have been interpolated onto 025 mesh

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
import mod_regmom as mrmom
importlib.reload(mrmom)

regn = 'north'
xtrp = 1   # extrapolate gaps in snow fields extending to coast or lat_xmin
lat_xmin = 65. # extend snow to this lat for ease of snow reconstruction in the IC files

parser = argparse.ArgumentParser()
parser.add_argument("--field", help=f"Field to interpolate: snow thkciness or ice thickn",
                    choices=['sndpth','ithkn'], required=True, type=str)
parser.add_argument("--xtrp", help=f"extrapolate snow fields to coast/min lat, default={xtrp}", 
                   choices=[0,1], type=int)
args = parser.parse_args()
  
field_name = args.field if args.field else None
xtrp_snow = xtrp == 1


syst_info = os.uname() 
machine = syst_info.nodename
  
if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
elif 'an' in machine:
  print("Running on PPAN node:", machine)
  node_nm = "ppan"
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

jdm, idm = HH.shape
LMsk = np.where(HH<0, 1, 0)

# For box-averaging:
jS = np.min(np.where(hlat >= 60)[0])
jE = jdm


A3d = np.zeros((12,jdm,idm))
for imo in range(1,13):
  print(f"Processing MM={imo:02d} ...")
  MM = imo
 
  if MM >= 10:
    yrS = 2018
    yrE = 2020
  else:
    yrS = 2019
    yrE = 2021

  if MM >=10 or MM <= 4:
    use_nsidc = True
  else:
    use_nsidc = False

  # Average monthly NSIDC fields:
  if use_nsidc:
    AA = None
    icc = 0
    for YR in range(yrS,yrE+1):

      # NSIDC snow/ice fields on mesh025:
      pthifld = os.path.join(pthdata,'CryoSat_arctic_ice_snow_thkn/interp_NSIDC_monthly')
      if field_name == 'sndpth':
        varnm = 'snow_depth'
        fliceout = f'hsnow_NSIDC_CryoSat_arctic_interp_mesh025_1080x1440_{YR}{MM:02d}.nc'
      elif field_name == 'ithkn':
        varnm = 'ice_thkn'
        fliceout = f'ithkn_NSIDC_CryoSat_arctic_interp_mesh025_1080x1440_{YR}{MM:02d}.nc'
  
      dflnsidc = os.path.join(pthifld,fliceout)

      with xarray.open_dataset(dflnsidc) as dsn:
        A2d = dsn[varnm].data.squeeze()

      if AA is None:
        AA = A2d.copy()
      else:
        AA += A2d
      icc += 1

    AA = AA / float(icc)

  else:
    if field_name == 'sndpth':
      varnm = 'snow_depth'
      flewg = 'hsnow_EWGatlas_interp_mesh025_1080x1440_north.nc'
    elif field_name == 'ithkn':
      varnm = None 
      fliceout = None  # not ready yet

    pthewg = os.path.join(pthdata,'Warren_snow_clim_EWG_atlas/snow_clim_interp')
    dflewg = os.path.join(pthewg, flewg)

    with xarray.open_dataset(dflewg) as dsn:
      AA = dsn[varnm].isel(time=imo-1).data.squeeze()

  # Cleanup missing data --> 0
  eps0 = 1.e-8
  hsnow_min = 0.01
  AA[AA <= eps0] = 0.
  mask_missed = (HH<0) & np.isnan(AA)
  AA[mask_missed] = 0.
  AA[HH>=0] = np.nan

  if xtrp_snow:
    AA[(AA > eps0) & (AA < hsnow_min)] = hsnow_min
    AAi = mrmom.extrapolate_to_lat_arctic(AA, hlon, hlat, HH, hlat0=65, Npnts=5, Rsearch=20., fill_land=False)
    #AA = mrmom.extrapolate_to_lat_arctic(AA, hlon, hlat, HH, hlat0=65, Npnts=5, Rsearch=20., fill_land=False)
    AAf = mrmom.box_averaging(AAi, HH, box_size = 15, jS=jS, jE=jE, pole_wrap = True)

  A3d[imo-1,:,:] = AA


jdim, idim = HH.shape
A3d  = A3d.astype('float32')
hlon = hlon.astype('float32')
hlat = hlat.astype('float32')
JD   = np.arange(jdim, dtype='int32')
ID   = np.arange(idim, dtype='int32')
time_months = np.arange(1,13, dtype='int32')
 
darr_hs = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
                   coords={"time": time_months,\
                           "jdim": JD,\
                           "idim": ID,})
darr_lon = xarray.DataArray(hlon, dims=("jdim","idim"),
                 coords={"jdim": JD,\
                         "idim": ID,})
darr_lat = xarray.DataArray(hlat, dims=("jdim","idim"),
                 coords={"jdim": JD,\
                         "idim": ID,})

if field_name == 'sndpth':
  dset_hs = xarray.Dataset({
    "snow_depth": darr_hs,
    "lon": darr_lon,
    "lat": darr_lat,
  })
  dset_hs['snow_depth'].attrs.update({
    "long_name": "snow depth on ice",
    "units": "m",
  })
  attr_str = 'Snow depth'
elif field_name == 'ithkn':
  dset_hs = xarray.Dataset({
    "ithkn": darr_hs,
    "lon": darr_lon,
    "lat": darr_lat,
  })
  dset_hs['ithkn'].attrs.update({
    "long_name": "ice thickness, cell mean",
    "units": "m",
  })
  attr_str = 'Ice thickness'

dset_hs['time'].attrs.update({
  "long_name": "months"
})
dset_hs['lon'].attrs.update({
  "long_name": "Longitudes",
  "units": "degrees_east",
})
dset_hs['lat'].attrs.update({
  "long_name": "Latitudes",
  "units": "degrees_north",
})

dset_hs.attrs.update({
    "title": f"{attr_str} climatology from CryoSat and EWG Atlas interpolated onto mesh025 grid",
    "institution": "NOAA NWS NCEP MDC",
    "source": "derive_hsnow_arctic_monthclim_mesh025.py",
    "region": regn,
})

if field_name == 'sndpth':
  fliceout = f'CryoSat_EWG_hsnow_mnthclim_mesh025_{idm}x{jdm}_{regn}.nc'
elif field_name == 'ithkn':
  fliceout = f'CryoSat_EWG_ithkn_mnthclim_mesh025_{idm}x{jdm}_{regn}.nc'

pthclim = os.path.join(pthdata,'CryoSat_arctic_ice_snow_thkn/clim')
dfliceout = os.path.join(pthclim,fliceout)
print(f'Dumping climatology fields  --> {dfliceout}')
dset_hs.to_netcdf(dfliceout,
        encoding={var: {'_FillValue': 1e30} for var in dset_hs.data_vars},
        format='NETCDF4')

f_chck = False
if f_chck:
  clrmp = mclrmps.colormap_temp(skip_dark=0.12)
  rmin = 0.
  rmax = 0.4
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  m = Basemap(projection='npstere',boundinglat=60,lon_0=-45,resolution='l')
  xh, yh = m(hlon,hlat) # GFS coords

  parallels = np.arange(40,89,10.)
  meridians = np.arange(-360,359.,45.)


  print("Plotting ...")

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.1, 0.15, 0.8, 0.8])
  m.drawparallels(parallels,labels=[0,0,0,0])
  # draw meridians
  meridians = np.arange(-360,359.,45.)
  m.drawmeridians(meridians,labels=[0,0,0,0])

  img = ax1.pcolormesh(xh,yh,AAf, cmap=clrmp, vmin=rmin, vmax=rmax)
  img = ax1.pcolormesh(xh,yh,AA, cmap=clrmp, vmin=rmin, vmax=rmax)
  #img = ax1.pcolormesh(xh,yh,AAi, cmap=clrmp, vmin=rmin, vmax=rmax)

  ax3 = fig1.add_axes([0.1, 0.08, 0.8, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)



