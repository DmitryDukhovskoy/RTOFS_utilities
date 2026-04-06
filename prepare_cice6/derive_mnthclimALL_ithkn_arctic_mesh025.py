"""
  Derive monthly clim iithkn Arctic region
  combining:
    AVHRR monthly clim (2016-2025)
    combined CryoSat AWI mnthly clim (2015-2024 Jan-Apr, Oct-Dec) + NSIDC EASE 100 for summer months (May - Sept)
    ICESat-2: laser monthly summer (May-Aug 2019-2020/2021)
    
    Use weights for averaging with higher weights for longer time periods used for clim. 

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
import mod_regmom as mrmom
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)  

fld_name = 'ithkn'  # albedo or ithkn
regn = 'north'
save_final = True
box_fltr = True

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

def read_ice(MM, dflname, varnm='ice_thkn'):
  file_nm = os.path.basename(dflname)
  with xarray.open_dataset(dflname) as dsice:
    months = dsice['time'].data

    if MM not in months:
      print(f"MM={MM:02d} is not found in months in {file_nm}")
      return None

    A2d = dsice[varnm].sel(time=MM).data
  
  A2d[np.isnan(A2d)] = 0.
  return A2d

# Combined ice thickness clims:
# AVHRR monthly clim (2016-2025)
# CryoSat AWI mnthly clim (2015-2024 Jan-Apr, Oct-Dec) + NSIDC EASE 100 for summer months (May - Sept)
# ICESat-2: laser monthly summer (May-Aug 2019-2020/2021)

pthavhrr = os.path.join(pthdata,'AVHRR_albedo_ithkn','clim')
dflavhrr = os.path.join(pthavhrr,'AVHRR_ithkn_mnthclim_2016-2025_1440x1080_north.nc')

pthcryo = os.path.join(pthdata,'CryoSat_arctic_ice_snow_thkn','clim')
dflcryo = os.path.join(pthcryo,'ithkn_CryoSat_arcticAWI_mnthclim_2015-2024_1080x1440.nc')

pthices = os.path.join(pthdata,'ICESat2_arctic_summer_ithkn_hsnow','clim')
dflices = os.path.join(pthices,'ithkn_ICESat2_arctic_mnth_2019-2021_1080x1440.nc')

A3d = np.zeros((12,jdim,idim))
for imo in range(12):
  MM = imo + 1
  Aavhrr = read_ice(MM, dflavhrr, varnm='ice_thkn')
  Acryo  = read_ice(MM, dflcryo, varnm='ithkn')
  Aices  = read_ice(MM, dflices, varnm='ice_thkn')

  # For ice thickness - ignore 0s
  ok_avhrr = (Aavhrr > 0.)
  ok_cryo  = (Acryo > 0.)

  if Aices is None:
    nrec = ok_avhrr.astype(int) + ok_cryo.astype(int)
    A2d = np.divide( (Aavhrr + Acryo), nrec, out=np.zeros_like(Aavhrr), where=nrec != 0)

  else:
    # lower weight fo ICESat-2 because - shorter time series
    ok_ices = (Aices > 0.)
    # base weights
    wI = 0.5
    w0 = (1 - wI) / 2

    wt_avhrr = ok_avhrr.astype(float) * w0
    wt_cryo  = ok_cryo.astype(float)  * w0
    wt_ices  = ok_ices.astype(float) * wI

    wt_tot = wt_avhrr + wt_cryo + wt_ices

    # normalize weights
    wt_avhrr = np.divide(wt_avhrr, wt_tot, out=np.zeros_like(wt_tot), where=wt_tot != 0)
    wt_cryo  = np.divide(wt_cryo,  wt_tot, out=np.zeros_like(wt_tot), where=wt_tot != 0)
    wt_ices  = np.divide(wt_ices,  wt_tot, out=np.zeros_like(wt_tot), where=wt_tot != 0)

    A2d =  Aavhrr*wt_avhrr + Acryo*wt_cryo + Aices*wt_ices

  A2d[LMsk==0] = np.nan
  if box_fltr:
    jS = np.min(np.where(hlat >= 50)[0])
    jE = jdim-1
    nfltr = 1
    bx_sz = 27
    AAf = A2d.copy()
    for ifltr in range(nfltr):
      AAf = mrmom.box_averaging(A2d, HH, box_size = bx_sz, jS=jS, jE=jE, pole_wrap = True, LAT=hlat, LON=hlon)
    A2d = AAf.copy()

  A3d[imo,:,:] = A2d

if save_final:
  pthice   = os.path.join(pthdata, 'ithkn_clim_combined')
  fliceout = f'ithkn_mnthclim_cryo_avhrr_ices_{idim}x{jdim}_{regn}.nc'
  dfliceout = os.path.join(pthice,fliceout)

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

  dset_hi = xarray.Dataset({
    "ice_thkn": darr_hi,
    "lon": darr_lon,
    "lat": darr_lat,
  })
  dset_hi['ice_thkn'].attrs.update({
    "long_name": "sea ice thickness fuesed climatology",
    "units": "m",
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
    "title": f"Ice thickness climatology fused: AVHRR, CryoSat AWI + NSIDC EASE100 + ICESat-2 on mesh025 grid",
    "institution": "NOAA NWS MDC",
    "source": "derive_mnthclimALL_ithkn_arctic_mesh025.py",
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

  MM = 8
  A2d = A3d[MM-1,:,:].squeeze()

  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
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

  sttl = f'ithkn clim MM={MM:02d}\ncombined(AVHRR+CryoSat AWI + NSIDC EASE100 + ICESat-2)'
  ax1.set_title(sttl, fontsize=10)


  ax3 = fig1.add_axes([0.1, 0.06, 0.8, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'derive_mnthclimALL_ithkn_arctic_mesh025.py'
  bottom_text(btx, pos=[0.08,0.02], fsz=8) 



