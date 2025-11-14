"""
  Derive monthly climatology of ice thickness in Antarctica region from
  Gridded estimates of Antarctic sea ice physical properties derived from 
  CryoSat-2 Baseline-D SAR and SARIn data spanning July 2010 through August 2021. 
  Data are generated using the CryoSat-2 Waveform-Fitting method for Antarctic sea ice (CS2WFA).

  Fons, S., Kurtz, N., & Bagnardi, M. (2022). 
  Antarctic Sea Ice Thickness Estimates from CryoSat-2: 2010-2021 (0.1.1) [Data set]. 
  Zenodo. https://doi.org/10.5281/zenodo.7327711

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
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6

regn = 'south'

yrE = 2020

parser = argparse.ArgumentParser()
parser.add_argument("--yrS", help="Start year to average, 2011,...,2020", type=int, required=True)
parser.add_argument("--yrE", help=f"End year to average, 2011,...,2020, default={yrE}", type=int)
args = parser.parse_args()
    
yrS    = args.yrS if args.yrS else None
yrE    = args.yrE if args.yrE else yrE
    
moS    = 1
moE    = 12
nmnths = moE-moS+1
 

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
pthice  = os.path.join(pthdata, 'CryoSat2_antarctic_ice_snow_thkn')

LON = LAT = None
H3D = None
kyrs = 0
NRECS = None  # count not nans records, they vary from year to year
for YR in range(yrS,yrE+1):
  kyrs += 1 

  for MM in range(moS, moE+1):
    print(f'Processing {YR}/{MM}')

    flice   = f"CS2WFA_25km_{YR}{MM:02d}.nc"
    dflice  = os.path.join(pthice,flice)
    print(f"Loading {dflice}")

    if LON is None: 
      with xarray.open_dataset(dflice) as dshi:
        LON = dshi['lon'].data.squeeze()
        LAT = dshi['lat'].data.squeeze()
      LON = np.where(LON > 180., LON-360., LON)

    # Note the NetCDF files have subgroups inside
    # Open the subgroup 'sea_ice_thickness' to read ice thikn.
    # other variables are in the root group:
    with xarray.open_dataset(dflice, group='sea_ice_thickness') as dshi:
      C2d = dshi['sea_ice_thickness'].data.squeeze()

    #C2d = np.where(np.isnan(C2d), 0., C2d)
    C2d = np.where(C2d < 0., 0., C2d)
    print(f'MM={YR}/{MM:02d}, Min/max ice thickn (m) = {np.nanmin(C2d):.1f}/{np.nanmax(C2d):.1f}')

    # Initialize arrays
    if H3D is None:
      jdim, idim = C2d.shape
      H3D = np.zeros((nmnths,jdim,idim))
      NRECS = np.zeros((nmnths, jdim, idim), dtype=int)
      DMsk = np.zeros((jdim,idim), dtype=bool)  # accumulate data points across all yrs, months

    data_pnts = ~np.isnan(C2d)   # valid data points
    DMsk |= data_pnts  # either DMsk or datapnts is True

    # Sum of N years by months:
    imo = MM-moS
    vals = np.where(data_pnts, C2d, 0.0)
    H3D[imo,:] += vals
    NRECS[imo,:] += data_pnts.astype(int)

# Permanent NaNs:
NAN_mask = ~DMsk 

# Average:
H3D_mean = np.divide(H3D, NRECS, out=np.full_like(H3D, np.nan), where=NRECS > 0)
H3D_mean[:, NAN_mask] = np.nan


# Save netcdf:
ntm, jdim, idim = H3D_mean.shape
  
# Double --> single precision
H3D_mean = H3D_mean.astype('float32')
LON = LON.astype('float32')
LAT = LAT.astype('float32')

Xgrid = np.arange(idim)
Ygrid = np.arange(jdim)

# Save to netCDF:
time_months = np.arange(1,13, dtype='int32')
darr_hi = xarray.DataArray(H3D_mean, dims=("time","y","x"),\
                 coords={"time": time_months,\
                         "y": Ygrid,\
                         "x": Xgrid})
darr_xx = xarray.DataArray(Xgrid, dims=("x"),
                 coords={"x": Xgrid})
darr_yy = xarray.DataArray(Ygrid, dims=("y"),
                 coords={"y": Ygrid})
darr_lon = xarray.DataArray(LON, dims=("y","x"),
                 coords={"y": Ygrid,\
                         "x": Xgrid})
darr_lat = xarray.DataArray(LAT, dims=("y","x"),
                 coords={"y": Ygrid,\
                         "x": Xgrid})
dset_hi = xarray.Dataset(
    {
        "ice_thickness": darr_hi,
        "lon": darr_lon,
        "lat": darr_lat,
        "x": darr_xx,
        "y": darr_yy,
    }
)

dset_hi['time'].attrs.update({
  "long_name": "months"
})
dset_hi['ice_thickness'].attrs.update({
  "long_name": "mean sea ice thickness in grid cell",
  "units": "m",
})
dset_hi['lon'].attrs.update({
  "long_name": "Longitudes",
  "units": "degrees",
})
dset_hi['lat'].attrs.update({
  "long_name": "Latitudes",
  "units": "degrees",
})

dset_hi.attrs.update({
  "title": f"Antarctic Sea Ice Thickness Estimates from CryoSat-2 monthly clim",
  "info": f"Monthly climatology fields on native grid years: {yrS}-{yrE}",
  "info2": "https://zenodo.org/records/7327711",
  "institution": "NOAA NWS NCEP MDC",
  "source": "derive_ithkn_clim_CryoSat_antarct.py",
  "contact": "dmitry.dukhovskoy@noaa.gov",
  "region": regn,
})


pthdump = os.path.join(pthice,'clim')
fhice  = f'CryoSat_ithkn_mnthly_clim_{idim}x{jdim}_{regn}.nc'
dfhice = os.path.join(pthdump, fhice)

print(f'Saving ice thickn climatology --> {dfhice}')
dset_hi.to_netcdf(dfhice, 
        encoding={var: {'_FillValue': 1e30} for var in dset_hi.data_vars},
        format='NETCDF4', engine='netcdf4')


f_chck = False
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:

  #clrmp = mclrmps.colormap_uv()
  #rmin = -1.
  #rmax = 1.

  #clrmp = mclrmps.colormap_conc()
  #rmin = 0.
  #rmax = 1.

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

  MM = 9
  C2d = H3D_mean[MM-1,:,:].squeeze()

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.12, 0.12, 0.8, 0.8])

  # Plot on grid:
  img = ax1.pcolormesh(C2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  # Check longitudes:
  cs = ax1.contour(LON1, lon_cntr1, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
  ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
  cs2 = ax1.contour(LON2, lon_cntr2, linestyles='solid', colors=[(0.,0.5,0.9)], linewidths=1)
  ax1.clabel(cs2, inline=True, fontsize=10, fmt="%.1f")
  cs3 = ax1.contour(LON1,[0], linestyles='solid', colors=[(1,0.,0.)], linewidths=2)
  ax1.clabel(cs3, inline=True, fontsize=12, fmt="%.1f")
  ax1.contour(LON3,[180], linestyles='solid', colors=[(0.6,0.,0.9)], linewidths=1)
  ax1.contour(LAT, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
  ax1.contour(LAT,[-75], linestyles='solid', colors=[(1,0.,0.)], linewidths=1)  
  ax1.axis('scaled')
  ax1.invert_yaxis() 
  ax1.set_ylabel('Inverted j index')
  ax1.set_ylabel('j index')
  ax1.set_xlabel('i index')

  ax1.set_title(f'CryoSat ice thickness {YR}/{MM:02d}')

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

  btx = 'derive_ithkn_clim_CryoSat_antarct.py'
  bottom_text(btx) 





