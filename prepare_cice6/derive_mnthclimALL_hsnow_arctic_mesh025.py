"""
  Derive monthly clim hsnow Arctic region
  combining:
    CryoSat AWI + EWG for summer
    ICESat-2: laser monthly summer (May-Aug 2019-2020/2021) - snow estimates 
    from Snow Models with ERA5, MERRA-2 forcing and Warren 1999 clim. 
    
    Use weights for averaging summer months

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

fld_name = 'hsnow'  # albedo or hsnow
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

# Combined snow depth clims:
# ICESat-2: model-based snow depth estimates and Warren 99 clim
#  
pthcryo = os.path.join(pthdata,'CryoSat_arctic_ice_snow_thkn','clim')
dflcryo = os.path.join(pthcryo,'CryoSat_EWG_hsnow_mnthclim_mesh025_1440x1080_north.nc')

pthices = os.path.join(pthdata,'ICESat2_arctic_summer_ithkn_hsnow','clim')
dflsme5 = os.path.join(pthices,'hsnow_sm_e5_ICESat2_arctic_mnth_2019-2021_1080x1440.nc')
dflsmm2 = os.path.join(pthices,'hsnow_sm_m2_ICESat2_arctic_mnth_2019-2021_1080x1440.nc')
dflw99r = os.path.join(pthices,'hsnow_w99r_ICESat2_arctic_mnth_2019-2021_1080x1440.nc')

A3d = np.zeros((12,jdim,idim))
for imo in range(12):
  MM = imo + 1
  Acryo  = read_ice(MM, dflcryo, varnm='snow_depth')
  Asme5  = read_ice(MM, dflsme5, varnm='snow_depth')
  Asmm2  = read_ice(MM, dflsmm2, varnm='snow_depth')
  Aw99r  = read_ice(MM, dflw99r, varnm='snow_depth')

  if Asme5 is None or Asmm2 is None or Aw99r is None:
    A2d = Acryo.copy()
  else:
    # For summer months only:
    # For snow depths - ignore 0s to avoid averaging ice / no ice grid points
    ok_cryo  = (Acryo > 0.)
    ok_sme5  = (Asme5 > 0.)
    ok_smm2  = (Asmm2 > 0.)
    ok_w99r  = (Aw99r > 0.)

    # change weights for Model-based estimates and Warn. clim 
    # base weights
    #wCryo = 0.4
    #wW99  = 0.3
    #wE5   = 0.15
    #wS2   = 0.15
    wCryo = 0.98
    wW99  = 0.01
    wE5   = 0.0
    wS2   = 0.01


    assert abs(1 - (wCryo + wW99 + wE5 + wS2)) < 1.e-6, "ERR: Check base weights"

    wt_cryo  = ok_cryo.astype(float) * wCryo
    wt_sme5  = ok_sme5.astype(float) * wE5
    wt_smm2  = ok_smm2.astype(float) * wS2
    wt_w99r  = ok_w99r.astype(float) * wW99

    wt_tot = wt_cryo + wt_sme5 + wt_smm2 + wt_w99r

    # normalize weights
    wt_cryo  = np.divide(wt_cryo,  wt_tot, out=np.zeros_like(wt_tot), where=wt_tot != 0)
    wt_sme5  = np.divide(wt_sme5,  wt_tot, out=np.zeros_like(wt_tot), where=wt_tot != 0)
    wt_smm2  = np.divide(wt_smm2,  wt_tot, out=np.zeros_like(wt_tot), where=wt_tot != 0)
    wt_w99r  = np.divide(wt_w99r,  wt_tot, out=np.zeros_like(wt_tot), where=wt_tot != 0)

    A2d = Acryo*wt_cryo + Asme5*wt_sme5 + Asmm2*wt_smm2 + Aw99r*wt_w99r

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
  pthice   = os.path.join(pthdata, 'hsnow_clim_combined')
  fliceout = f'hsnow_mnthclim_Cryo_ICESat2_{idim}x{jdim}_{regn}.nc'
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
    "snow_depth": darr_hi,
    "lon": darr_lon,
    "lat": darr_lat,
  })
  dset_hi['snow_depth'].attrs.update({
    "long_name": "snow depth fused climatology",
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
    "title": f"Snow depth climatology fused: CryoSat AWI + NSIWDC EWG + ICESat-2: SnowModel ERA5, MERRA2, Warren99",
    "institution": "NOAA NWS OMD",
    "source": "derive_mnthclimALL_hsnow_arctic_mesh025.py",
    "region": regn,
  })

  print(f'Dumping interpolated {fld_name} --> {dfliceout}')
  dset_hi.to_netcdf(dfliceout,
          encoding={var: {'_FillValue': np.float32(1e30)} for var in dset_hi.data_vars},
          format='NETCDF4')


f_chck = True
f_xy = False     # True - plot on X.Y grid. False - plot on index space
if f_chck:
  MM = 7
  A2d = A3d[MM-1,:,:].squeeze()

  clrmp = mclrmps.colormap_temp()
  rmin = 0.
  rmax = 0.2
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

  sttl = f'hsnow clim MM={MM:02d}\nfused: CryoSat + NSIDC EWG(summer) + ICESat-2(summer): sme5, smm2, w99r)'
  ax1.set_title(sttl, fontsize=10)


  ax3 = fig1.add_axes([0.1, 0.06, 0.8, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'derive_mnthclimALL_hsnow_arctic_mesh025.py'
  bottom_text(btx, pos=[0.08,0.02], fsz=8) 



