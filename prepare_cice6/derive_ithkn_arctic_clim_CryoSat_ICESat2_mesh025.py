"""
  Derive monthly clim of ice thickn in the Arctic region
  combining
  CryoSat AWI monthly fields for winter months
  and NSIDC on EASE100 grid monthly fields for May - Sept +
  ICESat-2 laser altim. ice thikn (2019-2020) for May - Sept
  
  All data sets have been interpolated onto 025 mesh

"""

NOT FINISHED

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib  
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
import mod_cice6_utils as mc6util
import mod_regmom as mrmom
importlib.reload(mrmom)

regn = 'north'
box_fltr = True

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


A3d = np.zeros((12,jdm,idm))
for imo in range(1,13):
  print(f"Processing MM={imo:02d} ...")
  MM = imo
 
  if MM >=10 or MM <= 4:
    use_awi = True
  else:
    use_awi = False

  varnm = 'ice_thkn'
  if use_awi:
    pthinp, flinp = mc6util.pathfname_icesnow_mesh025(fyaml, node_nm, "ithkn_AWI")
  else:
    pthinp, flinp = mc6util.pathfname_icesnow_mesh025(fyaml, node_nm, "ithkn_NSIDC")

  dflclim_inp = os.path.join(pthinp, flinp)

  with xarray.open_dataset(dflclim_inp) as dsn:
    Time = dsn['time'].values
    irec = np.where(Time == MM)[0]
    assert irec.size > 0, f"Couldnot find month {MM} in data set"
    A2d = dsn[varnm].isel(time=irec).data.squeeze()

  # Cleanup missing data --> 0
  eps0 = 1.e-8
  AA = A2d.copy()
  AA[AA <= eps0] = 0.
  AA[np.isnan(AA)] = 0.
  AA[HH >= 0] = np.nan

  if box_fltr:
    jS = np.min(np.where(hlat >= 50)[0])
    jE = jdm-1
    nfltr = 1
    bx_sz = 27
    if imo >= 5 and imo <= 9:
      nfltr = 2
      bx_sz = 51
    AAf = AA.copy()
    for ifltr in range(nfltr):
      AAf = mrmom.box_averaging(AAf, HH, box_size = bx_sz, jS=jS, jE=jE, pole_wrap = True, LAT=hlat, LON=hlon)
    AA = AAf.copy()
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

dset_hs = xarray.Dataset({
    "ithkn": darr_hs,
    "lon": darr_lon,
    "lat": darr_lat,
  })
dset_hs['ithkn'].attrs.update({
  "long_name": "ice thickness, cell mean",
  "units": "m",
})

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
    "title": f"Ice thickness climatology from CryoSat AWI (winter) and NSIDC EASE100 (summer) interpolated onto mesh025 grid",
    "institution": "NOAA NWS NCEP MDC",
    "source": "derive_ithkn_arctic_monthclim_mesh025.py",
    "region": regn,
})


pthclim, fliceout = mc6util.pathfname_icesnow_mesh025(fyaml, node_nm, 'ithkn_AWI')
#fliceout = f'CryoSat_AWI_NSIDC_ithkn_mnthclim_mesh025_{jdm}x{idm}_{regn}.nc'

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
  rmax = 4.
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



