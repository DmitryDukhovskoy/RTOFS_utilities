"""
  Use subset of SPEAR monthly climatology fields
  derived in derive_monthly_clim_spear.py

  Note SPEAR U,V are on MOM staggerd C-grid
 
  OB segments are on the supergrid
"""

import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
#import pickle
import matplotlib.pyplot as plt
from yaml import safe_load

import mod_utils_ob as mutob
importlib.reload(mutob)


PPTHN = []
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')
sys.path.append('./seasonal-workflow')
from boundary import Segment
import mod_time as mtime
import mod_mom6 as mmom6

# Climatology derived for these years, started at mstart
yr1    = 1993
yr2    = 2020
mstart = 1

indir = Path('/work/acr/spear/processed/ensmean')
outdir = Path('/work/acr/spear/climatology')

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

# MOM6 NEP topo/grid:
pthtopo    = gridfls['MOM6_NEP']['test']['pthgrid']
fgrid      = gridfls['MOM6_NEP']['test']['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"]["test"]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc')) 
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom  = os.path.join(pthtopo, fgrid)
# Hgrid lon. lat:
hlon, hlat  = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')


segments = [ Segment(1, 'north', hgrid, output_dir=outdir),
             Segment(2, 'east',  hgrid, output_dir=outdir),
             Segment(3, 'south', hgrid, output_dir=outdir),
             Segment(4, 'west',  hgrid, output_dir=outdir)]

nOB = len(segments)

# Time array for monthly clim:
time_days = np.zeros((12))
for mo in range(1,13):
  jday = mtime.date2jday([yr1,mo,15])
  time_days[mo-1] = jday

# segments[0].__dict__

#static = xarray.open_dataset('/work/acr/spear/analysis/ocean_z.static.nc')
grid_spear = xarray.open_dataset('/work/Dmitry.Dukhovskoy/data/SPEAR/ocean_z.static.nc')
icegrid_spear = xarray.open_dataset('/work/Dmitry.Dukhovskoy/data/SPEAR/ice.static.nc')

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)
spear_dir = os.path.join(config['filesystem']['spear_month_ens'], 'monthly_clim')

# Check if mapping indices exist, gmapi:
dirgmapi = config['filesystem']['spear_mom_gmapi']
flgmaph  = f'spear2mom_NEP_OB_gmapi_hpnt.nc'
flgmapu  = f'spear2mom_NEP_OB_gmapi_upnt.nc'
flgmapv  = f'spear2mom_NEP_OB_gmapi_vpnt.nc'
dflgmaph = os.path.join(dirgmapi, flgmaph)
dflgmapu = os.path.join(dirgmapi, flgmapu)
dflgmapv = os.path.join(dirgmapi, flgmapv)
# h-point indices
if not os.path.isfile(dflgmaph):
  print(f'Mapping indices hpnt are missing, {dflgmaph}, deriving ...')
  ds        = mutob.load_var(spear_dir, 'thetao', yr1=yr1, yr2=yr2, mstart=mstart)
  lon_spear, lat_spear = mutob.subset_spear_coord('h')
  mutob.spear2momOB_gmapi(segments, lon_spear, lat_spear, dflgmaph)
dsh = xarray.open_dataset(dflgmaph)

# u-point indices
if not os.path.isfile(dflgmapu):
  print(f'Mapping indices upnt are missing, {dflgmapu}, deriving ...')
  ds        = mutob.load_var(spear_dir, 'uo')
  lon_spear, lat_spear = mutob.subset_spear_coord('u')
  mutob.spear2momOB_gmapi(segments, lon_spear, lat_spear, dflgmapu)
dsu = xarray.open_dataset(dflgmapu)

# v-point indices
if not os.path.isfile(dflgmapv):
  print(f'Mapping indices vpnt are missing, {dflgmapv}, deriving ...')
  ds        = mutob.load_var(spear_dir, 'vo')
  lon_spear, lat_spear = mutob.subset_spear_coord('v')
  mutob.spear2momOB_gmapi(segments, lon_spear, lat_spear, dflgmapv)
dsv = xarray.open_dataset(dflgmapv)

icc = 0
dsetOB = xarray.Dataset()
for isgm in range(nOB):
  nsgm = isgm+1
  print(f'Processing lon/lat OB segment={nsgm}')
  dset   = mutob.derive_obsegm_lonlat(hgrid, segments, isgm)
  dsetOB = xarray.merge([dsetOB, dset])

for varnm in ['thetao', 'so']:
# Load monthly climatology SPEAR data subset for NEP
  ds = mutob.load_var(spear_dir, varnm, fzint=True)

  for isgm in range(nOB):
    nsgm   = isgm+1
    print(f'Processing {varnm} OB segment={nsgm}')
    INDX   = dsh[f'indx_segm{nsgm:03d}'].data
    JNDX   = dsh[f'jndx_segm{nsgm:03d}'].data
    dset   = mutob.derive_obsegm_3D(hgrid, ds, segments, isgm, varnm, 
                                    INDX, JNDX, time_steps=time_days)
#    dset   = xarray.Dataset({f"{varnm}_segment_{nsgm:03d}": darr})
    dsetOB = xarray.merge([dsetOB, dset])


# =============================
#  Plot OB points for checking
# =============================
f_chckOB = False
if f_chckOB: 
  varplt = 'u'
  nsgm = 2
  dset_segm = segments[nsgm-1]
  xOB   = dset_segm.coords.lon.data
  yOB   = dset_segm.coords.lat.data
  if varplt == 'h':
    INDX   = dsh[f'indx_segm{nsgm:03d}'].data
    JNDX   = dsh[f'jndx_segm{nsgm:03d}'].data
  elif varplt == 'u':
    INDX   = dsu[f'indx_segm{nsgm:03d}'].data
    JNDX   = dsu[f'jndx_segm{nsgm:03d}'].data
  elif varplt == 'v':
    INDX   = dsv[f'indx_segm{nsgm:03d}'].data
    JNDX   = dsv[f'jndx_segm{nsgm:03d}'].data
  HHM = dstopo_nep['depth'].data
  HHM = np.where(HHM < 1.e-20, np.nan, HHM)
  HHM = -HHM
  HHM = np.where(np.isnan(HHM), 1., HHM)

# SPEAR depths
  HHS, lon_topo, lat_topo = mutob.subset_spear_topo()
  LONs, LATs = subset_spear_coord(varplt)
  xOB = np.where(xOB<0., xOB+360., xOB)
  LONs = np.where(LONs<0., LONs+360., LONs)
  lon_topo = np.where(lon_topo<0., lon_topo+360., lon_topo)

  mutob.plot_OBpoints(xOB, yOB, INDX, JNDX, LONs, LATs, lon_topo, lat_topo, \
                      hlon, hlat, HHM, HHS, fgnmb=1)

#a=-1
#assert a==1

# Ssh - 1D sections
# Load monthly climatology SPEAR data subset for NEP
ds = mutob.load_var(spear_dir, 'ssh', fzint=False)

for isgm in range(nOB):
  nsgm  = isgm+1
  print(f'Processing ssh OB segment={nsgm}')
  INDX   = dsh[f'indx_segm{nsgm:03d}'].data
  JNDX   = dsh[f'jndx_segm{nsgm:03d}'].data
  dset   = mutob.derive_obsegm_ssh(hgrid, ds, segments, isgm, INDX, JNDX, time_steps=time_days)
#    dset   = xarray.Dataset({f"{varnm}_segment_{nsgm:03d}": darr})
  dsetOB = xarray.merge([dsetOB, dset])


# Plot rot angles:
f_plt = False
if f_plt: mutob.plot_rotangle(grid_spear, icegrid_spear)
f_pltNEP = False
if f_pltNEP: mutob.plot_rotangleNEP(hgrid, hmask)

# UV fields
# Derive rotation angle, rad and components of the rotation matrix
# for rotating vectors onto true N/E grid from SPEAR
r2d = 180./np.pi
theta_rot, cosrot, sinrot = mutob.get_rotangle(icegrid_spear, fconfig, grid_spear)
print(f"Rotation angle for SPEAR min/max: {np.min(theta_rot)*r2d:6.2f} " + \
      f"/ {np.max(theta_rot)*r2d:6.2f}")

# Load monthly climatology SPEAR data subset for NEP
ds_uo = mutob.load_var(spear_dir, 'uo', fzint=True)
ds_vo = mutob.load_var(spear_dir, 'vo', fzint=True)

for varnm in ['u', 'v']:
  for isgm in range(nOB):
    nsgm   = isgm+1
    print(f'Processing {varnm} OB segment={nsgm}')

    if varnm == 'u':
      INDX  = dsu[f'indx_segm{nsgm:03d}'].data
      JNDX  = dsu[f'jndx_segm{nsgm:03d}'].data
    elif varnm == 'v':
      INDX  = dsv[f'indx_segm{nsgm:03d}'].data
      JNDX  = dsv[f'jndx_segm{nsgm:03d}'].data

    dset   = mutob.derive_obsegm_uv(hgrid, ds_uo, ds_vo, segments, isgm, theta_rot,\
                                    varnm, INDX, JNDX, time_steps=time_days)
    dsetOB = xarray.merge([dsetOB, dset])

for varnm in ['thetao','so','u','v']:
  for segm in [1,2,3,4]:
    vv  = f"{varnm}_segment_{segm:03d}"
    vdz = f"dz_{varnm}_segment_{segm:03d}"
    dm1 = f'lat_segment_{segm:03d}'
    dm2 = f'lon_segment_{segm:03d}'
    dsetOB[vv].attrs["coordinates"] = f"{dm1} {dm2}"
    dsetOB[vdz].attrs["coordinates"] = f"{dm1} {dm2}"

varnm='zos'
for segm in [1,2,3,4]:
  vv  = f"{varnm}_segment_{segm:03d}"
  dm1 = f'lat_segment_{segm:03d}'
  dm2 = f'lon_segment_{segm:03d}'
  dsetOB[vv].attrs["coordinates"] = f"{dm1} {dm2}"

dsetOB.attrs["history"] = f"Created from SPEAR monthly clim {yr1}-{yr2}"
dsetOB.attrs["code"] = f"/home/Dmitry.Dukhovskoy/python/setup_seasonal_NEP/write_spearOB_climatology.py"

"""
  Add attributes for time var for climatology OBCs
"""
#dsetOB['time'] = np.arange(0, ntsteps, dtype='float')
dsetOB['time'].attrs['units'] = 'days since 0001-01-01'
dsetOB['time'].attrs['calendar'] = 'noleap'
dsetOB['time'].attrs['modulo'] = ' '
dsetOB['time'].attrs['cartesian_axis'] = 'T'


encode = {
  'time': dict(dtype='float64', _FillValue=1.0e20)
}
for varnm in ['lon','lat']:
  for segm in [1, 2, 3, 4]:
    encode.update({
      f'{varnm}_segment_{segm:03d}': dict(dtype='float64', _FillValue=1.0e20)
    })
for varnm in ['zos']:
  for segm in [1, 2, 3, 4]:
    encode.update({
    f'{varnm}_segment_{segm:03d}': dict(_FillValue=1.0e20)
    })
for varnm in ['thetao', 'so', 'u', 'v']:
  for segm in [1, 2, 3, 4]:
    encode.update({
    f'{varnm}_segment_{segm:03d}': dict(_FillValue=1.0e20),
    f'dz_{varnm}_segment_{segm:03d}': dict(_FillValue=1.0e20)
    })

pthoutp = gridfls['MOM6_NEP']['seasonal_fcst']['pthoutp']
fobc_out = os.path.join(pthoutp,f'OBCs_spear_clim_{yr1}-{yr2}_mstart{mstart:02d}.nc')
print(f'Saving OBCs ---> {fobc_out}')
dsetOB.to_netcdf(fobc_out, 
                 format='NETCDF3_64BIT', 
                 engine='netcdf4', 
                 encoding=encode, 
                 unlimited_dims='time')


