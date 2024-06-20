"""
  Use subset of SPEAR monthly climatology fields
  derived in derive_monthly_clim_spear.py

  boundary/Segments functions are based  on Andrew C. Ross 

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
import mod_mom6 as mom6util
#importlib.reload(mom6util)

indir = Path('/work/acr/spear/processed/ensmean')
outdir = Path('/work/acr/spear/climatology')

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

# MOM6 NEP topo/grid:
pthtopo = gridfls['MOM6_NEP']['test']['pthgrid']
fgrid   = gridfls['MOM6_NEP']['test']['fgrid']
hgrid = xarray.open_dataset(os.path.join(pthtopo,fgrid))

segments = [ Segment(1, 'north', hgrid, output_dir=outdir),
             Segment(2, 'east',  hgrid, output_dir=outdir),
             Segment(3, 'south', hgrid, output_dir=outdir),
             Segment(4, 'west',  hgrid, output_dir=outdir)]

nOB = len(segments)

# segments[0].__dict__

static = xarray.open_dataset('/work/acr/spear/analysis/ocean_z.static.nc')

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
  ds        = mutob.load_var(spear_dir, 'thetao')
  lon_spear = ds.xh.data
  lat_spear = ds.yh.data
  mutob.spear2momOB_gmapi(segments, lon_spear, lat_spear, dflgmaph)
dsh = xarray.open_dataset(dflgmaph)

# u-point indices
if not os.path.isfile(dflgmapu):
  print(f'Mapping indices upnt are missing, {dflgmapu}, deriving ...')
  ds        = mutob.load_var(spear_dir, 'uo')
  lon_spear = ds.xq.data
  lat_spear = ds.yh.data
  mutob.spear2momOB_gmapi(segments, lon_spear, lat_spear, dflgmapu)
dsu = xarray.open_dataset(dflgmapu)

# v-point indices
if not os.path.isfile(dflgmapv):
  print(f'Mapping indices vpnt are missing, {dflgmapv}, deriving ...')
  ds        = mutob.load_var(spear_dir, 'vo')
  lon_spear = ds.xh.data
  lat_spear = ds.yq.data
  mutob.spear2momOB_gmapi(segments, lon_spear, lat_spear, dflgmapv)
dsv = xarray.open_dataset(dflgmapv)

icc = 0
dsetOB = xarray.Dataset()
for varnm in ['thetao', 'so']:
# Load monthly climaotology SPEAR data subset for NEP
  ds = mutob.load_var(spear_dir, varnm)

  for isgm in range(nOB):
    nsgm   = isgm+1
    print(f'Processing {varnm} OB segment={nsgm}')
    INDX   = dsh[f'indx_segm{nsgm:03d}'].data
    JNDX   = dsh[f'jndx_segm{nsgm:03d}'].data
    dset   = mutob.derive_obsegment(hgrid, ds, segments, isgm, varnm, INDX, JNDX)
#    dset   = xarray.Dataset({f"{varnm}_segment_{nsgm:03d}": darr})
    dsetOB = xarray.merge([dsetOB, dset])

# Ssh - 1D sections

# UV fields


  climo = compute_climatology(ds)
  climo = add_coords(climo)
  # this wont work if using daily SSH
  for mon in [3, 6, 9, 12]:
    print(f'  month {mon:02d}')
    mon_data = climo.sel(month=mon)
    mon_data = add_valid_time(mon_data)
    for seg in segments:
      print('    ' + seg.segstr)
      regridded = seg.regrid_tracer(
          mon_data,
          regrid_suffix='spear_tracer',
          write=False
      )
      seg.to_netcdf(regridded, f'spear_{varnm}_i{mon:02d}_climo')

print('uv')
u_climo = add_coords(compute_climatology(load_var('uo')))
v_climo = add_coords(compute_climatology(load_var('vo')))
for mon in [3, 6, 9, 12]:
    print(f'  month {mon:02d}')
    u_mon = add_valid_time(u_climo.sel(month=mon))
    v_mon = add_valid_time(v_climo.sel(month=mon))
    for seg in segments:
        print('    ' + seg.segstr)
        regridded = seg.regrid_velocity(
            u_mon, v_mon,
            regrid_suffix='spear_velocity',
            write=False
        )
        seg.to_netcdf(regridded, f'spear_uv_i{mon:02d}_climo')

