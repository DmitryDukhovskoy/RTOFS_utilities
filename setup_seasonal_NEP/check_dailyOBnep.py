"""
  Check NEP OBC fields created from SPEAR monthly fields
  and interpoalted into daily data
  Check 2D scalar fields (vertical sections along OBs)
  see: create_daily_from_monthly_spear.py
"""
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
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
from mod_utils_fig import bottom_text
import mod_utils as mutil
importlib.reload(mutil)


varnm  = 'thetao' # so thetao uspd
iday   = 1  # f/cast day sarting from 1, ...., 365(6)
nsegm  = 3  # OB segments, 1 - North, 2 - East, ...
run_name = 'seasonal_fcst_daily'

# OB derived for f.cast initialized on these yr/mo
ens_spear = 1       # ens run used to create OB
yr_start = 1993
mo_start = 4
dd_start = 1
dnmb_start = mtime.datenum([yr_start, mo_start, dd_start])
dv_start = mtime.datevec(dnmb_start)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

# MOM6 NEP topo/grid:
pthtopo     = gridfls['MOM6_NEP']['test']['pthgrid']
fgrid_mom   = gridfls['MOM6_NEP']['test']['fgrid']
ftopo_mom   = gridfls["MOM6_NEP"]["test"]["ftopo"]
dfgrid_mom  = os.path.join(pthtopo, fgrid_mom)
dftopo_mom  = os.path.join(pthtopo, ftopo_mom)
LONM, LATM  = mmom6.read_mom6grid(dfgrid_mom)
HHM         = mmom6.read_mom6depth(dftopo_mom)
jdm, idm    = np.shape(HHM)
hgrid       = xarray.open_dataset(os.path.join(pthtopo,fgrid_mom))
hmask       = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))

date_init = f'{dv_start[0]}{dv_start[1]:02d}{dv_start[2]:02d}'
pthoutp = gridfls['MOM6_NEP'][run_name]['pthoutp']
fobc_out = os.path.join(pthoutp,f'OBCs_spear_daily_init{date_init}_e{ens_spear:02d}.nc')
#fobc_out = os.path.join(pthoutp,f'OBCs_spear_daily_init{date_init}.nc')
print(f'Loading OBCs <--- {fobc_out}')
dsetOB = xarray.open_dataset(fobc_out)

sfx   = f"segment_{nsegm:03d}"
vards = f"{varnm}_{sfx}"
dzvar = f"dz_{varnm}_{sfx}"

dset_segm = mutob.segm_topo(nsegm, HHM, hgrid)

A2d    = dsetOB[vards].isel(time=iday-1).data.squeeze()  # 2D section
dZ     = dsetOB[dzvar].isel(time=iday-1).data.squeeze()  
ZZ, ZM = mmom6.zz_zm_fromDZ(dZ)


rmin, rmax = mutob.minmax_clrmap(A2d, cpnt=0)


if varnm == 'so':
  clrmp = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
elif varnm == 'thetao':
  clrmp = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[1,1,1])


lon1     = dset_segm['lon_segm'].data[0]
lon2     = dset_segm['lon_segm'].data[-1]
lat1     = dset_segm['lat_segm'].data[0]
lat2     = dset_segm['lat_segm'].data[-1]
X        = dset_segm['dist_supergrid'].data
Xbtm     = dset_segm['dist_grid'].data
Hbtm     = dset_segm['topo_segm'].data
Hbtm     = np.where(Hbtm > 0, 0., Hbtm)
segm_nm  = dset_segm['segm_name'].data[0]

stxt = f"X-axis: Distance (km) along segment from {lon1:5.2f}W, {lat1:5.2f}N to {lon2:5.2f}W, {lat2:5.2f}N"

Xbtm[-1] = X[-1]
if segm_nm == 'north':
  xl1 = 1900.
  if varnm == 'thetao':
    rmin = -2.
    rmax = 1.
  if varnm == 'so':
    rmin = 29.
    rmax = 35.
elif segm_nm == 'east':
  xl1 = 5800
  if varnm == 'thetao':
    rmin = -2.
    rmax = 2.
  if varnm == 'so':
    rmin = 29.
    rmax = 35.
elif segm_nm == 'south':
  xl1 = 0
  if varnm == 'so':
    rmin = 34.0
    rmax = 34.9
  else:
    rmin = 0.
    rmax = 28.
elif segm_nm == 'west':
  xl1 = 0
  if varnm == 'so':
    rmin = 34.0
    rmax = 34.9
  else:
    rmin = 0.
    rmax = 25.0
else:
  xl1 = 0

xl2  = max(Xbtm)

sttl = f"NEP OB: {varnm} fcast day={iday} OB={segm_nm}"
mutob.plot_xsection(A2d, X, ZM, Hbtm, Xbtm, clrmp, rmin, rmax, xl1, xl2, sttl=sttl, stxt=stxt, fgnmb=1)

btx = 'check_dailyOBnep.py'
bottom_text(btx)




