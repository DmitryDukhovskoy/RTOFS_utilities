"""
  Check NEP OBC UV fields created from SPEAR monthly climatology
  Check 2D scalar fields (vertical sections along OBs)
  see: write_spearOB_climatology.py
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

imonth = 1
nsegm  = 3

# Climatology derived for these years, started at mstart
yr1    = 1993
yr2    = 2020
mstart = 1

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

pthoutp = gridfls['MOM6_NEP']['seasonal_fcst']['pthoutp']
fobc_out = os.path.join(pthoutp,f'OBCs_spear_clim_{yr1}-{yr2}_mstart{mstart:02d}.nc')
print(f'Loading OBCs <--- {fobc_out}')
dsetOB = xarray.open_dataset(fobc_out)

sfx   = f"segment_{nsegm:03d}"
uds   = f"u_{sfx}"
vds   = f"v_{sfx}"
dzvar = f"dz_u_{sfx}"

dset_segm = mutob.segm_topo(nsegm, HHM, hgrid)

U2d    = dsetOB[uds].isel(time=imonth-1).data.squeeze()  # 2D section
V2d    = dsetOB[vds].isel(time=imonth-1).data.squeeze()  # 2D section
dZ     = dsetOB[dzvar].isel(time=imonth-1).data.squeeze()  
ZZ, ZM = mmom6.zz_zm_fromDZ(dZ)

S2d = np.sqrt(U2d**2 + V2d**2)

#rmin, rmax = mutob.minmax_clrmap(S2d, cpnt=0)

clrmp = mutil.colormap_temp(clr_ramp=[0.95,0.92,1])
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
  rmin = 0.
  rmax = 0.1
elif segm_nm == 'east':
  xl1 = 5800
  rmin = 0.
  rmax = 0.1
elif segm_nm == 'south':
  xl1 = 0
  rmin = 0.
  rmax = 0.2
elif segm_nm == 'west':
  xl1 = 0
  rmin = 0.
  rmax = 0.2
else:
  xl1 = 0

xl2  = max(Xbtm)

sttl = f"NEP OB: |U| M={imonth} OB={segm_nm}"
mutob.plot_xsection(S2d, X, ZM, Hbtm, Xbtm, clrmp, rmin, rmax, xl1, xl2, sttl=sttl, stxt=stxt, fgnmb=1)

btx = 'check_uvOBnep.py'
bottom_text(btx)




