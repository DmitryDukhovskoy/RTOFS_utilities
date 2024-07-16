"""
  Check NEP OBC UV vectors created from SPEAR monthly climatology
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
import matplotlib as mpl

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

klr    = 1  # model layer to plot
imonth = 1
nsegm  = 4  # 1,2,3,4

# Climatology derived for these years, started at mstart
yr1    = 1993
yr2    = 2020
mstart = 1

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)
spear_dir = os.path.join(config['filesystem']['spear_month_ens'], 'monthly_clim')
spear_topo_dir = config['filesystem']['spear_topo_grid']

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

# Load mapping indices gmapi:
dirgmapi = config['filesystem']['spear_mom_gmapi']
flgmaph  = f'spear2mom_NEP_OB_gmapi_hpnt.nc'
flgmapu  = f'spear2mom_NEP_OB_gmapi_upnt.nc'
flgmapv  = f'spear2mom_NEP_OB_gmapi_vpnt.nc'
dflgmaph = os.path.join(dirgmapi, flgmaph)
dflgmapu = os.path.join(dirgmapi, flgmapu)
dflgmapv = os.path.join(dirgmapi, flgmapv)
dsh = xarray.open_dataset(dflgmaph)

sfx   = f"segment_{nsegm:03d}"
uds   = f"u_{sfx}"
vds   = f"v_{sfx}"
dzvar = f"dz_u_{sfx}"

dset_segm = mutob.segm_topo(nsegm, HHM, hgrid)
U2d       = dsetOB[uds].isel(time=imonth-1).data.squeeze()  # 2D section
V2d       = dsetOB[vds].isel(time=imonth-1).data.squeeze()  # 2D section
segm_nm   = dset_segm['segm_name'].data[0]
Hbtm      = dset_segm['topo_segm'].data


if segm_nm == 'north':
  uscl = 700.
  xl1  = 180
  xl2  = idm+20
  yl1  = jdm-80
  yl2  = jdm+50
  jB   = jdm-1 
elif segm_nm == 'east':
  uscl = 700.
  xl1  = 300
  xl2  = idm+20 
  yl1  = 580
  yl2  = jdm 
  iB   = idm-1
elif segm_nm == 'south':
  uscl = 500.
  xl1  = -80
  xl2  = idm
  yl1  = -20
  yl2  = 100
  jB  = 0
elif segm_nm == 'west':
  uscl = 800.
  xl1  = -250
  xl2  = 120
  yl1  = 0
  yl2  = jdm 
  iB  = 0


cland = mpl.colormaps['PuBu_r']
cland.set_over(color=[0.3, 0.3, 0.3])
#clrmp.set_bad(color=[1,1,1])
rmin = -5000
rmax = 0.
plt.ion()
fgnmb = 1
fig = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
im  = ax1.pcolormesh(HHM, cmap = cland,
                     vmin=rmin, vmax=rmax)

ax1.contour(LONM, [x for x in range(150,290,5)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
ax1.contour(LATM, [x for x in range(5, 89, 5)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])

ax1.axis('scaled')

# Plot vectors
dstep = 20
varnmi = f"indx_segm{nsegm:03d}"
varnmj = f"jndx_segm{nsegm:03d}"
INDX   = dsh[varnmi].data
JNDX   = dsh[varnmj].data
lindx  = INDX.shape[0]
II     = INDX[:,0].squeeze()
JJ     = JNDX[:,0].squeeze()
for kk in range(1, lindx, dstep):
  k0 = klr-1
  kh = int(np.floor(kk/2.))  # closest h-point for topo
  if kh > len(Hbtm)-1: kh = leb(Hbtm)-1
  hh = Hbtm[kh]
  if hh >= 0.: continue
  uu = U2d[k0, kk]*uscl  # cm/ s
  vv = V2d[k0, kk]*uscl  # cm/s

  if segm_nm == 'north' or segm_nm == 'south':
    j0 = jB
    i0 = kk/2.
  else:
    j0 = kk/2.
    i0 = iB

  ax1.quiver(i0, j0, uu, vv, angles='xy', scale_units='xy', \
             scale=1, width=0.002, headlength=4, color=[0.,1.,0.8])


ax1.set_xlim([xl1, xl2])
ax1.set_ylim([yl1, yl2])

sttl = f"NEP OB: U z_lr={klr} M={imonth} OB={nsegm} {segm_nm}"
ax1.set_title(sttl)
ax1.set_xlabel('NEP I grid points')
ax1.set_ylabel('NEP J grid points')

btx = 'check_uvOBnep_vectors.py'
bottom_text(btx)




