"""
  Plot SPEAR and NEP bathmetries
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
import mod_utils as mutil
importlib.reload(mutil)
from boundary import Segment
import mod_time as mtime
import mod_mom6 as mmom6
from mod_utils_fig import bottom_text
import mod_misc1 as mmisc


fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)
spear_dir = os.path.join(config['filesystem']['spear_month_ens'], 'monthly_clim')
spear_topo_dir = config['filesystem']['spear_topo_grid']

# SPEAR topo, subset to NEP domain similar to data files:
ftopo_spear = os.path.join(spear_topo_dir,'ocean_z.static.nc')
ds_grid_spear  = xarray.open_dataset(ftopo_spear)
LONS = ds_grid_spear['geolon'].data
LATS = ds_grid_spear['geolat'].data
HHS = ds_grid_spear.variables['depth_ocean'].data
HHS = -abs(HHS)
HHS = np.where(np.isnan(HHS), 0.5, HHS)

# Load mapping indices gmapi:
dirgmapi = config['filesystem']['spear_mom_gmapi']
flgmaph  = f'spear2mom_NEP_OB_gmapi_hpnt.nc'
flgmapu  = f'spear2mom_NEP_OB_gmapi_upnt.nc'
flgmapv  = f'spear2mom_NEP_OB_gmapi_vpnt.nc'
dflgmaph = os.path.join(dirgmapi, flgmaph)
dflgmapu = os.path.join(dirgmapi, flgmapu)
dflgmapv = os.path.join(dirgmapi, flgmapv)
# h-point indices
dsh = xarray.open_dataset(dflgmaph)


# MOM6 NEP topo/grid:
nsegm = 1
pthtopo     = gridfls['MOM6_NEP']['test']['pthgrid']
fgrid_mom   = gridfls['MOM6_NEP']['test']['fgrid']
ftopo_mom   = gridfls["MOM6_NEP"]["test"]["ftopo"]
dfgrid_mom  = os.path.join(pthtopo, fgrid_mom)
dftopo_mom  = os.path.join(pthtopo, ftopo_mom)
LONM, LATM  = mmom6.read_mom6grid(dfgrid_mom)
HHM         = mmom6.read_mom6depth(dftopo_mom)
jdm, idm    = np.shape(HHM)
hgrid       = xarray.open_dataset(os.path.join(pthtopo,fgrid_mom))
dsegm_nep   = mutob.segm_topo(nsegm, HHM, hgrid)
lon1        = dsegm_nep['lon_segm'].data[0]
lon2        = dsegm_nep['lon_segm'].data[-1]
lat1        = dsegm_nep['lat_segm'].data[0]
lat2        = dsegm_nep['lat_segm'].data[-1]

# OB segments:
for nsegm in range(1,5):
  print(f'nsegm={nsegm}')
  dset_segm_nep = mutob.segm_topo(nsegm, HHM, hgrid)
  lon_segm      = dset_segm_nep['lon_segm']
  lat_segm      = dset_segm_nep['lat_segm']
  npnts         = len(lon_segm)

  II = []
  JJ = []
  for ik in range(0,npnts,10):
    x0 = lon_segm[ik]
    y0 = lat_segm[ik]
    i0, j0 = mmisc.find_closest_indx([x0],[y0], LONS, LATS)
    II.append(i0)
    JJ.append(j0) 

  II = np.array(II).squeeze()
  JJ = np.array(JJ).squeeze()
  ni = len(II)
  darri = xarray.DataArray(II, dims=(f"npnts_{nsegm:03d}"),\
                           coords={f"npnts_{nsegm:03d}": np.arange(ni)})
  darrj = xarray.DataArray(JJ, dims=(f"npnts_{nsegm:03d}"),\
                           coords={f"npnts_{nsegm:03d}": np.arange(ni)})
 
  if nsegm == 1:
    dset_ij = xarray.Dataset({f"indx_segm{nsegm:03d}": darri, 
                              f"jndx_segm{nsegm:03d}": darrj})
  else:
    dstmp = xarray.Dataset({f"indx_segm{nsegm:03d}": darri, 
                            f"jndx_segm{nsegm:03d}": darrj})
    dset_ij = xarray.merge([dset_ij, dstmp])


# Find NEP OB segment on SPEAR grid:
#IIv,JJv = mmisc.find_closest_indx([lon1,lon2],[lat1,lat2], lon_spear, lat_spear)
#II, JJ  = mmisc.xsect_indx(IIv, JJv)
varnmi = f"indx_segm{nsegm:03d}"
varnmj = f"jndx_segm{nsegm:03d}"
INDX = dsh[varnmi].data
JNDX = dsh[varnmj].data
II   = INDX[:,0].squeeze()
JJ   = JNDX[:,0].squeeze()

# Keep only non-repeated grid points 
II, JJ    = mutob.discard_repeated_indx(II, JJ, keep_last=True)
dset_segm = mutob.segm_topo_spear(II, JJ, LONS, LATS, HHS, nsegm=nsegm)

# Det vertical grid and 2D section from SPEAR:
ZZ  = -ds.variables['z_interf'].data
dZ  = np.diff(ZZ)
_, ZM = mmom6.zz_zm_fromDZ(dZ)
dmm = A3d[:,JJ,II]
A2d =  mmom6.fill_land3d(dmm)

rmin, rmax = mutob.minmax_clrmap(A2d, cpnt=0)

if varnm == 'so':
  clrmp = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
elif varnm == 'thetao':
  clrmp = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[1,1,1])

#X        = dset_segm['dist_supergrid'].data
Xbtm     = dset_segm['dist_grid'].data
Hbtm     = dset_segm['topo_segm'].data
Hbtm     = np.where(Hbtm > 0, 0., Hbtm)
segm_nm  = dset_segm['segm_name'].data[0]

# Scale SPEAR distance to match NEP:
npnts = len(Xbtm)
Xbtm  =  np.linspace(0, 1., num=npnts, endpoint=True)*Xmax
X     = Xbtm.copy()

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
    rmax = 1.
  if varnm == 'so':
    rmin = 29.
    rmax = 35.
elif segm_nm == 'south':
  xl1 = 0
  if varnm == 'so':
    rmin = 34.0
    rmax = 34.9
else:
  xl1 = 0

xl2  = max(Xbtm)

sttl = f"SPEAR along NEP OB: {varnm} M={imonth} OB={segm_nm}"
mutob.plot_xsection(A2d, X, ZM, Hbtm, Xbtm, clrmp, rmin, rmax, xl1, xl2, sttl=sttl, stxt=stxt, fgnmb=1)

btx = 'check_spearOB.py'
bottom_text(btx)



cland = mpl.colormaps['PuBu_r']
cland.set_over(color=[0.3, 0.3, 0.3]) 
rmin = -5000
rmax = 0.
plt.ion()
fgnmb = 1
fig = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.25, 0.8, 0.7])
im  = ax1.pcolormesh(HHS, cmap = cland, 
                     vmin=rmin, vmax=rmax)
ax1.contour(LATS,[x for x in range(5,89,2)], linestyles='solid',
              linewidths=0.5, colors=[(0.7,0.7,0.7)])
ax1.contour(LATS,[70,80], linestyles='solid',
              linewidths=0.5, colors=[(1, 0, 0)])
ax1.contour(LONS,[x for x in range(0,359,5)], linestyles='solid',
              linewidths=0.5, colors=[(0.7,0.7,0.7)])
ax1.contour(LONS,[x for x in range(-359,-1,5)], linestyles='solid',
              linewidths=0.5, colors=[(0.7,0.7,0.7)])
ax1.contour(LONS,[-180,-90, 0, 90], linestyles='solid',
              linewidths=0.5, colors=[(1, 0, 0)])

for nsegm in range(1,5):
  II = dset_ij[f'indx_segm{nsegm:03d}'].data
  JJ = dset_ij[f'jndx_segm{nsegm:03d}'].data
  ax1.plot(II,JJ,'-', color=[0.0,0.7,0.9])

ax1.axis('scaled')
ax1.set_ylim([160, 320])

ax1.set_title('SPEAR grid & bathymetry')

btx = 'compare_spear_nep_bathymetry.py'
bottom_text(btx)




ax1.plot(II,JJ,'-', color=[0.9,0.4,0])

ax1.set_title('SPEAR grid & bathymetry and NEP OBs')


bottom_text(btx)
