"""
  Plot UV from SPEAR monthly clim along NEP OB on SPEAR grid
  2D vertical sections
  use subset domains
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

imonth = 1
nsegm  = 4   # segment # 1,2,3,4
klr    = 1   # model layer to plot

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

# SPEAR topo, subset to NEP domain similar to data files:
ftopo_spear = os.path.join(spear_topo_dir,'ocean_z.static.nc')
spear_topo = xarray.open_dataset(ftopo_spear)
lonW = config['domain']['west_lon']
lonE = config['domain']['east_lon']
LONh = xarray.open_dataset(ftopo_spear).variables['geolon'].data

if lonW > np.max(LONh):
  lonW = lonW-360.
elif lonW < np.min(LONh):
  lonW = lonW+360.

if lonE > np.max(LONh):
  lonE = lonE-360.
elif lonE < np.min(LONh):
  lonE = lonE+360.

lon_slice = slice(lonW, lonE)
lat_slice = slice(config['domain']['south_lat'], config['domain']['north_lat'])
ds_grid_spear  = xarray.open_dataset(ftopo_spear).sel(xh=lon_slice, xq=lon_slice, \
                 yh=lat_slice, yq=lat_slice)
HHS = ds_grid_spear.variables['depth_ocean'].data
HHS = -abs(HHS)
HHS = np.where(np.isnan(HHS), 0.5, HHS)

# Load SPEAR data
ds = mutob.load_var(spear_dir, 'uo', fzint=True)
U3d = ds.variables['uo'].data[imonth-1,:,:,:].squeeze()

ds = mutob.load_var(spear_dir, 'vo', fzint=True)
V3d = ds.variables['vo'].data[imonth-1,:,:,:].squeeze()

# V-points - add 2 jslices to have same dims as H after collocation
V3d = np.append(V3d, np.expand_dims(V3d[:,-1,:], axis=1), axis=1)
V3d = np.append(V3d, np.expand_dims(V3d[:,-1,:], axis=1), axis=1)
#LONs, LATs = mutob.subset_spear_coord('u')  # coordinates of the u-points
LONs, LATs = mutob.subset_spear_coord('h')  # coordinates of the u-points
# Make Longitudes same range:
LONs = np.where(LONs<0., LONs+360., LONs)


# Collocate U V components at h-points:
kdim, jdim, idim = U3d.shape
U3dh = np.zeros((kdim, jdim, idim-1))
kdimv, jdimv, idimv = V3d.shape
V3dh = np.zeros((kdimv, jdimv-1, idimv))
for kk in range(kdim):
  dmm = mmom6.collocateU2H(U3d[kk,:,:].squeeze(), 'symmetr', f_land0 = False)
  U3dh[kk,:,:] = dmm

for kk in range(kdim):
  dmm = mmom6.collocateV2H(V3d[kk,:,:].squeeze(), 'symmetr', f_land0 = False)
  V3dh[kk,:,:] = dmm


# Load mapping indices gmapi:
dirgmapi = config['filesystem']['spear_mom_gmapi']
flgmaph  = f'spear2mom_NEP_OB_gmapi_hpnt.nc'
flgmapu  = f'spear2mom_NEP_OB_gmapi_upnt.nc'
flgmapv  = f'spear2mom_NEP_OB_gmapi_vpnt.nc'
dflgmaph = os.path.join(dirgmapi, flgmaph)
dflgmapu = os.path.join(dirgmapi, flgmapu)
dflgmapv = os.path.join(dirgmapi, flgmapv)
# h-point indices because U V will be collocated into h-points
dsh = xarray.open_dataset(dflgmaph)

# MOM6 NEP topo/grid:
pthtopo     = gridfls['MOM6_NEP']['test']['pthgrid']
fgrid_mom   = gridfls['MOM6_NEP']['test']['fgrid']
ftopo_mom   = gridfls["MOM6_NEP"]["test"]["ftopo"]
dfgrid_mom  = os.path.join(pthtopo, fgrid_mom)
dftopo_mom  = os.path.join(pthtopo, ftopo_mom)
#LONM, LATM  = mmom6.read_mom6grid(dfgrid_mom)
HHM         = mmom6.read_mom6depth(dftopo_mom)
jdm, idm    = np.shape(HHM)
hgrid       = xarray.open_dataset(os.path.join(pthtopo,fgrid_mom))
dsegm_nep   = mutob.segm_topo(nsegm, HHM, hgrid)

# Scale SPEAR distances to NEP distance for comparison
dset_segm_nep = mutob.segm_topo(nsegm, HHM, hgrid)
Xnep          = dset_segm_nep['dist_supergrid'].data
Xmax          = np.max(Xnep)

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
dset_segm = mutob.segm_topo_spear(II, JJ, LONs, LATs, HHS, nsegm=nsegm)
lon_spear = dset_segm['lon_segm'].data
lat_spear = dset_segm['lat_segm'].data

# 2D section from SPEAR:
#dmm = U3dh[:,JJ,II]
#U2d = mmom6.fill_land3d(dmm)
#dmm = V3dh[:,JJ,II]
#V2d = mmom6.fill_land3d(dmm)


lon1      = dsegm_nep['lon_segm'].data[0]
lon2      = dsegm_nep['lon_segm'].data[-1]
lat1      = dsegm_nep['lat_segm'].data[0]
lat2      = dsegm_nep['lat_segm'].data[-1]
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

lon_mom = dsegm_nep['lon_segm'].data
lat_mom = dsegm_nep['lat_segm'].data
if segm_nm == 'north':
  uscl = 250
  xl1 = -5
  xl2 = 20
  yl1 = 86
  yl2 = 118
elif segm_nm == 'east':
  uscl = 250
  xl1 = 5
  xl2 = 70
  yl1 = 105
  yl2 = 125
elif segm_nm == 'south':
  uscl = 100
  xl1 = 50
  xl2 = 110
  yl1 = -5
  yl2 = 30
elif segm_nm == 'west':
  uscl = 100
  xl1 = -10
  xl2 = 80
  yl1 = 0
  yl2 = 60

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
im  = ax1.pcolormesh(HHS, cmap = cland,
                     vmin=rmin, vmax=rmax)
ax1.contour(LONs, [x for x in range(150,290,5)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])
ax1.contour(LATs, [x for x in range(5, 89, 5)], linestyles='solid',
                linewidths=0.5, colors=[(0.7,0.7,0.7)])

ax1.axis('scaled')
for ksgm in range(1,5):
  varnmi = f"indx_segm{ksgm:03d}"
  varnmj = f"jndx_segm{ksgm:03d}"
  INDX = dsh[varnmi].data
  JNDX = dsh[varnmj].data
  II   = INDX[:,0].squeeze()
  JJ   = JNDX[:,0].squeeze()
  II, JJ    = mutob.discard_repeated_indx(II, JJ, keep_last=True)
  ax1.plot(II,JJ,'-', color=[0.9,0.4,0])

# Plot vectors
dstep = 20
varnmi = f"indx_segm{nsegm:03d}"
varnmj = f"jndx_segm{nsegm:03d}"
INDX   = dsh[varnmi].data
JNDX   = dsh[varnmj].data
lindx  = len(lon_mom) 
II     = INDX[:,0].squeeze()
JJ     = JNDX[:,0].squeeze()
for kk in range(1,lindx, dstep):
  k0 = klr-1
  i0 = II[kk]
  j0 = JJ[kk]
  hh = HHS[j0, i0]
  if hh >= 0.: continue
  uu = U3dh[k0, j0, i0]*uscl  # cm/ s
  vv = V3dh[k0, j0, i0]*uscl  # cm/s
  ax1.quiver(i0, j0, uu, vv, angles='xy', scale_units='xy', \
             scale=1, width=0.002, headlength=4, color=[0.,1.,0.8])

ax1.set_xlim([xl1, xl2])
ax1.set_ylim([yl1, yl2])
ax1.set_title(f'SPEAR U vectors, M={imonth} layer = {klr} at NEP OB={nsegm} {segm_nm}')
ax1.set_xlabel('SPEAR I grid points')
ax1.set_ylabel('SPEAR J grid points')

btx = 'plot_spear_uvOB_vectors.py'
bottom_text(btx)


