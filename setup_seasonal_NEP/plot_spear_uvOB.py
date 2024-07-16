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
nsegm  = 2

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
lon1        = dsegm_nep['lon_segm'].data[0]
lon2        = dsegm_nep['lon_segm'].data[-1]
lat1        = dsegm_nep['lat_segm'].data[0]
lat2        = dsegm_nep['lat_segm'].data[-1]

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

# Det vertical grid and 2D section from SPEAR:
ZZ  = -ds.variables['z_interf'].data
dZ  = np.diff(ZZ)
_, ZM = mmom6.zz_zm_fromDZ(dZ)
dmm = U3dh[:,JJ,II]
U2d = mmom6.fill_land3d(dmm)
dmm = V3dh[:,JJ,II]
V2d = mmom6.fill_land3d(dmm)
S2d = np.sqrt(U2d**2 + V2d**2)

#rmin, rmax = mutob.minmax_clrmap(A2d, cpnt=0)
clrmp = mutil.colormap_temp(clr_ramp=[0.95,0.92,1])
clrmp.set_bad(color=[1,1,1])

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
if segm_nm == 'north':
  xl1 = 1900.
  rmin = 0.
  rmax = 0.1
elif segm_nm == 'east':
  xl1 = 5200
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

sttl = f"SPEAR along NEP OB: |U|  M={imonth} OB={segm_nm}"
mutob.plot_xsection(S2d, X, ZM, Hbtm, Xbtm, clrmp, rmin, rmax, xl1, xl2, sttl=sttl, stxt=stxt, fgnmb=1)

btx = 'plot_spear_uvOB.py'
bottom_text(btx)


f_plt = False
if f_plt:
  from matplotlib import cm
  cland = cm.get_cmap('PuBu_r',200)
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
  ax1.axis('scaled')
  ax1.plot(II,JJ,'-', color=[0.9,0.4,0])
 
  ax1.set_title('SPEAR grid & bathymetry and NEP OBs')


  bottom_text(btx)
