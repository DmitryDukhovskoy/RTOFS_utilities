"""
  Plot NEP domain bathymetry
"""
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from yaml import safe_load
import mod_mom6 as mmom6

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

from mod_utils_fig import bottom_text

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


#from matplotlib import cm
cland = mpl.colormaps['PuBu_r']
cland.set_over(color=[0.3, 0.3, 0.3]) 
rmin = -5000
rmax = 0.
plt.ion()
fgnmb = 1
fig = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.25, 0.8, 0.7])
im  = ax1.pcolormesh(HHM, cmap = cland, 
                     vmin=rmin, vmax=rmax)
ax1.axis('scaled')

ax1.set_title('NEP bathymetry')

btx = 'plot_NEPtopo.py'
bottom_text(btx)



