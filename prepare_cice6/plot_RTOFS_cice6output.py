"""
  Plot CICE6 output fields
  fomr a run with restart
  created from RTOFS CICE4
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray as xr 
import matplotlib.colors as colors
from yaml import safe_load
from mpl_toolkits.basemap import Basemap, cm
import argparse

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
import mod_read_hycom as mhycom
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

def grid_rad2dgr(ulat,ulon, f180 = True):
  """ 
  Convert CICE grid coordinates from radians to degrees.
  If f180 is True, convert longitude to the range
  -180 <= lon < 180.
  """
  rdn2dgr = 180.0 / np.pi
  ulat = ulat * rdn2dgr
  ulon = ulon * rdn2dgr

  # Normalize longitude to [0, 360)
  ulon = np.mod(ulon, 360.0)

  # Optionally convert to [-180, 180]
  if f180:
      ulon = np.where(ulon > 180.0, ulon - 360.0, ulon)

  ulat = np.clip(ulat, -89.99999, 89.99999)

  return ulat, ulon

def check_depth_file(pthdpth, dpthfl):
  """
  Checks if input topography file is netcdf or unformatted binary *.a.
  """
  fldptha  = fldpthb = None
  topo_nc = topo_ab = False

  if dpthfl.endswith('.nc'):
      fldpthnc = pthdpth
      fdpthin = os.path.join(pthdpth, dpthfl)
      topo_nc = True
  elif dpthfl.endswith('.a'):
      fldptha = dpthfl
      fldpthb = fldptha.replace('.a', '.b')
      ftopo   = fldptha.removesuffix('.a')
      topo_ab = True
  else:
      raise ValueError(f"topo file {dpthfl} not recognized, expected *.a or *.nc")

  return fldptha, fldpthb, topo_nc, topo_ab


flcice = 'iceh_01h.2025-12-15-03600.nc'
pthcice = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS_output/cice6'
fyaml    = 'restart_cice6.yaml'


parser = argparse.ArgumentParser()
parser.add_argument("--flcice", help=f"CICE6 output file name, default={flcice}", default=flcice)
parser.add_argument(
    "--varnm",
    help="variable to plot",
    choices=['aice','ithkn'],
    required=True
)
parser.add_argument(
           "--pthcice", 
           help=f"Directory of the output file, default={pthcice}", 
           default=pthcice
)
parser.add_argument("--regn", help=f"where icon incerted: south, north, global", default='north')
args = parser.parse_args()

flcice  = args.flcice
pthcice = args.pthcice 
varnm   = args.varnm 
regn    = args.regn 

with open(fyaml) as ff:
  PATHS = safe_load(ff)

# Get RTOFS CICE6 grid
pthgrd  = PATHS["grid_topo"]["cice6"]["pthgrid"]
grdfl   = PATHS["grid_topo"]["cice6"]["filegrid"]
fgrdin  = os.path.join(pthgrd, grdfl)
pthtopo = PATHS["grid_topo"]["cice6"]["pthtopo"]
fldepth = PATHS["grid_topo"]["cice6"]["filedepth"]

# Read CICE6 
with xr.open_dataset(fgrdin) as dset:
  hlat_rad = dset['ulat'].data.squeeze()
  hlon_rad = dset['ulon'].data.squeeze()

hlat, hlon = grid_rad2dgr(hlat_rad, hlon_rad, f180 = True)

# Normalize hlon to avoid bougs longitude values for the projections:
# The near-pole grid points may result in infinity for stereographic projection
hlon2 = ((hlon + 180) % 360) - 180

JDIM, IDIM = hlon.shape
JDIM = JDIM + 1   # ocean grid has + 1 row

# Check type of topo file:
fldptha, fldpthb, topo_nc, topo_ab = check_depth_file(pthtopo, fldepth)
if topo_ab:
  ftopo = fldptha.removesuffix('.a')
elif topo_nc:
  ftopo = fldepth

# Read RTOFS topo:
# Note that RTOFS grid has +1 row at the top compared to CICE6  <--- old
if topo_ab:
  HH = mhycom.read_topo(pthtopo, ftopo, IDIM, JDIM)
elif topo_nc:
  with xr.open_dataset(os.path.join(pthtopo, ftopo)) as dset:
    HH = dset["depth"].values

# Convert to negative depths, land assigned arbitrary positive value
if np.nanmin(HH) >= 0.:
  LMsk = (HH < 1e-10) | np.isnan(HH)
  HH = -HH
  HH[LMsk] = 100

assert HH.shape == hlon.shape, f"Topo shape {HH.shape} and lon shape mismatch {hlon.shape}"
#HH = HH[:-1,:]     # discard the extra row, obsolete


# Read output
dflcice = os.path.join(pthcice, flcice)
print(f"Reading cice6 output: {dflcice}")

if varnm == 'aice' or varnm == 'iconc':
  varnc = 'aice_h'
elif varnm == 'ithkn' or varnm == 'hice':
  varnc = 'hi_h'

with xr.open_dataset(dflcice) as ds_out:
  aice = ds_out['aice_h'].values.squeeze()  # partial area by cats
  A2d  = ds_out[varnc].values.squeeze()


match varnm:
  case 'aice' | 'icon':
    clrmp = mclrmps.colormap_conc_blue()
    rmin = 0.
    rmax = 1.
  case 'ithkn':
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 3.
  case _:
    raise ValueError(f"Unknown variable name: '{varnc}'")

clrmp.set_bad(color=[0.2, 0.2, 0.2])

sttl = f"{varnm} {varnc} \nCICE6: {dflcice}"

# Subset fields to avoid inifinity in projected coordinates:
jsub = 1650

if regn == 'south':
  m = Basemap(projection='spstere', boundinglat=-55, lon_0=180, resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
  Asub = A2d[:jsub,:]
  aice_sub = aice[:jsub,:]
  lons = hlon[:jsub,:]
  lats = hlat[:jsub,:]
  HHs  = HH[:jsub,:]

elif regn == 'north':
  m = Basemap(projection='npstere', boundinglat=60, lon_0=-10, resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)
  Asub = A2d[jsub:,:]
  aice_sub = aice[jsub:,:]
  lons = hlon[jsub:,:]
  lats = hlat[jsub:,:]
  HHs  = HH[jsub:,:]

xh, yh = m(lons, lats)
dx = np.diff(xh, axis=1)
if np.max(abs(dx)) > 1e12:
  print(f"WARN: Basemap produced pathological projected coordinates max(dx) = {np.max(abs(dx))} ...")
 

plt.ion()

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
img1 = ax1.pcolormesh(xh, yh, Asub, cmap=clrmp, vmin=rmin, vmax=rmax)
# Overlay land
#ax1.pcolormesh(xh, yh, land_overlay, cmap=land_cmap, zorder=5)

ax1.contour(xh, yh, aice_sub, [0.15], linestyles='solid', colors=[(0.4,0.4,0.4)], linewidths=1)
ax1.contour(xh, yh, HHs, [0.], linestyles='solid', colors=[(0.,0.,0.)], linewidths=1)

ax1.set_title(sttl)


# Colorbars
ax3 = fig1.add_axes([0.1, 0.05, 0.8, 0.02])
clb = plt.colorbar(img1, cax=ax3, orientation='horizontal', extend='both')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_RTOFS_cice6output.py'
bottom_text(btx, pos=[0.1,0.01])



