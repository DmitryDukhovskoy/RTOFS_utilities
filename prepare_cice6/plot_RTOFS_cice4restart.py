"""
  Plot CICE4 restart fields
  from RTOFS CICE4
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray as xr 
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
import mod_read_hycom as mhycom
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

class CICE: 
  def __init__(self, nx, ny, ncat, nilyr, nslyr):
    self.ncat = ncat
    self.nilyr = nilyr            
    self.nslyr = nslyr
    self.nx = nx
    self.ny = ny

    self.ntilyr = self.ncat * self.nilyr
    self.ntslyr = self.ncat * self.nslyr

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

def read_rest_cice4(fid, nx, ny):
    """               
    Read a 2D field from an open CICE4 restart file.
    CICE4 restart files are unformatted sequential binary records
    in big-endian format. 
    """         
    recS = np.fromfile(fid, dtype='>i4', count=1)[0]
    A     = np.fromfile(fid, dtype='>f8', count=nx*ny)
    A     = np.reshape(A,(ny,nx), order='C')
    recE = np.fromfile(fid, dtype='>i4', count=1)[0]
    if recS != recE:
       raise ValueError(f"Record length mismatch: {recS} != {recE}")
    
    return A

def print_minmax(sfld, A):
    """Print min/max statistics of a numpy array."""
    print(f'   {sfld} min/max:  {np.nanmin(A)} / {np.nanmax(A)}')
    return

def read_cice4_layers(fid, nlrs, nx, ny, label):
    """ Read CICE4 fields by layers. """
    print(f'\n Reading {label}:')
    fld = np.zeros((nlrs, ny, nx), dtype=np.float64)

    for k in range(nlrs):
        A = read_rest_cice4(fid, nx, ny)
        fld[k, :, :] = A
        print_minmax(f"{k+1} {label}", A)

    return fld

def read_cice4_2D(fid, nx, ny, varnm):
    """Read CICE4 2D fields."""
    print(f'\nReading {varnm}:')
    A = read_rest_cice4(fid, nx, ny)
    print_minmax(varnm, A)

    return A


rest_date = 20250103
rest_hr   = 0
hunits    = 'cm'
regn = 'south'
#flrst_in = 'rtofs_glo.20001216_00000.restart_cice.nc'
flrst_in = 'rtofs_glo.t00z.n-24.restart_cice'
fyaml    = 'restart_cice6.yaml'
spval     = 1.e30             

# alvl the fraction of the level ice area
# vlvl the volume of the level ice area
# apnd -  the fraction of ponds of ice area, 
#         i.e. apnd = 0.9 and aice=0.25 --> pond fraction of a grid cell = apnd*aice[*alvl] - depending on 
#         pond parameterization, alvl - area of level ice
# hpnd - and depth of the ponds in a cell
# iceumask - The mask where there is ice in a gridcell, dynamics ice extent mask (U-cell)
# dhs - local diff. in snow depth on sea ice and pond ice (see: icepack_meltpond_lvl.F90)
# ipnd - melt pond refrozen lid thickness (see: icepack_therm_vertical.F90)
parser = argparse.ArgumentParser()
parser.add_argument("--flrst_in", help=f"rest file name, default={flrst_in}", type=str)
parser.add_argument(
    "--varnm",
    help="variable to plot",
    choices=['aicen','vicen','vsnon','Tsfcn','alvl','vlvl','apnd','hpnd','ipnd','iceumask'],
    type=str,
    required=True
)
parser.add_argument("--cat", help="category to plot: 1, ..., ncat, =0 - aggregated", type=int, required=True)
parser.add_argument("--regn", help=f"where icon incerted: south, north, global, default={regn}", type=str)
args = parser.parse_args()

flrst_in  = args.flrst_in if args.flrst_in else flrst_in
varnc = args.varnm if args.varnm else None
icat  = args.cat if args.cat is not None else -999
regn  = args.regn if args.regn else regn

# Derive dates assuming file nameing is cice_restart.res.YYYYMMDD.XX[XXX]
#yrR, mmR, ddR, hrR, mintR = mc6util.get_date_filename(flrst_in)
#rest_date = int(yrR*1e4 + mmR*100 + ddR)
#rest_hr = hrR

with open(fyaml) as ff:
  PATHS = safe_load(ff)

cicerst4 = PATHS["rest_names"]["cice4"]["flnm"]
cicerstT = PATHS["rest_names"]["tmplt"]["flnm"]  # template CICE6
pthrst4  = PATHS["cice_paths"]["cice4"]["pth"]
pthrstT  = PATHS["cice_paths"]["tmplt"]["pth"]

fl_restart4 = os.path.join(pthrst4, cicerst4)
fl_restartT = os.path.join(pthrstT, cicerstT)

if flrst_in is None:
  flrst_in = PATHS["rest_names"]["cice4"]["flnm"]

pthrest = pthrst4
#pthrest  = PATHS["cice_paths"]["cice6"]["pth"]
#pthrest = os.path.join(pths_ufs[node_nm]["MOM6"]["pthrest"],'new')
#pthrest = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/cice_restarts/cice6'
#pthrest = '/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_da'


# CICE parameters:
puny      = 1.e-11
c0        = 0.0
c1        = 1.0
c2        = 2.0
p5        = 0.5
Lsub      = 2.835e6    # latent heat sublimation fw (J/kg)
Lvap      = 2.501e6    # latent heat vaporization fw (J/kg)
Lfresh    = Lsub - Lvap # latent heat of melting of fresh ice (J/kg)
cp_ice    = 2106.       # specific heat of fresh ice (J/ kg/K)
rhos      = 330.        # density of snow (kg/m3)
hs_min    = 1.e-4       # min snow thickness for computing Tsno (m)
nsal      = 0.407
msal      = 0.573
min_salin = 0.1      # threshold for brine pocket treatment
saltmax   = 3.2        # max S at ice base
hg        = 1.e20    # bad values, land mask, etc.
nslyr     = 1   # snow layers
Tmin      = -100.   # minimum snow T


dflrst_in = os.path.join(pthrest, flrst_in)
print(f"Reading restart: {dflrst_in}")
# Grid CICE4 unformatted binary file.
pthgrd4 = PATHS["grid_topo"]["cice4"]["pthgrid"]
grdfl4  = PATHS["grid_topo"]["cice4"]["filegrid"]
fgrdin4 = os.path.join(pthgrd4, grdfl4)

# Create object with CICE4 grid parameters.
nx    = PATHS["cice_params"]["cice4"]["nx"]
ny    = PATHS["cice_params"]["cice4"]["ny"]
ny = 3298
ncat  = PATHS["cice_params"]["cice4"]["ncat"]
nilyr = PATHS["cice_params"]["cice4"]["nilyr"]
nslyr = PATHS["cice_params"]["cice4"]["nslyr"]
cice4 = CICE(nx, ny, ncat, nilyr, nslyr)



# Read fields from the CICE4 restart file.
if not os.path.exists(fl_restart4):
  raise FileNotFoundError(f"Does not exist: {fl_restart4}")
    
print(f'Reading restart: {fl_restart4}')
with open(fl_restart4, 'rb') as fid:
  # Read the 1st sequential record of CICE4 restart file.
  # recS:     4-byte record-length marker (start marker)
  # istep:    current model step
  # runtime:  total elapsed model time (s)
  # frtime:   elapsed time since the last forcing update (s)
  # recE:     record-length marker (end marker)
  recS    = np.fromfile(fid, dtype='>i4', count=1)[0]
  istep   = np.fromfile(fid, dtype='>i4', count=1)[0]
  runtime = np.fromfile(fid, dtype='>f8', count=1)[0]
  frtime  = np.fromfile(fid, dtype='>f8', count=1)[0]
  recE    = np.fromfile(fid, dtype='>i4', count=1)[0]

  if recS != recE:
    raise ValueError(f"Record length mismatch: {recS} != {recE}")

  print(
  f"Restart: step={istep}, "
  f"total time={runtime/(3600*24*365.25):.2f} yr, "
  f"forcing update={frtime/3600:.1f} hr ago"
  ) 

  nx     = cice4.nx
  ny     = cice4.ny
  ncat   = cice4.ncat
  ntilyr = cice4.ntilyr  # total # of icelrs * cat 
  ntslyr = cice4.ntslyr

  aicen = np.zeros((ncat,ny,nx), dtype='float64')
  vicen = np.zeros((ncat,ny,nx), dtype='float64')
  vsnon = np.zeros((ncat,ny,nx), dtype='float64')
  trcrn = np.zeros((ncat,ny,nx), dtype='float64')

  for n in range(ncat):
      print(f" Category {n+1}")
      for varname, arr in [
          ("ice area", aicen),
          ("ice vol",  vicen),
          ("snow vol", vsnon),
          ("surf T",   trcrn),
      ]:
          A = read_rest_cice4(fid, nx, ny)
          arr[n, :, :] = A
          print_minmax(varname, A)

  eicen = read_cice4_layers(fid, ntilyr, nx, ny, "eicen")
  esnon = read_cice4_layers(fid, ntslyr, nx, ny, "esnon")
  uvel  = read_cice4_2D(fid, nx, ny, 'uvel')
  vvel  = read_cice4_2D(fid, nx, ny, 'vvel')
  uvelE = None
  vvelN = None
  fsnow = None
  iage  = None
  alvl  = None
  vlvl  = None
  apnd  = None
  hpnd  = None
  ipnd  = None
  dhs   = None
  ffrac = None
  coszen = None

  scale_factor = read_cice4_2D(fid, nx, ny, 'scale factor')
  swvdr        = read_cice4_2D(fid, nx, ny, 'sh/wave vis direct')
  swvdf        = read_cice4_2D(fid, nx, ny, 'sh/wave vis diff')
  swidr        = read_cice4_2D(fid, nx, ny, 'sh/wave IR dir')
  swidf        = read_cice4_2D(fid, nx, ny, 'sh/wave IR diff')
  strocnxT     = read_cice4_2D(fid, nx, ny, 'ocean stress x-comp')
  strocnyT     = read_cice4_2D(fid, nx, ny, 'ocean stress y-comp')

  stressp = {}
  for fld in ["stressp_1", "stressp_3", "stressp_2", "stressp_4"]:
      stressp[fld] = read_cice4_2D(fid, nx, ny, fld)
  stressm = {}
  for fld in ["stressm_1", "stressm_3", "stressm_2", "stressm_4"]:
      stressm[fld] = read_cice4_2D(fid, nx, ny, fld)
  stress12 = {}

  for fld in ["stress12_1", "stress12_3", "stress12_2", "stress12_4"]:
      stress12[fld] = read_cice4_2D(fid, nx, ny, fld)

  iceumask = read_cice4_2D(fid, nx, ny, 'ice umask')             
  sst      = read_cice4_2D(fid, nx, ny, 'ocean mixed layer sst')
  frzmlt   = read_cice4_2D(fid, nx, ny, 'frzmlt')


# Mask out land points:
print(' Masking out fields ')
maskval = 0.5 * spval
for A in [
    aicen, vicen, vsnon, trcrn, eicen, esnon,
    uvel, vvel, scale_factor, swvdr, swvdf, swidr, swidf,
    strocnxT, strocnyT, sst, frzmlt
]:
    A[A > maskval] = np.nan

for A in stressp.values():
    A[A > maskval] = np.nan

for A in stressm.values():
    A[A > maskval] = np.nan

for A in stress12.values():
    A[A > maskval] = np.nan


# Collect fields with names and corresponding arrays
cice4_vars = {
'uvel': uvel,
'vvel': vvel,
'scale_factor': scale_factor,
'swvdr': swvdr,
'swvdf': swvdf,
'swidr': swidr,
'swidf': swidf,
'strocnxT': strocnxT,
'strocnyT': strocnyT,
'iceumask': iceumask,
'fsnow': fsnow,
'aicen': aicen,
'vicen': vicen,
'vsnon': vsnon,
'iage': iage,
'alvl': alvl,
'vlvl': vlvl,
'apnd': apnd,
'hpnd': hpnd,
'ipnd': ipnd,
'dhs': dhs,
'ffrac': ffrac,
'Tsfcn': trcrn,
'coszen': coszen,
}

# Add all stress fields
cice4_vars.update(stressp)
cice4_vars.update(stressm)
cice4_vars.update(stress12)


A3d = cice4_vars[varnc]
if A3d.ndim == 3:
  ncat, jdim, idim = A3d.shape
  if icat == 0:
    A2d = np.nansum(A3d, axis=0).squeeze()
  else:
    A2d = A3d[icat-1,:].squeeze()
elif A3d.ndim == 2:
  A2d = A3d.copy()
  jdim, idim = A3d.shape
  ncat = 1

# Get aice:
A3d = cice4_vars['aicen']  
aice = np.nansum(A3d, axis=0).squeeze()  # aggreagetd ice conc


# Read  CICE4 grid coordinates, Bu points.
hlat_rad = mc6util.read_cice4_grid(fgrdin4, 'ulati', IDM=cice4.nx, JDM=cice4.ny)
hlon_rad = mc6util.read_cice4_grid(fgrdin4, 'uloni', IDM=cice4.nx, JDM=cice4.ny)


# Get RTOFS CICE6 grid
pthgrd  = PATHS["grid_topo"]["cice6"]["pthgrid"]
grdfl   = PATHS["grid_topo"]["cice6"]["filegrid"]
fgrdin  = os.path.join(pthgrd, grdfl)
pthtopo = PATHS["grid_topo"]["cice6"]["pthtopo"]
fldepth = PATHS["grid_topo"]["cice6"]["filedepth"]

hlat, hlon = grid_rad2dgr(hlat_rad, hlon_rad, f180 = True)

JDIM, IDIM = hlon.shape
JDIM = JDIM + 1   # ocean grid has + 1 row

# Check type of topo file:
fldptha, fldpthb, topo_nc, topo_ab = check_depth_file(pthtopo, fldepth)
if topo_ab:
  ftopo = fldptha.removesuffix('.a')
elif topo_nc:
  ftopo = topo_nc

# Read RTOFS topo:
# Note that RTOFS grid has +1 row at the top compared to CICE6
HH = mhycom.read_topo(pthtopo, ftopo, IDIM, JDIM)
HH = HH[:-1,:]     # discard the extra row

match varnc:
  case 'aicen' | 'alvl' | 'vlvl' | 'apnd' | 'hpnd' | 'ipnd' | 'dhs' | 'iceumask':
    clrmp = mclrmps.colormap_conc()
    rmin = 0.
    rmax = 1.
  case 'vicen':
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 3.
  case 'vsnon':
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 1.
  case 'Tsfcn':
    clrmp = mclrmps.colormap_cold_warm()
    rmin = -3.
    rmax = 0.
    clrmp.set_over([1.,0.,1])
  case _:
    raise ValueError(f"Unknown variable name: '{varnc}'")

clrmp.set_bad(color=[0.2, 0.2, 0.2])

sttl = f"{varnc} cat={icat}:  CICE4 restart {flrst_in}"

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

btx = 'plot_cicerestart.py'
bottom_text(btx, pos=[0.1,0.01])



