"""
  Plot some CICE6 restart fields
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
from copy import copy
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
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_sis2_relax as msisrlx
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

rest_date = 20250103
rest_hr   = 0
hunits    = 'cm'
regn = 'south'
flrst_in = 'cice_model.res.20250103.00.iconc_thkn.snow.nc'

yrR = mmR = ddR = hrR = None

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
    choices=['aicen','vicen','vsnon','Tsfcn','alvl','vlvl','apnd','hpnd','ipnd','dhs','iceumask'],
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
yrR, mmR, ddR, hrR, mintR = mc6util.get_date_filename(flrst_in)
rest_date = int(yrR*1e4 + mmR*100 + ddR)
rest_hr = hrR

syst_info = os.uname()
machine = syst_info.nodename

if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
elif 'an' in machine:
  print("Running on PPAN node:", machine)
  node_nm = "ppan"
else:
  print("Unknown machine:", machine)

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

pthrest = os.path.join(pths_ufs[node_nm]["MOM6"]["pthrest"],'new')
#pthrest = '/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_da'
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]

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
ds_in = xarray.open_dataset(dflrst_in)
ds_out = ds_in.copy(deep=True)
ds_in.close()

with xarray.open_dataset(dflrst_in) as ds_out:
  aicen = ds_out['aicen'].data  # partial area by cats
  A3d    = ds_out[varnc].data 
aice = np.nansum(aicen, axis=0).squeeze()  # aggreagetd ice conc

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


# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
    
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

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

sttl = f"{varnc} cat={icat}:  CICE6 restart {flrst_in}"

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
#lons, lats = m.makegrid(idim, jdim) # get lat/lons of ny by nx evenly spaced grid.
parallels = np.arange(-80,-10,10.)
meridians = np.arange(-360,359.,45.)
xl1 = -8.5e6
xl2 = -0.9e6
yl1 = xl1
yl2 = xl2

xh, yh = m(hlon,hlat)

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
img1 = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
# Overlay land
#ax1.pcolormesh(xh, yh, land_overlay, cmap=land_cmap, zorder=5)

ax1.contour(xh, yh, aice, [0.15], linestyles='solid', colors=[(0.4,0.4,0.4)], linewidths=1)
ax1.contour(xh, yh, HH, [0.], linestyles='solid', colors=[(0.,0.,0.)], linewidths=1)

ax1.set_title(sttl)
ax1.set_xlim([xl1, xl2])
ax1.set_ylim([yl1, yl2])
ax1.invert_yaxis()
ax1.invert_xaxis()


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



