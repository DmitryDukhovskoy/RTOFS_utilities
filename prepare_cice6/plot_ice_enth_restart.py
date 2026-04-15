"""
  Plot ice enthalpy  CICE6 restart fields
  q(i,n,k) - ice layer enthalpy (J/m3)
  e(i,n,k) - ice layer energy (J/m2)
  e(i,n,k) = v(i,n) / Nice * q(i,n,k), v(i,n) - ice volume, n - cat, Nice - N ice layers
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
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
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

rest_hr   = 0
flrst_in = '20240701.00.cice_model.res.nc'

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
parser.add_argument("--flrst", help=f"rest file name, default={flrst_in}", type=str)
parser.add_argument("--ilr", help="Ice layer to plot, if None - aggregated", 
                    choices=[1,2,3,4,5,6,7], type=int)
parser.add_argument("--fld", help="Plot energy J/m2_ice, or J/m2_cell, or enthalpy J/m3",
                    choices=['eJm2cell', 'eJm2ice', 'qJm3'], required=True, type=str)
parser.add_argument("--regn", help=f"region: south, north", required=True, type=str)
args = parser.parse_args()

flrst_in  = args.flrst if args.flrst else flrst_in
regn  = args.regn if args.regn else regn
ilr   = args.ilr if args.ilr else None
plt_fld = args.fld

# Derive dates assuming file nameing is cice_restart.res.YYYYMMDD.XX[XXX]
yrR, mmR, ddR, hrR, mintR = mc6util.get_date_filename(flrst_in)

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

#pthrest = '/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_da'
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthrest = '/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_sfs_C192mx025/ice/cice6_global'

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

with xarray.open_dataset(dflrst_in) as ds:
  aicen = ds['aicen'].data  # partial area by cats
  vicen = ds['vicen'].data  # ice volume m3/m2 cell
  qice = [
      ds['qice001'].data,
      ds['qice002'].data,
      ds['qice003'].data,
      ds['qice004'].data,
      ds['qice005'].data,
      ds['qice006'].data,
      ds['qice007'].data,
  ]

ncat = aicen.shape[0]
Nilrs = len(qice)

# Initialize total ice energy (J/m2)
Eice = np.zeros_like(aicen[0])   # J/m2 (grid cell or ice area)
Qnum = np.zeros_like(aicen[0])   # numerator for volume-weighted enthalpy
Vsum = np.zeros_like(aicen[0])   # denominator (volume sum)

# Total ice concentration
aice = np.sum(aicen, axis=0)

# Select layers
if ilr is None:
  layers = range(Nilrs)
else:
  if ilr < 1 or ilr > Nilrs:
    raise ValueError(f"ilr must be between 1 and {Nilrs}")
  layers = [ilr - 1]

# Loop over layers and categories
j0 = 956
i0 = 476
for il in layers:
  q_layer = qice[il]   # (ncat, nj, ni)

  if j0 is not None and i0 is not None:
    print(f"j0={j0} i0={i0}, Layer {il+1}: ")

  # Debug: sum over cats:
  V_lr = np.zeros_like(aicen[0])
  E_lr = np.zeros_like(aicen[0])
  for icat in range(ncat):
    # volume per layer (m3/m2)
    v_layer = vicen[icat] / Nilrs

    # energy (J/m2_cell)
    E_layer = v_layer * q_layer[icat]

    #if plt_fld in ['eJm2cell', 'eJm2ice']:
    Eice += E_layer

    #elif plt_fld == 'qJm3':
    # volume-weighted enthalpy
    # sum over all cats and ice layers:
    Qnum += E_layer
    Vsum += v_layer

    # Sum over cats in 1 lr:
    V_lr += v_layer
    E_lr += E_layer
    
    if j0 is not None and i0 is not None:
      print(
          f"  cat {icat+1}: "
          f"q_lr={q_layer[icat,j0,i0]:.3e}, "
          f"v_lr={v_layer[j0,i0]:.3e}, "
          f"E_lr={E_layer[j0,i0]:.3e}"
      )

  if j0 is not None and i0 is not None:
    Q_lr = np.sum(q_layer[:,j0,i0])
    vice_lr = V_lr[j0,i0]
    Eice_lr = E_lr[j0,i0]
    print(
        f"  Sum over cat Layer: Q_lr={Q_lr:.3e}, "
        f"V_lr={vice_lr:.3e}, "
        f"E_lr={Eice_lr:.3e}"
    )


# Final fields
cff = 1e-8

Eice_ice = np.divide(Eice, aice,
                 out=np.zeros_like(aice),
                 where=aice != 0)

Qice = np.divide(Qnum, Vsum,
                 out=np.zeros_like(Qnum),
                 where=Vsum != 0)

if plt_fld == 'eJm2cell':
  A2d = Eice * cff
  stxt = 'Eice (J/m2_cell)'

elif plt_fld == 'eJm2ice':
  A2d = Eice_ice * cff
  stxt = 'Eice (J/m2_ice)'

elif plt_fld == 'qJm3':
  A2d = Qice * cff
  stxt = 'Qice (J/m3)'

if j0 is not None and i0 is not None:
  print(
    f"Vol-wgt Qice={Qice[j0,i0]:.3e} J/m3, " 
    f"E_cell={Eice[j0,i0]:.3e} J/m2, "
    f"E_ice={Eice_ice[j0,i0]:.3e}\n"
  )


# Mask where no ice
A2d = np.where(aice > 0, A2d, 100.)
A2d = np.where(HH >= 0, np.nan, A2d)

#clrmp = mclrmps.colormap_conc()
#clrmp = mclrmps.colormap_ice_thkn()
#clrmp = mclrmps.colormap_ice_thkn()
clrmp = mclrmps.colormap_cold_warm()
rmin = -5.
rmax = 0.
clrmp.set_over([1,1,1])
clrmp.set_bad(color=[0.2, 0.2, 0.2])

if ilr is not None:
  sttl = f"{stxt}*1e8 ilayer={ilr}:  CICE6 restart {flrst_in}"
else:
  sttl = f"{stxt}*1e8 all layers:  CICE6 restart {flrst_in}"


if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)


xh, yh = m(hlon,hlat)

plt.ion()

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
img1 = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
# Overlay land
#ax1.pcolormesh(xh, yh, land_overlay, cmap=land_cmap, zorder=5)

ax1.contour(xh, yh, aice, [0.15], linestyles='solid', colors=[(0.4,0.4,0.4)], linewidths=1)
ax1.contour(xh, yh, HH, [0.], linestyles='solid', colors=[(0.,0.,0.)], linewidths=1)

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



