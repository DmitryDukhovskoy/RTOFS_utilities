"""
  Check restart CICE6
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
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
  
rest_date = 20250103
rest_hr   = 0
regn = 'north'
#flrst = 'cice_model.res.20250103.00.iconc.snow.nc'
#flrst = 'cice_model.res.20250103.00.iconc.nc'
#flrst = '20250103.030000.cice_model.res.nc'       # original ice restart
#flrst = 'cice_model.res.20250103.00.iconc_thkn.snow.nc'
flrst = 'cice_restart.20250103.00.iconc_ithkn.hsnow.snphys.nc' 

# ithkn - mean ice thkn per m2 of grid cell area
# ithkn_ice - mean ice thkn per m2 of ice area
# hsnow - mean snow depth per grid cell area
# hsnow_ice - mean snow depth per ice area 
parser = argparse.ArgumentParser()
parser.add_argument("--flrst", help=f"restart file name, default={flrst}", type=str)
parser.add_argument("--regn", help=f"region to plot, default={regn}", type=str)
parser.add_argument("--fplt", help="Field to plot", 
                   choices=['iconc','ithkn','ithkn_ice','hsnow','hsnow_ice','none'], type=str)
args = parser.parse_args()

flrst = args.flrst if args.flrst else flrst
regn  = args.regn if args.regn else regn
fld_plt = args.fplt if args.fplt else None
if fld_plt == 'none':
  fld_plt = None

#rest_date     = args.rdate if args.rdate else rest_date
#rest_hr       = args.rhr if args.rhr else rest_hr
#rest_date_out = args.rdate_out if args.rdate_out else rest_date
#rest_hr_out   = args.rhr_out if args.rhr_out else rest_hr

# Get date numbers:
# Input restart file
#dnmbR = mtime.rdate2datenum(rest_date*100+rest_hr)  # restart day nmb
#yrR,mmR,ddR,hrR = mtime.datevec(dnmbR, round_hrs=True)[:4]
#nsecR = hrR*3600

# Dates of the output fields in the new restart:
#dnmbN = mtime.rdate2datenum(rest_date_out*100+rest_hr_out)
#yrN, mmN, ddN, hrN = mtime.datevec(dnmbN, round_hrs=True)[:4]
#nsecN = hrN*3600 

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
# Check ice_in what ITD is used
# 0.00, 0.64, 1.39, 2.47, 4.57
hicat = np.array([0., 0.64, 1.39, 2.47, 4.57, 1000.])


fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

#flrst = f"cice_model.res.{yrN}{mmN:02d}{ddN:02d}.{hrN:02d}.iconc.nc"
#pthrest = os.path.join(pths_ufs[node_nm]["MOM6"]["pthrest"],'new')
#pthrest = os.path.join(pths_ufs[node_nm]["MOM6"]["pthrest"])
pthrest = os.path.join(pths_ufs[node_nm]["MOM6"]["pthrest"],'cice6_global')
#pthrest = os.path.join(pths_ufs[node_nm]["MOM6"]["pthrest"],'cice6_north')
dflrst = os.path.join(pthrest,flrst)

print(f"Reading {dflrst}")

# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

LON, LAT = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')


ds_rst = xarray.open_dataset(dflrst)
aicen = ds_rst['aicen'].data
vicen = ds_rst['vicen'].data
vsnon = ds_rst['vsnon'].data
ncat = vsnon.shape[0]

# Total ice vol per grid cell area, m3_ice / m2_cell:
vice = np.sum(vicen, axis=0).squeeze()

# Aggreg iconc:
aice = np.sum(aicen, axis=0).squeeze()

# grid-cell mean ice thickness:
hice_cell = vice.copy()

# ice thickness over ice only, m3_ice / m2_ice:
hice_ice = np.divide(vice, aice, out=np.zeros_like(aice), where=aice != 0)

# grid-cell snow depth:
hsnow_cell = np.sum(vsnon, axis=0)

# mean snow depth over ice area:
hsnow_ice = np.divide(hsnow_cell, aice, out=np.zeros_like(aice), where=aice != 0)

# Check hice(n) as it is caclulated in icepack_therm_vertical.F90
# hice(n) = vice(n) / aice(n) 
print(' =========  ICE  =========')
for k in range(1,ncat+1):
  aice_n = aicen[k-1,:].squeeze()
  vice_n = vicen[k-1,:].squeeze()
  hice_n = np.divide(vice_n, aice_n, out=np.zeros_like(aice), where=aice_n != 0)
  Mvalid = (~np.isnan(hice_n)) & (hice_n > puny)  # valid values
  hice_min = np.where(Mvalid, hice_n, 1.e20)
  hice_max = np.where(Mvalid, hice_n, -1.e20)

  jmin, imin = np.unravel_index(hice_min.argmin(), hice_n.shape)
  jmax, imax = np.unravel_index(hice_max.argmax(), hice_n.shape)

  hbin_min = hicat[k-1]
  hbin_max = hicat[k]
  Mcrct  = (hice_n >= hbin_min) & (hice_n < hbin_max) # correct
  Merr   = (~Mcrct) & (Mvalid)
  Jerr, Ierr = np.where(Merr)

  # ice thicknesses within the cats cannot cross-over!!!
  # Check hice[k+1] > hice[k], see icepack_therm_itd.F90 ITD thermodyn 
  Jdh, Idh = [], []
  if k > 1:
    dlt_hice = hice_n - hice_km1
    dlt_hice = np.where(Mvalid, dlt_hice, np.inf)
    Merr_dh = (dlt_hice <= 0)
    Jdh, Idh = np.where(Merr_dh)

  hice_km1 = hice_n.copy()

  print(f"Cat {k}: ")
  print(f"  j={jmin}, i={imin}, min hice(n): {hice_n[jmin,imin]}, "+\
        f"aice(n): {aice_n[jmin,imin]}, vice(n): {vice_n[jmin,imin]}")
  print(f"  j={jmax}, i={imax}, max hice(n): {hice_n[jmax,imax]}, "+\
        f"aice(n): {aice_n[jmax,imax]}, vice(n): {vice_n[jmax,imax]}")
  print(f"  found {len(Jerr)} points violating: {hbin_min:.3f} <= hice < {hbin_max:.3f}")
  print(f"  found {len(Jdh)} points violating hice[k] > hice[k-1]") 

 
print('\n =========  SNOW =========')
for k in range(1,ncat+1):
  aice_n = aicen[k-1,:].squeeze()
  vsno_n = vsnon[k-1,:].squeeze()
  hsno_n = np.divide(vsno_n, aice_n, out=np.zeros_like(aice), where=aice_n != 0)
  jmin, imin = np.unravel_index(hsno_n.argmin(), hsno_n.shape)
  jmax, imax = np.unravel_index(hsno_n.argmax(), hsno_n.shape)
  print(f"Cat {k}: ")
  print(f"  j={jmin}, i={imin}, min hsnow(n): {np.nanmin(hsno_n)}, "+\
        f"aice(n): {aice_n[jmin,imin]}, vsno(n): {vsno_n[jmin,imin]}")
  print(f"  j={jmax}, i={imax}, max hsnow(n): {np.nanmax(hsno_n)}, "+\
        f"aice(n): {aice_n[jmax,imax]}, vsno(n): {vsno_n[jmax,imax]}")

print('\n =========  ICE CONCENTRATION =========')
jmax, imax = np.unravel_index(aice.argmax(), aice.shape)
print(f"Max ice conc: {aice[jmax,imax]}, j={jmax}, i={imax}")
for k in range(1,ncat+1):
  print(f"Cat {k}: ")
  aice_n = aicen[k-1,:].squeeze()
  print(f"  j={jmax}, i={imax}, aicen(n): aice(n): {aice_n[jmax,imax]}")

f_plt = fld_plt is not None

if f_plt:
  j0 = i0 = None

  print(" Plotting ...")
  if fld_plt == 'iconc':
    clrmp = mclrmps.colormap_conc()
    rmin = 0.
    rmax = 1.
    A2d = aice.copy()

  if fld_plt == 'ithkn':
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 3.
    A2d = hice_cell.copy()

  if fld_plt == 'ithkn_ice':
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 3.
    A2d = hice_ice.copy()

  if fld_plt == 'hsnow':
    clrmp = mclrmps.colormap_temp()
    rmin = 0.
    rmax = 0.4
    clrmp.set_under(color=[1,1,1])
    A2d = hsnow_tot.copy()
    
  if fld_plt == 'hsnow_ice':
    clrmp = mclrmps.colormap_temp()
    rmin = 0.
    rmax = 0.4
    clrmp.set_under(color=[1,1,1])
    A2d = hsnow_tot.copy()
    
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  hlon = LON
  hlat = LAT

  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
    parallels = np.arange(-80,-10,10.)
    meridians = np.arange(-360,359.,45.)
  elif regn == 'north':
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
    parallels = np.arange(50, 90, 5)
    meridians = np.arange(-360, 359., 45.)

  xh, yh = m(hlon,hlat) # GFS coords

  plt.ion()
  fig1 = plt.figure(1, figsize=(8,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.13, 0.8, 0.8])

  m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
  img1 = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)

  ax1.contour(xh,yh,HH,[0], linestyles='solid', colors=[(0.,0.,0.)], linewidths=1)

  ax1.set_title(f'{flrst}, {fld_plt}\n {pthrest}')

  # Plot pnt:
  if j0 is not None:
    x0 = hlon[j0,i0]
    y0 = hlat[j0,i0]
    xh0 = xh[j0,i0]
    yh0 = yh[j0,i0]
    ax1.plot(xh0,yh0,'o')

  # Colorbars
  ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
  clb = plt.colorbar(img1, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'check_cice6_restart.py'
  bottom_text(btx, pos = [0.02,0.02])




