"""
  After correcting CICE6 initial fields
  Updated MOM6 surface T:
  under sea ice T = Tfreeze(S)
  for aice < 1: weighted update 

  The algorithm checks water column stability after SST adjustment 

   Need MOM restart with surf T and CICE6 restart with corrected aice
   to use for updating SST under the sea ice

  Note: if "eta" is not in the restart options are:
  (1) Use warm start (restart) to skip reading "eta" in MOM.res.nc
      that is specified in MOM_input

      THICKNESS_CONFIG = "file"       ! default = "uniform"
      THICKNESS_FILE = MOM.res.nc     !
      INTERFACE_IC_VAR = "eta"        ! default = "eta"

      To use warm restart change:
      in input.nml:  input_filename = 'r'  

  (2) Change to construct eta from layer depths
      THICKNESS_CONFIG = "thickness_file"       ! default = "uniform"
      THICKNESS_FILE = MOM.res.nc     !
      INTERFACE_IC_VAR = "h"     

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
import mod_swstate as msws
import mod_mom6 as mmom6
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

yrR = mmR = ddR = hrR = None
yrN = mmN = ddN = hrN = None

# Current code assumes 1 MOM restart file
# No date / time change 
# TODO: add logic for N MOM restart files: MOM_*.res.nc, MOM_*.res_1.nc, ... 
# YAML is used if input file name or vartmp is None
#flmom_in = '20240701.000000.MOM.res.nc'  # MOM restart file or specify in YAML 
flrst_out = None        # MOM modified restart

parser = argparse.ArgumentParser()
parser.add_argument("--flmom_in", help=f"MOM6 rest file name, original")
parser.add_argument("--flice", help="CICE6 restart file to use for MOM6 surf T", required=True)
parser.add_argument("--rdate_out", help="Restart date for modified MOM6 restart", required=True)
parser.add_argument(
  "--fyaml",
  help=f"YAML with local directories, filenames, restart dates",
  required=True
)
args = parser.parse_args()

flmom_in   = args.flmom_in if args.flmom_in else None 
fyaml      = args.fyaml
flrst_cice = args.flice
rdate_out  = args.rdate_out

print(f"Reading YAML with restart info: {fyaml}\n")
with open(fyaml) as ff:
  config_rest = safe_load(ff)

if flmom_in is None:
  flmom_in = config_rest["mom6_names"]["flrst_in"]
  
assert flmom_in is not None, f"flmom_in has to be provided either --flmom_in or in {fyaml}"


# MOM restart path: input and output restart dir
pthrest_in  = config_rest["mom6_paths"]["pth_in"]
pthrest_out = config_rest["mom6_paths"]["pth_out"]

if pthrest_out is None:
  pthrest_out = pthrest_in 

# Get MOM6 grid
pthgrid = config_rest["mom6_paths"]["pth_grid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

if np.nanmin(HH) >= 0:
  HH = np.where(HH < 1.e-20, 1, -HH)
  HH[~np.isfinite(HH)] = 1
jdm, idm = HH.shape 

LON, LAT = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

# CICE parameters:
puny      = 1.e-11

# CICE6 restart file:
pthcice_in = config_rest["cice_paths"]["pth_out"]
dflcice = os.path.join(pthcice_in, flrst_cice)
with xarray.open_dataset(dflcice) as ds_ice:
  aicen = ds_ice["aicen"].data  # partial ice area by cats

# Aggregated ice partial area:
aice = np.sum(aicen, axis=0).squeeze()

# MOM input restart:
dflmom_in = os.path.join(pthrest_in, flmom_in) 
print(f"Reading restart: {dflmom_in}")
ds_in = xarray.open_dataset(dflmom_in, decode_times=False)
ds_out = ds_in.copy(deep=True)

# Determine variable name, first use the one from yaml if provided:
vartmp = config_rest["mom6_names"]["varnm_temp"]
varlr  = config_rest["mom6_names"]["varnm_intrf"]  # interf depth or lr thicknesses
varsal = config_rest["mom6_names"]["varnm_sal"]

varlr  = 'eta' if varlr is None else varlr
varsal = 'Salt' if varsal is None else varsal

if vartmp is None:
  # Try temp - case insensitive
  vartmp = next((v for v in ds_in.data_vars if v.lower() == 'temp'), None)

for vv in [vartmp, varlr, varsal]:
  if vv is None or vv not in ds_in.variables:
    raise KeyError(
      f"variable {vv} not found in restart file, "
      f"provide a valid varnm in {fyaml}"
    )

# Check if varlr is lr thickness or lr depth:
long_name = ds_in[varlr].attrs.get("long_name")
intrf_depth = None
if long_name is not None and "layer thickn" in long_name.lower():
  intrf_depth = False

elif long_name is not None and "interface height" in long_name.lower():
  intrf_depth = True

else:
  # Fall back to variable name if long_name is missing or unrecognized
  intrf_depth = varlr.lower() == "eta"

  print(
      f"WARNING: {'No long_name attribute' if long_name is None else f'Unrecognized long_name: {long_name!r}'} "
      f"for variable '{varlr}',\n"
      f"decisions based on variable name: Interface depth={intrf_depth}"
  )

assert intrf_depth is not None, f"Failed to determine {varlr!r} type (interf depth or lr. thickness) in netcdf"

# Assumed only 1 time level in MOM restart:
temp_in = ds_in[vartmp].values.squeeze().copy()  # 3D pot temp
salt_in = ds_in[varsal].values.squeeze().copy()
lr_in   = ds_in[varlr].values.squeeze().copy() 

ds_in.close()

assert temp_in.ndim == 3, f"Check temp array, expected 3D, ndim={temp_in.ndim}"
assert salt_in.ndim == 3, f"Check temp array, expected 3D, ndim={salt_in.ndim}"


# Construct lr depths from lr thickness
if intrf_depth:
  lrthkn_in = None
  lrdpth_in = -np.abs(lr_in)
  assert np.nanmin(lrdpth_in[1, :, :]) < 0
else:
  lrthkn_in = lr_in
  lrdpth_in = None

# Check signs:
assert np.nanmin(HH) < 0, f"Bottom topo has to be negative"

# Compute freezing T:
S2d = salt_in[0,:,:].squeeze()
S2d[HH >= 0] = np.nan
Tfrz = msws.freezing_temp_linear(S2d,0)

# Update surf SST in MOM6:
T2d = temp_in[0,:,:].squeeze()

T2d_new = T2d*(1. - aice) + aice * Tfrz
T2d_new = np.where(aice < puny, T2d, T2d_new)

T3d_new = temp_in.copy()
T3d_new[0,:,:] = T2d_new
S3d_new = salt_in.copy()

def drho_lrdpth(ilr, S3d, T3d, aice, puny, HH, lrdpth_in):
  rho1 = msws.sw_dens0(S3d[ilr,:,:], T3d[ilr,:,:])
  rho2 = msws.sw_dens0(S3d[ilr+1,:,:], T3d[ilr+1,:,:])
  drho = rho1 - rho2   
  # Mask not ice regions:
  drho[aice < puny] = np.nan
  # Mask regions where layer depth is shallower than local bottom
  # HH < 0 for ocean grid points
  zzLr = lrdpth_in[ilr,:,:]  # upper interface
  drho[zzLr < HH] = np.nan 

  return drho

def drho_lrthkn(ilr, S3d, T3d, aice, puny, HH, lrthkn_in, dhmin=0.05):
  rho1 = msws.sw_dens0(S3d[ilr,:,:], T3d[ilr,:,:])
  rho2 = msws.sw_dens0(S3d[ilr+1,:,:], T3d[ilr+1,:,:])
  drho = rho1 - rho2   

  dz1 = lrthkn_in[ilr]
  dz2 = lrthkn_in[ilr+1]

  invalid = (
      (aice < puny)
      | (HH >= 0)
      | (~np.isfinite(dz1))
      | (~np.isfinite(dz2))
      | (dz1 < dhmin)
      | (dz2 < dhmin)
  )

  drho[invalid] = np.nan

  return drho

def mix_layers_lrdpth(ilr, lrdpth_in, HH, T3d_new, S3d_new, JJ, II):
  # Interface depths:
  zlr1 = lrdpth_in[ilr,:,:]
  zlr2 = lrdpth_in[ilr+1,:,:]
  zlr3 = lrdpth_in[ilr+2,:,:]
  
  # mask out below-bottom vanished layers
  # last interf = bottom depth
  zlr1 = np.where(zlr1 < HH, HH, zlr1)
  zlr1[HH >= 0] = np.nan
  zlr2 = np.where(zlr2 < HH, HH, zlr2)
  zlr2[HH >= 0] = np.nan
  zlr3 = np.where(zlr3 < HH, HH, zlr3)
  zlr3[HH >= 0] = np.nan

  # Lyaer thickness:
  dz1 = np.abs(zlr2 - zlr1)
  dz2 = np.abs(zlr3 - zlr2)
 
  # mix two layers where drho > eps:
  T1 = T3d_new[ilr,:,:]
  T2 = T3d_new[ilr+1,:,:]
  S1 = S3d_new[ilr,:,:]
  S2 = S3d_new[ilr+1,:,:]
  den = dz1 + dz2
  den = np.where(den == 0, np.nan, den)

  assert np.all(np.isfinite(dz1[JJ, II]))
  assert np.all(np.isfinite(dz2[JJ, II]))
  assert np.all(dz1[JJ, II] > 0)
  assert np.all(dz2[JJ, II] > 0)

  tmix =  (dz1*T1 + dz2*T2) / den
  smix =  (dz1*S1 + dz2*S2) / den

  T3d_new[ilr,JJ,II]   = tmix[JJ,II]
  T3d_new[ilr+1,JJ,II] = tmix[JJ,II]
  S3d_new[ilr,JJ,II]   = smix[JJ,II]
  S3d_new[ilr+1,JJ,II] = smix[JJ,II]

  return T3d_new, S3d_new
    
def mix_layers_lrthkn(ilr, lrthkn_in, HH, T3d_new, S3d_new, JJ, II, dhmin=0.05):
  # Lr thickness:
  dz1 = lrthkn_in[ilr,:,:]
  dz2 = lrthkn_in[ilr+1,:,:]
  
  # mix two layers where drho > eps:
  T1 = T3d_new[ilr,:,:]
  T2 = T3d_new[ilr+1,:,:]
  S1 = S3d_new[ilr,:,:]
  S2 = S3d_new[ilr+1,:,:]
  den = dz1 + dz2
  den = np.where(den < dhmin, np.nan, den)

  assert np.all(np.isfinite(dz1[JJ, II]))
  assert np.all(np.isfinite(dz2[JJ, II]))
  assert np.all(dz1[JJ, II] > 0)
  assert np.all(dz2[JJ, II] > 0)

  tmix =  (dz1*T1 + dz2*T2) / den
  smix =  (dz1*S1 + dz2*S2) / den

  T3d_new[ilr,JJ,II]   = tmix[JJ,II]
  T3d_new[ilr+1,JJ,II] = tmix[JJ,II]
  S3d_new[ilr,JJ,II]   = smix[JJ,II]
  S3d_new[ilr+1,JJ,II] = smix[JJ,II]

  return T3d_new, S3d_new

def compare_finite(name, original, modified):
  """
  Identify if new nans were introduced during
  SST adjustment
  """
  orig_finite = np.isfinite(original)
  mod_finite = np.isfinite(modified)

  introduced = orig_finite & ~mod_finite
  repaired = ~orig_finite & mod_finite

  print(
    f"{name}: introduced non-finite = {introduced.sum()}, "
    f"original non-finite made finite = {repaired.sum()}"
  )

  if np.any(introduced):
    inds = np.argwhere(introduced)
    raise ValueError(
        f"{name}: modification introduced non-finite values; "
        f"first 10 locations: {inds[:10].tolist()}"
    )

# Next: check water column stability in the upper ocean:
# Potential density ref = 0, assumed that only upper ocean needs adjustment
# Note that in the original fields, some layers can be unstable which seems to be ok for MOM6
# Therefore, do not mix all whatercolumn but only the upper few layers if needed
eps_drho = 1e-2    # threshold for delta rho
mixlrs = 3     # the upper layers to mix if unstable after Tfrz 
nlrs = salt_in.shape[0]
ilr = 0
max_iter = 50
print(f"Adjusting unstable layers under sea ice, max_iter={max_iter}")
for iter in range(max_iter):
  unstb = False

  print(f"iter = {iter+1}")
  for ilr in range(mixlrs):
    if intrf_depth:
      drho = drho_lrdpth(ilr, S3d_new, T3d_new, aice, puny, HH, lrdpth_in)
    else:
      drho = drho_lrthkn(ilr, S3d_new, T3d_new, aice, puny, HH, lrthkn_in)

    print(f"  Lr={ilr+1:02d}, max drho={np.nanmax(drho):.6f}")
 
    JJ, II = np.where(drho > eps_drho)
    if JJ.size:
      unstb = True

      # mix layers ilr and ilr+1
      if intrf_depth:
        T3d_new, S3d_new = mix_layers_lrdpth(ilr, lrdpth_in, HH, T3d_new, S3d_new, JJ, II)
      else:
        T3d_new, S3d_new = mix_layers_lrthkn(ilr, lrthkn_in, HH, T3d_new, S3d_new, JJ, II)

  if not unstb:
    print(f"Converged in {iter+1} iterations")
    break

if unstb:
  raise RuntimeError(f"Instability adjustment did not converge after {iter+1} iterations")

# Check if any  new nans has been created
compare_finite("Temperature", temp_in, T3d_new)
compare_finite("Salinity", salt_in, S3d_new)


def fill_bottom_nans(arr):
  """
  Extend the deepest finite T/S value downward through trailing
  NaNs or Infs in each water column.
  Assumes axis 0 is depth, surface -> bottom.
  Internal gaps and leading invalid values are not filled.
  """
  arr = arr.copy()

  if arr.ndim != 3:
      raise ValueError(f"Expected 3D array, got shape {arr.shape}")

  # Dim 0 is depth
  nz = arr.shape[0]

  for j in range(arr.shape[1]):
    for i in range(arr.shape[2]):
      col = arr[:, j, i]
      valid = np.flatnonzero(np.isfinite(col))

      if valid.size == 0:
        continue

      last = valid[-1]
      col[last + 1:] = col[last]

  return arr

# Check for NaNs below bottom
# Not really needed in MOM restart, already filled
if np.isnan(T3d_new).any():
  print("T3d: Filling NaNs in the bottom cells ")
  T3d_new = fill_bottom_nans(T3d_new)

if np.isnan(S3d_new).any():
  print("S3d: Filling Nans in the bottom cells")
  S3d_new = fill_bottom_nans(S3d_new)


T3d_new = np.expand_dims(T3d_new, axis=0)
S3d_new = np.expand_dims(S3d_new, axis=0)
ds_out[vartmp].values[:] = T3d_new
ds_out[varsal].values[:] = S3d_new

# Get rid of Nans in layer thicknesses:
# Not recommended
#if varlr.lower() == "h":
#  ds_out[varlr].values[:] = 1.e-3


ds_out.attrs.update({
    "title": f"MOM6 restart with adjusted SST under ice to Tfrz weighted by aice {flrst_cice}",
    "source": "correct_surfTwclmn_stab_mom6_restart.py",
})

# Save:flmom_in
# Construct output file name if not provided
# use cice_restart used for aice
# assumed *.res.[sfx].nc  where sfx is iconc, hsnow, etc
# if missing - control run
if flrst_out is None:
  base = os.path.splitext(flmom_in)[0]
  cice_parts = flrst_cice.split(".")
  cice_sfx = cice_parts[-2]
  if cice_sfx == 'res':
    cice_sfx = 'cntrl'
  flrst_out = f"{base}.{cice_sfx}.{rdate_out}.Tfrz.nc"
  #flrst_out = f"{base}.{rdate_out}.Tfrz.nc"

dflrst_out = os.path.join(pthrest_out, flrst_out)
print(f"Saving MOM6 restart --> {dflrst_out}")
ds_out.to_netcdf(dflrst_out, encoding={var: {'_FillValue': None} for var in ds_out.data_vars}, format='NETCDF3_64BIT')
ds_out.close()


f_chck = False
if f_chck:
  import mod_colormaps as mclrmps
  
  plt.ion()

  units = 'm'
  clrmp = mclrmps.colormap_conc()
  rmin = 0.
  rmax = 1.

  clrmp = mclrmps.colormap_cold_warm()
  rmin = -1.8
  rmax = 0.
  clrmp.set_over([0.9,0.9,0.9])

  # Plot diff. btw MOM6 surf T and Tfrz
  #clrmp = mclrmps.colormap_uv()
  clrmp = mclrmps.colormap_warm(clrS=[1,1,1])
  rmin = 0
  rmax = 1.8
  
  #dltT = T2d_new - T3d_new[0,:,:]
  dltT = T2d_new - Tfrz
  dltT = T3d_new[0,:,:] - Tfrz
  dltT = Tfrz_new - Tfrz

  dltT[HH>=0] = np.nan
  clrmp.set_under([0.9,0.0,0.9])

  clrmp.set_bad(color=[0.1, 0.1, 0.1])

  regn = 'north'

  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
    parallels = np.arange(-80,-10,10.)
    meridians = np.arange(-360,359.,45.)
  elif regn == 'north':
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
    parallels = np.arange(50, 90, 5)
    meridians = np.arange(-360, 359., 45.)

  hlon = LON
  hlat = LAT   

  xh, yh = m(hlon,hlat)

  sttl = f'{flmom_in}, MOM6 sst - Tfrz\nCICE6 rest: {flrst_cice}'

  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)

  #img = ax1.pcolormesh(xh, yh, aice, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')
  #img = ax1.pcolormesh(xh, yh, Tfrz, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')
  img = ax1.pcolormesh(xh, yh, dltT, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')
  ax1.contour(xh, yh, aice, [0.15], linestyles='solid', colors=[(0.2,1,0.8)], linewidths=1)

  ax1.set_title(sttl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  if rmin < 0:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  else:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')

  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)


  btx = 'correct_surfT_mom6_restart.py'
  bottom_text(btx, pos=[0.2, 0.01])

  #==== debug
  with xarray.open_dataset(dflmom_in, decode_times=False) as ds:
    eta = ds["eta"].values.squeeze()
    lrh = ds["h"].values.squeeze()
    ssh = ds["sfc"].values.squeeze()

  # Deep
  i0 = 500
  j0 = 1050
  
  # Shallow:
  i0=261
  j0=958

  kdm = lrh.shape[0]
  for kk in range(kdm+1):
    if kk < kdm:
      print(f"eta: {eta[kk,j0,i0]:.7f} lrh: {lrh[kk,j0,i0]:.7f} zz+1: {eta[kk,j0,i0]-lrh[kk,j0,i0]:.7f}")
    else:
      print(f"eta: {eta[kk,j0,i0]:.7f} lrh: -----")

  print(f"Sum lr thickness: {np.sum(lrh[:,j0,i0]):.7f}, btm depth: {HH[j0,i0]:.7f}, ssh={ssh[j0,i0]:.7f}")








