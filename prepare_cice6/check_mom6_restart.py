"""
  After correcting CICE6 initial fields
  Updated MOM6 surface T:
  under sea ice T = Tfreeze(S)
  for aice < 1: weighted update 

   Need MOM restart with surf T and CICE6 restart with corrected aice
   to use for updating SST under the sea ice


Check NaNs in the MOM6 restart files - this will cause model to blow up
444: FATAL from PE     0: NaN in input field of reproducing_EFP_sum(_2d).


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
import mod_swstate as msws
import mod_mom6 as mmom6
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

# YAML is used if input file name or vartmp is None
flmom_in = '20240701.000000.MOM.res.nc'  # MOM restart file or specify in YAML 
fyaml = 'cice6rest_files_SFS.yaml'
vartmp = 'Temp'      # variable name in the restart file


print(f"Reading YAML with restart info: {fyaml}\n")
with open(fyaml) as ff:
  config_rest = safe_load(ff)

if vartmp is None:
  vartmp = config_rest["mom6_names"]["vartmp_temp"]
  assert vartmp is not None, "vartmp is not provided"

if flmom_in is None:
  flmom_in = config_rest["mom6_names"]["flmom_in"]
  assert flmom_in is not None, "flmom_in is not provided"


# MOM restart path: input and output restart dir
pthrest_in  = config_rest["mom6_paths"]["pth_in"]

# Get MOM6 grid
pthgrid = config_rest["mom6_paths"]["pth_grid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape 

LON, LAT = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')


# MOM input restart:
dflmom_in = os.path.join(pthrest_in, flmom_in) 
print(f"Reading restart: {dflmom_in}")
ds_in = xarray.open_dataset(dflmom_in, decode_times=False)

# 3D fields:
varlr  = 'eta'   # interface depth
varsal = 'Salt'  # salinity

T3d  = ds_in[vartmp].isel(Time=0).data  # 3D pot temp
L3d  = ds_in[varlr].isel(Time=0).data  # interf. depths, lr0 = surf interf
S3d  = ds_in[varsal].isel(Time=0).data  

ds_in.close()

# Check signs:
assert np.nanmin(HH) < 0, f"Bottom topo has to be negative"
assert np.nanmin(lrthk_in[1,:,:]) < 0, f"Layer interf depths expected <0"

def check_lrs(A3d, fname)
  nlrs = 75
  for ilr in range(nlrs):
    A2d = A3d[ilr,:,:]

    if np.isnan(A2d).any():
      print(f"  {fname} Lr={ilr+1} NaNs")
    else:
      print(f"  {fname} Lr={ilr+1} ok")



# Check for NaNs:
if np.isnan(T3d).any():
  fname="T3d"
  print(f"{fname} NaNs, checking layers ...")
  check_lrs(T3d, fname)




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





