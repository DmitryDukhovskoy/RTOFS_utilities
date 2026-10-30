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
import xarray as xr
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
#flmom_in = '20240701.000000.MOM.res.nc'  # MOM restart file or specify in YAML 
#flmom_in = '20240701.000000.MOM.res.iconc.nc'  # MOM restart file or specify in YAML 
pthrest_in = '/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_sfs_C192mx025/ocean'
flmom_in = '20250701.000000.MOM.res.nc'

#pthrest_in = '/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/RUNDIRS/restart_sfs_C192mx025/ocean/mom6_Tfrz_ML_GDAS'
#flmom_in = 'gdas.t18z.inst.f006.snphys.20250701.Tfrz.nc'
#pthrest_in = '/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/GFSv17/ocean.20250630'
#flmom_in = 'gdas.t18z.inst.f006.nc'


fyaml = 'cice6rest_files_SFS.yaml'
vartmp = 'Temp'      # variable name in the restart file

# 3D fields to check in restart file:
#varlr  = 'eta'   # interface depth
#varsal = 'Salt'  # salinity
VARS = ['u','v','h','Salt','Temp']


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
if pthrest_in is None:
  pthrest_in  = config_rest["mom6_paths"]["pth_in"]

# Get MOM6 grid
pthgrid = config_rest["mom6_paths"]["pth_grid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

with xr.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, 1., -HH)
jdm, idm = HH.shape 

LON, LAT = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')


# MOM input restart:
dflmom_in = os.path.join(pthrest_in, flmom_in) 
print(f"Reading restart: {dflmom_in}")


#ds_in = xarray.open_dataset(dflmom_in, decode_times=False)

FLD3D = []
with xr.open_dataset(dflmom_in, decode_times=False) as ds_in: 
  for varnm in VARS:
    A3d = ds_in[varnm].values.squeeze()
    FLD3D.append(A3d)

#ds_in.close()


def check_lrs(A3d, fname):
  nlrs = 75
  for ilr in range(nlrs):
    A2d = A3d[ilr,:,:]

    if np.isnan(A2d).any():
      print(f"  {fname} Lr={ilr+1} NaNs")
    else:
      print(f"  {fname} Lr={ilr+1} ok")

print(f"MOM6 restart: {dflmom_in}")
# Check for NaNs:
for icc, varnm in enumerate(VARS):
  A3d = FLD3D[icc]
  if np.isnan(A3d).any():
    print(f"{varnm} NaNs, checking layers ...")
    check_lrs(A3d, varnm)
  else:
    print(f"{varnm} is ok")
  print(" ")




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


  # ======= DEBUG =======
  with xr.open_dataset(dflmom_in, decode_times=False) as ds_in:
    lrh = ds_in["h"].values.squeeze()
    if "ave_ssh" in ds_in:
      ssh = ds_in["ave_ssh"].values.squeeze()
    else:
      ssh = ds_in["sfc"].values.squeeze()
    U = ds_in["u"].values.squeeze()
    T = ds_in["Temp"].values.squeeze()
    S = ds_in["Salt"].values.squeeze()

  # Deep
  i0 = 500
  j0 = 1050
  
  # Shallow:
  i0=261
  j0=958

  kdm = lrh.shape[0]
  dLr = lrh[:,j0,i0]
  temp = T[:,j0,i0]
  sal = S[:,j0,i0] 

  zz = np.zeros((kdm+1))
  zm = np.zeros((kdm))
  zz[0] = ssh[j0,i0]
  zz[1:] = zz[0] - np.cumsum(dLr)
  zm = zz[:kdm] - dLr/2

  for kk in range(kdm):
      print(f"lrh: {dLr[kk]:.7f} u: {U[kk,j0,i0]:.7f}")

  print(f"Sum lr thkn: {np.nansum(-lrh[:,j0,i0]):.7f}, btm depth: {HH[j0,i0]:.7f}, ssh={ssh[j0,i0]:.7f}")
  #print(f"Sum lr thkn + ssh = {(np.nansum(-lrh[:,j0,i0]) + ssh[j0,i0]):.7f}")

  plt.clf()
  ax1 = plt.axes([0.08, 0.1, 0.4, 0.8])
  ax1.plot(temp, zm, '.-')
  
  ax2 = plt.axes([0.5, 0.1, 0.4, 0.8])
  ax2.plot(sal, zm, '.-')
  



