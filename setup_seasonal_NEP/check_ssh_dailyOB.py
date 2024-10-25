"""
  Check daily or monthly SSH OB created from SPEAR runs
"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
#import pickle
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
from boundary import Segment
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
from mod_utils_fig import bottom_text

# Climatology derived for these years, started at mstart
# Inidicate start of the SPEAR forecast:
ens_spear  = 1      # ens run used to create OB's
yr_start   = 1993
mo_start   = 4
dnmb_start = mtime.datenum([yr_start,mo_start,1])
dv_start   = mtime.datevec(dnmb_start)
itime  = 0  # time index to plot
varnm  = 'ssh'

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)

# MOM6 NEP topo/grid:
run_name   = 'seasonal_fcst_daily'
pthtopo    = gridfls['MOM6_NEP'][run_name]['pthgrid']
fgrid      = gridfls['MOM6_NEP'][run_name]['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"][run_name]["ftopo"]
outdir     = gridfls['MOM6_NEP'][run_name]['pthoutp']
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
dftopo_mom  = os.path.join(pthtopo, ftopo_mom)
LONM, LATM  = mmom6.read_mom6grid(dfgrid_mom)
HHM         = mmom6.read_mom6depth(dftopo_mom)

segments = [ Segment(1, 'north', hgrid, output_dir=outdir),
             Segment(2, 'east',  hgrid, output_dir=outdir),
             Segment(3, 'south', hgrid, output_dir=outdir),
             Segment(4, 'west',  hgrid, output_dir=outdir)]

nOB = len(segments)

# Load mapping indices exist, gmapi:
dirgmapi = config['filesystem']['spear_mom_gmapi']
flgmaph  = f'spear2mom_NEP_OB_gmapi_hpnt.nc'
flgmapu  = f'spear2mom_NEP_OB_gmapi_upnt.nc'
flgmapv  = f'spear2mom_NEP_OB_gmapi_vpnt.nc'
dflgmaph = os.path.join(dirgmapi, flgmaph)
dflgmapu = os.path.join(dirgmapi, flgmapu)
dflgmapv = os.path.join(dirgmapi, flgmapv)
# h-point indices
dsh = xarray.open_dataset(dflgmaph)

spear_dir = config['filesystem']['nep_spear_subset'].\
                   format(year=dv_start[0], ens=ens_spear)


dltx = 0.07
dlty = 0.1
dx = 0.4
dy = 0.3
xl = 0.06
yb = 0.15
FPOS = [[xl, yb+dy+dlty, dx, dy],
        [xl+dx+dltx, yb+dy+dlty, dx, dy],
        [xl, yb, dx, dy],
        [xl+dx+dltx, yb, dx, dy]]

# Ssh - 1D sections
# Load ssh daily fields for NEP subset SPEAR 
flnm_spear = f'NEP_spear_{dv_start[0]}{dv_start[1]:02d}.ssh_daily.nc'
ds_spear = mutob.read_spear_output(spear_dir, varnm, flnm_spear, fzint=True)

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()

for isgm in range(nOB):
  nsgm  = isgm+1
  print(f'Processing ssh OB segment={nsgm}')
  INDX   = dsh[f'indx_segm{nsgm:03d}'].data
  JNDX   = dsh[f'jndx_segm{nsgm:03d}'].data
  dset_segm = mutob.segm_topo(nsgm, HHM, hgrid)
  distOB   = dset_segm['dist_supergrid'].data
  Xbtm     = dset_segm['dist_grid'].data
  Hbtm     = dset_segm['topo_segm'].data
  Hbtm     = np.where(Hbtm > 0, 0., Hbtm)
  segm_nm  = dset_segm['segm_name'].data[0]

#  # Spatial interpolation from SPEAR --> NEP OB supergrid
#  dset   = mutob.derive_obsegm_ssh(hgrid, ds, segments, isgm, INDX, JNDX, time_steps=time_days)

# Segment coordinates SPEAR
  dset_sgmspear = segments[isgm]
  xOB_spear = dset_sgmspear.coords.lon.data
  yOB_spear = dset_sgmspear.coords.lat.data
  npnts     = len(xOB_spear)
  nx_spear  = dset_sgmspear.nx
  ny_spear  = dset_sgmspear.ny

  # Calculate distance along the OB segment:
  distOB_spear, _ = mutob.calculate_dist_section(xOB_spear, yOB_spear)

  ssh_spear = ds_spear['ssh'].isel(time=itime).data
  sshOB_spear = ssh_spear[JNDX[:,0],INDX[:,0]].squeeze()

# Debugging ssh interpolation:
  f_debug = False
  if f_debug:
    import mod_mom6 as mom6util
    A = ssh_spear.copy() 
    # Fill missing values (bottom/land):
    dmm =  mom6util.fill_land3d(A)

    # 4 vertices chosen for SSH interpolation onto MOM grid:
    sshS1 = ssh_spear[JNDX[:,0],INDX[:,0]].squeeze()
    sshS2 = ssh_spear[JNDX[:,1],INDX[:,1]].squeeze()
    sshS3 = ssh_spear[JNDX[:,2],INDX[:,2]].squeeze()
    sshS4 = ssh_spear[JNDX[:,3],INDX[:,3]].squeeze()

    # Filled land OB segment of ssh at 4 indices for interpolation onto MOM grid
    sshF1 = dmm[JNDX[:,0],INDX[:,0]].squeeze()
    sshF2 = dmm[JNDX[:,1],INDX[:,1]].squeeze()
    sshF3 = dmm[JNDX[:,2],INDX[:,2]].squeeze()
    sshF4 = dmm[JNDX[:,3],INDX[:,3]].squeeze()

    plt.ion()
    fig1 = plt.figure(1,figsize=(9,8))

    fig1.clf()
    ax1  = plt.axes([0.1, 0.3, 0.8, 0.6])
    ax1.plot(sshF1,'-')
    ax1.plot(sshS1)
 

  # Read saved and interpolated OBs:
  date_init = f'{dv_start[0]}{dv_start[1]:02d}{dv_start[2]:02d}'
  pthoutp = gridfls['MOM6_NEP'][run_name]['pthoutp']
  fobc_out = os.path.join(pthoutp,f'OBCs_spear_daily_init{date_init}_e{ens_spear:02d}.nc')
  dsetOB = xarray.open_dataset(fobc_out) 
  sshOB  = dsetOB[f'zos_segment_{nsgm:03d}'].isel(time=itime).data.squeeze()
  xOB    = dsetOB[f'lon_segment_{nsgm:03d}'].data
  yOB    = dsetOB[f'lat_segment_{nsgm:03d}'].data
#  distOB, _ = mutob.calculate_dist_section(xOB, yOB)

  sttl = f'SSH SPEAR and OB file,  OB={segm_nm} F/cast day={itime+1}'

  axpos = FPOS[isgm]
  ax1 = plt.axes(axpos)
  ax1.plot(distOB_spear, sshOB_spear)
  ax1.plot(distOB, sshOB)
  # Show land:
  if np.max(Hbtm) >=0.0:
    Iland = np.where(Hbtm >= 0.0)[0]
    Yland = np.zeros((len(Iland)))
    ax1.plot(Xbtm[Iland],Yland,'b-')
  ax1.grid('on')
  ax1.set_title(sttl)
  ax1.set_xlabel('Dist., km')

btx='check_ssh_dailyOB.py'
bottom_text(btx)



