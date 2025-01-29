"""
  Create inverse relaxation e-folding time scale
  Time scale can vary spatially allowing different relaxation
  rates for different parts of the domain

"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
from yaml import safe_load

PPTHN = '/home/Dmitry.Dukhovskoy/python'
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc

# Select max and min relaxation time scalres, hrs
# max relaxation - strongest, typically along the OBs
# min relaxation - somewhere in the domain where sea ice presents
#                  e.g. Bering strait
# rx_min, ry_min - approximate # of i, j pnts from the ice OBs
#                  i.e., from i=imax to Ber. Str. (342-200)
# relaxation time scales will be going to 0 away from the ice OBs
rate_max_hrs = 2.  # max relaxation time
rate_min_hrs = 96. # lower relax. 
rx_min = 140.       # approximate range for rate_min from rate_max along X
ry_min = 145.       # approximate range for rate_min from rate_max along Y
f_save = True
rlx_name = 'relax_rate' # name of the variable, should be the same in the SIS_input

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

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
# Hgrid lon. lat:
hlon, hlat  = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

# Basis fns
rate_fract = rate_max_hrs / rate_min_hrs
gauss_dnm_x = -(rx_min**2)/np.log(rate_fract)
gauss_dnm_y = -(ry_min**2)/np.log(rate_fract)

# 2D Gaussian:
nbndry = 7  # number of near-boundary points to keep max relaxataion
PSI = np.zeros((jdm,idm))
for jj in range(jdm):
  for ii in range(idm):
#    ay = gauss_dnm_y*(ii/(idm-rx_min))**2  
#    ax = gauss_dnm_x*(jj/(jdm-ry_min))**2
    ay = gauss_dnm_y*(ii/50.)**2  
    ax = gauss_dnm_x*(jj/100)**2
    if ax==0:
      ax = 1.e-8
    if ay == 0:
      ay = 1.e-8
    PSI[jj,ii] = np.exp(-(ii+1-idm)**2/ax -(jj+1-jdm)**2/ay)
    if ii >= idm-nbndry:
      PSI[jj,ii] = 1.
    if jj >= jdm-nbndry:
      PSI[jj,ii] = 1.
#Psi_i = np.exp(-(X-max(X))**2/gauss_dnm)
#Psi_j = np.expand_dims(Psi_j, axis=0)
#PSIJ = np.tile(Psi_j.transpose(),(1,idm))
#PSII = np.tile(Psi_i,(jdm,1))
#PSI = mmisc.box_fltr(PSI, nbx=25)

RLXIS = 1./(rate_max_hrs*3600)*PSI # relaxation, 1/sec
RLXIS[:,-1] = 1./(rate_max_hrs*3600)
RLXIS[-1,:] = 1./(rate_max_hrs*3600)

# Add land mask and southern domains = 0
RLXIS = np.where(HH>=0, 0.0, RLXIS)
RLXIS[:585,:] = 0.
RLXIS[:,:187] = 0.

# For checking, relaxation time, hrs:
RLXHR = RLXIS.copy()
RLXHR = np.where(RLXHR==0., np.nan, RLXHR)
RLXHR = 1./RLXHR * 1/3600.

# Write relax time scale:
dflstat = os.path.join(pthtopo, 'ocean_static.nc')
ds_stat = xarray.open_dataset(dflstat)
ds_rlx  = ds_stat['wet'].copy()
ds_rlx.name = rlx_name
ds_rlx *= RLXIS
ds_rlx = ds_rlx.to_dataset()
ds_rlx[rlx_name].attrs['units'] = 's-1'
ds_rlx[rlx_name].attrs['cell_method'] = 'time: point'


if f_save:
  encoding = {rlx_name: {'_FillValue': None}}
  flout = f'relax_rate_{int(rate_max_hrs):03d}hrs.nc'
  pthsis = gridfls['MOM6_NEP'][run_name]['pthsis']
  dflout = os.path.join(pthsis, flout)

  print(f'Saving SIS2 relaxation time scale --> {dflout}')
  ds_rlx.to_netcdf(
      dflout,
      format='NETCDF3_64BIT',
      engine='netcdf4',
      encoding=encoding
  )


check_rlx = False
if check_rlx:
  plt.ion()

  clrmp = mclrmps.colormap_temp2()
  rmin = 2.
  rmax = 20.

  # Stereographic Map projection:
  from mpl_toolkits.basemap import Basemap, cm
  m = Basemap(width=5000*1.e3,height=5000*1.e3, resolution='l',\
              projection='stere', lat_ts=50, lat_0=62, lon_0=-165)

  xR, yR = m(hlon, hlat)

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = ax1.pcolormesh(xR, yR, RLXHR, cmap=clrmp, vmin=rmin, vmax=rmax)
#  img = ax1.pcolormesh(RLXHR, cmap=clrmp)

  ax1.set_title('Relaxation time, hrs')

  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  # extend: min, max, both
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'relax_timescale.py'
  bottom_text(btx, pos=[0.2, 0.01])


