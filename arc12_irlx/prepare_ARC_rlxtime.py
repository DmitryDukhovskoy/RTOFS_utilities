"""
  Prepare relaxation rate (s-1) 
  for ARC12 
  Relaxation is uniform over the domain

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
import argparse

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
from mod_utils_fig import bottom_text
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--trlx", help="relaxation time scale, hours", type=float, required=True)
args = parser.parse_args()

rate_max_hrs = args.trlx if args.trlx else None
Irate_max_sec = 1./(rate_max_hrs*3600.)  # relaxation rate, s-1

lat_stop = 60. # no relaxation south of this lat

f_save    = False         # Save netcdf relax file
check_rlx = True         # Plot relaxation field

rlx_name = 'relax_rate' # name of the variable, should be the same in the SIS_input


pthdata = '/work/Dmitry.Dukhovskoy/ARC12/irlx'
pthtopo = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'

dflarc  = os.path.join(pthtopo,'ocean_hgrid.nc')
hlon, hlat = mmom6.read_mom6grid(dflarc, grdpnt='hgrid')

dtopo = os.path.join(pthtopo,'ocean_topog.nc')
ds_topo = xarray.open_dataset(dtopo)
HH = -ds_topo['depth'].data
jdm, idm = HH.shape

assert HH[300,200] < 0., f'Check sign of topography, ocean pnts should be < 0'
RLX = np.where(HH>=0, 0.0, Irate_max_sec)  # relax rate, s-1

# No relaxation south of lat_stop:
RLX[hlat<lat_stop] = 0.0

# Rewrite old relax file:
varnm = 'Idamp'
frlx_old = 'damping_full_t_30.nc'
dfrlx_old = os.path.join(pthdata,frlx_old)
ds_irlx = xarray.open_dataset(dfrlx_old)
darray = ds_irlx[varnm]
darray.values[:,:] = RLX
ds_irlx[varnm] = darray

btx = 'prepare_ARC_rlxtime.py'
cwd = os.getcwd()
ds_irlx.attrs["info"]=f"Time relaxation scale for ARC12, no nudging south of {lat_stop:.1f}N"
ds_irlx.attrs["history"] = f"Created {cwd}/{btx}"

if f_save:
  pthoutp = pthdata
  os.makedirs(pthoutp, exist_ok=True)

  flnew = f'ARC12_rlx_timescale_{int(rate_max_hrs):03d}hrs.nc'

  dflout = os.path.join(pthoutp, flnew)
  print(f"Saving to {dflout}")
  ds_irlx.to_netcdf(
    dflout,
    format='NETCDF3_64BIT',
    engine='netcdf4',
  )

ds_irlx.close()



if check_rlx:
  # For checking, relaxation time, hrs:
  RLXHR = RLX.copy()
  RLXHR = np.where(RLXHR==0., np.nan, RLXHR)
  RLXHR = 1./RLXHR * 1/3600.
  max_rlx = np.nanmax(RLXHR)

  cff = 1.e7
  AP = RLX.copy()*cff
  AP = np.where(HH>=0., np.nan, AP)

  clrmp = mclrmps.colormap_temp2()
  clrmp = mclrmps.colormap_conc()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin=0
  rmax=3.5

  latcntrs = [x for x in range(40,90,10)]
  loncntrs1 = [x for x in range(10,350,10)]
  loncntrs2 = [x for x in range(-10,10,10)]
  hlon1 = np.where(hlon<0., hlon+360., hlon)
  hlon1 = np.where(hlon1 > 350., np.nan, hlon1)
  hlon1[hlat>80.] = np.nan
  hlon2 = np.where(hlon > 150., np.nan, hlon)
  hlon2 = np.where(hlon < -50., np.nan, hlon2)
  hlon2[hlat>80.] = np.nan 


  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  img = ax1.pcolormesh(AP, cmap=clrmp, vmin=rmin, vmax=rmax)
  #ax1.contour(HH,[0],linestyles='solid', linewidths=1, colors=[(0., 0., 0.)])
  ax1.axis('scaled')
  ax1.contour(hlat, latcntrs, linestyles='solid', linewidths=1, colors=[(0.5, 0.5, 0.5)]) 
  ax1.contour(hlon1, loncntrs1, linestyles='solid', linewidths=1, colors=[(0.5, 0.5, 0.5)])
  ax1.contour(hlon2, loncntrs2, linestyles='solid', linewidths=1, colors=[(0.5, 0.5, 0.5)])

  sttl = (f'Relaxation rate (s-1), strongest rlx {max_rlx:.1f} hrs')
  ax1.set_title(sttl)

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

  btx = 'plot_rlxtime.py'
  bottom_text(btx, pos=[0.2, 0.01])



