"""
  Plot relaxation field for ARC12

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
parser.add_argument("--mm", help="cal. month to plot: 1,..., 12, ...", type=int)
parser.add_argument("--varnm", help="field to plot: ithkn or iarea/iconc", type=str)
args = parser.parse_args()

if args.mm:
  MM = args.mm
if args.varnm:
  varnm = args.varnm
  if varnm=='iconc' or varnm=='iarea':
    varnm_nc = 'sic_rg'
  elif varnm=='ithkn' or varnm=='ithk':
    varnm_nc = 'sit_rg'

if varnm == 'iarea':
  varnm = 'iconc'

# Test point:
jF0 = 377
iF0 = 24
i0  = iF0-1
j0  = jF0-1


print(f'Plotting 1993/{MM} {varnm}')
pthdata = '/work/Dmitry.Dukhovskoy/ARC12/irlx'
pthtopo = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'
frlx = 'nudging_ice.nc'
dfrlx = os.path.join(pthdata,frlx)
ds_irlx = xarray.open_dataset(dfrlx)

itime = MM-1
A2d = ds_irlx[varnm_nc].isel(time=itime).data

dtopo = os.path.join(pthtopo,'ocean_topog.nc')
ds_topo = xarray.open_dataset(dtopo)
HH = -ds_topo['depth'].data

LON = ds_irlx['lon'].data
LAT = ds_irlx['lat'].data

match varnm:
  case('ithkn'):
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 5.
  case('iconc'):
    clrmp = mclrmps.colormap_conc()
    rmin = 0.
    rmax = 1.

clrmp.set_bad(color=[0.2, 0.2, 0.2])

# Mask Land:
A2d = np.where(HH>=0, np.nan, A2d)

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
img = ax1.pcolormesh(A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.contour(HH,[0],linestyles='solid', linewidths=1, colors=[(0., 0., 0.)])

if varnm == 'ithkn':
  ax1.contour(A2d,[1,2,3,4,5], linestyles='solid', colors=[(0.95, 0.95, 0.95)])
  
# Show test pnt:
if j0 >= 0 and i0 >=0:
  ax1.plot(i0,j0,'o')
ax1.axis('scaled')

sttl = (f'ARC12 irlx target field: {varnm} 1993/{MM}')
if j0 >= 0 and i0 >= 0:
  sttl = sttl + f"\n Test pnt iF0/jF0 = {iF0}/{jF0} {varnm}={A2d[j0,i0]:.6f}"

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

btx = 'plot_arc_rlx.py'
bottom_text(btx, pos=[0.2, 0.01])



