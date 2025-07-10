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


pthdata = '/work/Dmitry.Dukhovskoy/ARC12/irlx'
pthtopo = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'
frlx = 'damping_full_t_30.nc'
dfrlx = os.path.join(pthdata,frlx)
ds_irlx = xarray.open_dataset(dfrlx)
RLXIS = ds_irlx['Idamp'].data

dtopo = os.path.join(pthtopo,'ocean_topog.nc')
ds_topo = xarray.open_dataset(dtopo)
HH = -ds_topo['depth'].data


# For checking, relaxation time, hrs:
RLXHR = RLXIS.copy()
RLXHR = np.where(RLXHR==0., np.nan, RLXHR)
RLXHR = 1./RLXHR * 1/3600.
max_rlx = np.nanmax(RLXHR)

cff = 1.e7
AP = RLXIS.copy()*cff
AP = np.where(HH>=0., np.nan, AP)

clrmp = mclrmps.colormap_temp2()
clrmp = mclrmps.colormap_conc()
clrmp.set_bad(color=[0.2, 0.2, 0.2])
rmin=0
rmax=3.5


fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
img = ax1.pcolormesh(AP, cmap=clrmp, vmin=rmin, vmax=rmax)
#ax1.contour(HH,[0],linestyles='solid', linewidths=1, colors=[(0., 0., 0.)])
ax1.axis('scaled')
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



