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
import argparse
import pandas as pd

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

from mod_utils_fig import bottom_text
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_interp1D as mint1d
#importlib.reload(mutob)


parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="year to plot", type=int)
parser.add_argument("--mm", help="month to plot", type=int)
args = parser.parse_args()

if args.yr:
  YR = args.yr
if args.mm:
  MM = args.mm

pthoras = '/work/Dmitry.Dukhovskoy/data/ORAS5'
flnm = f'sossheig_control_monthly_highres_2D_{YR}{MM:02d}_CONS_v0.1.nc'
dfl = os.path.join(pthoras,flnm)

# Beaufort Sea region:
iBG1 = 158
iBG2 = 670
jBG1 = 800
jBG2 = 1020

print(f'Opening {dfl}')
varnm='sossheig'
dset = xarray.open_dataset(dfl)
ssh = dset[varnm].data.squeeze()


clrmp = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBuGn_r')
rmin = -0.5
rmax = 0.5

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

img = ax1.pcolormesh(ssh, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.axis('scaled')
ax1.set_xlim([iBG1,iBG2])
ax1.set_ylim([jBG1,jBG2])

sttl = f'SSH demeaned, ORAS5, {YR}/{MM:02d}'
ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)



