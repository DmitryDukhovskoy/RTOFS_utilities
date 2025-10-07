"""
  Plot snow thickness in Antarctica 
  derived from AMSR 19 and 37 GHz microwave brightness temperatures
  https://earth.gsfc.nasa.gov/cryo/data/antarctic-snow-depth-sea-ice
  available data: 1993-2008
  daily data
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

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps

parser = argparse.ArgumentParser()
parser.add_argument("--yrS", help="Start year to average, 1993,..., 2007", type=int, required=True)
parser.add_argument("--yrE", help="End year to average, default=yrS", type=int)
parser.add_argument("--moS", help="month to start averaging", type=int, required=True)
parser.add_argument("--moE", help="month to end averaging, default = moS", type=int)
args = parser.parse_args()

yrS    = args.yrS if args.yrS else None
yrE    = args.yrE if args.yrE else yrS
moS    = args.moS if args.moS else None
moE    = args.moE if args.moE else moS

assert moS <= moE, f'ERR: end month should be same or later than {moS}'

# Find year days for specified months:
def jdays(YR,moS,moE):
  jd1 = int(mtime.date2jday([YR,moS,1]))
  mo_days = mtime.month_days(moE,YR)
  jd2 = int(mtime.date2jday([YR,moE,mo_days]))

  return jd1, jd2

def read_hsnow(dflin):
  # Grid dimensions
  idim, jdim = 316, 332

  # Read the binary data
  with open(dflin, "rb") as fid:
    data = np.fromfile(fid, dtype=np.uint8)

  # Check size
  if data.size != jdim * idim:
    raise ValueError(f"Unexpected file size: expected {jdim * idim}, got {data.size}")

  A2D = data.reshape((jdim,idim))

  return A2D  

icc = 0
ASUM = None
for YR in range(yrS,yrE+1):
  pthdata=f'/work/Dmitry.Dukhovskoy/data/snow_nasa/{YR}'
  jd1, jd2 = jdays(YR,moS,moE)

  for jday in range(jd1,jd2+1):
    flnm = f's{YR}{jday:03d}.hs'
    dflin = os.path.join(pthdata,flnm)
    print(f'Reading  {dflin}')

    hs2D = read_hsnow(dflin)
    hs2D = hs2D.astype(float)
    # Flip the data to have correct orientation of Antarctica
    hs2D = np.flipud(hs2D)
    if ASUM is None:
      ASUM = hs2D.copy()
    else:
      ASUM = ASUM + hs2D
    icc += 1

#ASUM = ASUM.astype(float)
if icc > 1:
  HS2D = ASUM / float(icc)
else:
  HS2D = ASUM.copy()

# Land = 200
HS2D[HS2D>190.]=np.nan   # thickns in cm !, 200 - land

sttl = f'S.Ocean snow (cm), from SSM/I, avrg. {yrS}-{yrE}, {moS}-{moE}'


clrmp = mclrmps.colormap_ice_thkn()
clrmp.set_bad(color=[0.2, 0.2, 0.2])
rmin = 0.
rmax = 60.


plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

ax1.set_xlim([25,300])
ax1.set_ylim([10,280])
img = ax1.pcolormesh(HS2D, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(sttl)

# extend: min, max, both
ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')

ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_hsnow_AMSRobs_Antarctic.py'
bottom_text(btx, pos=[0.2, 0.01])








