"""
  Plot monhtly clim of snow thickness in Antarctica 
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
parser.add_argument("--yrS", help="Start year to average, 1998,..., 2007", type=int, required=True)
parser.add_argument("--yrE", help="End year to average, default=yrS", type=int)
args = parser.parse_args()

yrS    = args.yrS if args.yrS else None
yrE    = args.yrE if args.yrE else yrS
moS    = 1
moE    = 12
nmnths = moE-moS+1
ncol   = 4
nrow   = 3

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

kyrs = 0
H3D = None
for YR in range(yrS,yrE+1):
  kyrs += 1
  pthdata=f'/work/Dmitry.Dukhovskoy/data/snow_nasa/{YR}'
  for MM in range (moS,moE+1):
    print(' ')
    print(f'Processing {YR}/{MM}')

    jd1, jd2 = jdays(YR,MM,MM)
    icc = 0
    ASUM = None
    for jday in range(jd1,jd2+1):
      flnm = f's{YR}{jday:03d}.hs'
      dflin = os.path.join(pthdata,flnm)
      #print(f'Reading  {dflin}')

      # Some files may be missing, skip those
      if not os.path.isfile(dflin):
        print(f"File not found: {dflin}, skipping ...")
        continue  

      hs2D = read_hsnow(dflin).astype(float)
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

    print(f'MM={MM}, Min/max snow (cm) = {np.nanmin(HS2D):.4f}/{np.nanmax(HS2D):.4f}')

    if H3D is None:
      jdim, idim = HS2D.shape
      H3D = np.zeros((nmnths,jdim,idim))

    # Sum of N years:
    H3D[MM-moS,:,:] = H3D[MM-1,:,:]+HS2D


if kyrs > 1:
  H3D = H3D/float(kyrs)


plt.ion()
#plt.close('all')
#fig1, axes = plt.subplots(nrow, ncol, figsize=(15, 10))  # 3 rows, 4 columns
fig1 = plt.figure(1,figsize=(12, 10))
fig1.clf()  # Clear the figure
axes = fig1.subplots(nrows=nrow, ncols=ncol)

# Shift subplots to the left
# and up for colorbar and text
# also keep subplots close to each other : wspace, hspace
fig1.subplots_adjust(
    left=0.05,
    right=0.85,  # More right-side room
    top=0.95,
    bottom=0.1,
    wspace=0.05,
    hspace=0.1
)


units = 'cm'
clrmp = mclrmps.colormap_temp()
rmin = 0.
#rmax = 20.
rmax = 50.
#clrmp = mclrmps.colormap_ice_thkn()
#rmin = 0.
#rmax = 50.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

iplt = 0

for MM in range(moS,moE+1):
  sttl=f'hsnow MM={MM:02d}, {yrS}-{yrE}'

  irow = iplt // ncol
  icol = iplt % ncol
  iplt += 1

  ax1 = axes[irow, icol]

  hs2d = H3D[MM-moS,:,:].squeeze()
  img = ax1.pcolormesh(hs2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.axis('scaled')
  ax1.set_xlim([10,310])
  ax1.set_ylim([20,330])
  ax1.tick_params(
      axis='both',     # 'x', 'y', or 'both'
      which='both',    # 'major', 'minor', or 'both'
      bottom=False,    # ticks on bottom axis
      top=False,       # ticks on top axis
      left=False,      # ticks on left axis
      right=False,     # ticks on right axis
      labelbottom=False,  # remove x-axis tick labels
      labelleft=False     # remove y-axis tick labels
  )
  ax1.set_title(sttl)

# Colorbar
# extend: min, max, both
ax2 = fig1.add_axes([0.9,0.12,0.013,0.8])
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')

ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

sinfo = 'NASA, Snow thickness on sea ice, monthly clim, from SSM/I satellite obs\n'
sinfo = sinfo + 'Calculations based on 19 and 37 GHz microwave brightness temperatures\n'
sinfo = sinfo + ' (vertical polarization) as well as sea ice concentration measurements.'
ax3 = plt.axes([0.05, 0.04, 0.9,0.07])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')

btx='plot_hsnow_Nmnths_AMSRobs_Antarctic.py'
bottom_text(btx, pos=[0.1, 0.02])








