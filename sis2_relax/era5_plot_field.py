"""
  Plot fields from ERA5 atm. forcing files

  ERA5 atm. forcing files:
  ERA5_u10_2000_padded.nc
  ERA5_v10_2000_padded.nc
  ERA5_lp_2000_padded.nc
  ERA5_msl_2000_padded.nc
  ERA5_sf_2000_padded.nc
  ERA5_sphum_2000_padded.nc
  ERA5_ssrd_2000_padded.nc
  ERA5_strd_2000_padded.nc
  ERA5_t2m_2000_padded.nc
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
parser.add_argument("--yrf", help="FileYear, for padded data where end record may be next year", type=int)
parser.add_argument("--date", help="Date & time to plot: YYYYMMDDHH, or YYYYMMDD-->HH=0", type=int)
parser.add_argument("--varnm", help="field to plot: u10,v10,lp,msl,sf,sphum,ssrd,strd,t2m", type=str)
args = parser.parse_args()

YRF = 0
if args.varnm:
  varnm = args.varnm
if args.date:
  date = int(args.date)
  if date // 10**8 < 1:
    date *= 100
  yr, mm, dd, hh = mtime.extract_yymmdd(date)
  dnmbR = mtime.datenum([yr,mm,dd,hh])
if args.yrf:
  YRF=args.yrf

YRR,MMR,DDR,hrr,minr = mtime.datevec(dnmbR, round_hrs=True)
if YRF == 0:
  YRF = YRR

pthera = '/archive/e1n/mom6/NEP/atmos_forcing/era5_padded'
ptherafld = os.path.join(pthera,varnm)
flera = f'ERA5_{varnm}_{YRF}_padded.nc'
dflera = os.path.join(ptherafld,flera)
# Corrected fields:
pthera=f'/work/Dmitry.Dukhovskoy/NEP_input/ERA5_padded_changed/{YRF}/'
flera = f'ERA5_{varnm}_{YRF}_corrected_padded.nc'
dflera = os.path.join(pthera,flera)

print(f'Processing {YRF} {varnm}, requested time={date}')
print(f'Opening {dflera}')
dset = xarray.open_dataset(dflera, mask_and_scale=False)  # keep missing values unmasked
#dset = xarray.open_dataset(dflera)
Time = dset['time'].data
tmP = pd.to_datetime(Time)
nrec = len(tmP)
TNEP = np.zeros((nrec,4), dtype=int)
years  = tmP.year.to_numpy()
months = tmP.month.to_numpy()
days   = tmP.day.to_numpy()
hours  = tmP.hour.to_numpy()
TM = np.zeros((nrec))
#TM= mtime.datenum([years,months,days,hours]) <-- need to change mtime.datenum to work with 1D arrays
for irec in range(nrec):
  yy,mm,dd,hh = years[irec],months[irec],days[irec],hours[irec]
  TM[irec] = mtime.datenum([yy,mm,dd,hh])

# Check if requested date is in the time range:
assert(dnmbR >= TM[0] and dnmbR <= TM[-1]), 'Requested date is outside the time range in the file'    
DD = abs(TM-dnmbR)
itime = np.argmin(DD)
assert(DD[itime] < 1.e-3),f'Could not find requested date  {dnmbR}'
A2d = dset[varnm].isel(time=itime).data

print(f'Min/max values of {varnm} = {np.min(A2d)}/{np.max(A2d)}')
  
import mod_colormaps as mclrmp
clrmp = mclrmps.colormap_temp2()
#clrmp = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBu_r')
clrmp.set_bad(color=[0.2,0.2,0.2])
rmin, rmax = mclrmp.minmax_clrmap(A2d, cpnt=0.0000001)
rmin = 0.
rmax = 0.002
clrmp.set_under((0,0,0))

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
img = ax1.pcolormesh(A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
#CS = ax1.contour(xR,yR,dltT,tscntrs, linestyles='solid', linewidths=1, colors=[(0., 0., 0.)])
#ax1.clabel(CS, tslabels,inline=1, fontsize=10)
sttl = f'ERA5 {varnm} {date} \n {dflera}'
ax1.set_title(sttl)

# extend: min, max, both
ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')

ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.3f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


btx = 'era5_plot_diff2rcrds.py'
bottom_text(btx, fsz=6, pos=[0.08, 0.03])








