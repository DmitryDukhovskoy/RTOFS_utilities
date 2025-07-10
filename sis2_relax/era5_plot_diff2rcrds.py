"""
  Compare 2 recrods from the atm. forinc file to investigate
  SIS2 blowups 
  Some ERA5 padded forcing files cause ice continuity 
  blow up (negative ice thickness)

  Happened on Dec. 31 2000, 24:00 (1st record in 2001 padded to the end of 2000)
  

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
parser.add_argument("--date1", help="Date 1: YYYYMMDDHH, or YYYYMMDD-->HH=0", type=int)
parser.add_argument("--date2", help="Date 2: YYYYMMDDHH, or YYYYMMDD-->HH=0", type=int)
parser.add_argument("--varnm", help="field to change: u10,v10,lp,msl,sf,sphum,ssrd,strd,t2m", type=str)
args = parser.parse_args()

#VARS = ['u10','v10','lp','msl','sf','sphum','ssrd','strd','t2m']
#VARS = ['msl','sf','sphum']
#VARS = ['ssrd','strd','t2m']
#VARS = ['u10','v10','lp']

f_save = True

dltHR = 1.
if args.varnm:
  varnm = args.varnm
if args.date1:
  date1 = int(args.date1)
  if date1 // 10**8 < 1:
    date1 *= 100
  yr1, mm1, dd1, hh1 = mtime.extract_yymmdd(date1)
  dnmb1 = mtime.datenum([yr1,mm1,dd1,hh1])
if args.date2:
  date2 = int(args.date2)
  if date2 // 10**8 < 1:
    date2 *= 100
  yr2, mm2, dd2, hh2 = mtime.extract_yymmdd(date2)
  dnmb2 = mtime.datenum([yr2,mm2,dd2,hh2])

# Assuming the plotted dates are in the same file
# file year is:
YR = np.min([yr1,yr2])

def read_ERA5_field(varnm,dnmbR):
  """
    Read ERA5 field for given date dnmbR
    Assumed hourly data
  """
  YRR,MMR,DDR,hrr,minr = mtime.datevec(dnmbR, round_hrs=True)

  pthera = '/archive/e1n/mom6/NEP/atmos_forcing/era5_padded'
  ptherafld = os.path.join(pthera,varnm)
  flera = f'ERA5_{varnm}_{YRR}_padded.nc'
  # Modified fields:
  #ptherafld = f'/work/Dmitry.Dukhovskoy/NEP_input/ERA5_padded_changed/{YRR}/'
  #flera = f'ERA5_{varnm}_{YRR}_cp049Bhrs_padded.nc'
  dflera = os.path.join(ptherafld,flera)

  print(f'Processing {YRR} {varnm}')
  print(f'Opening {dflera}')
  dset = xarray.open_dataset(dflera)
  Time = dset['time'].data
  tmP = pd.to_datetime(Time)
  nrec = len(tmP)
  TNEP = np.zeros((nrec,4), dtype=int)
  years  = tmP.year.to_numpy()
  months = tmP.month.to_numpy()
  days   = tmP.day.to_numpy()
  hours  = tmP.hour.to_numpy()
  TM = np.zeros((nrec))
  #TM     = mtime.datenum([years,months,days,hours]) <-- need to change mtime.datenum to work with 1D arrays
  for irec in range(nrec):
    yy,mm,dd,hh = years[irec],months[irec],days[irec],hours[irec]
    TM[irec] = mtime.datenum([yy,mm,dd,hh])

  # Check if requested date is in the time range:
  assert(dnmbR >= TM[0] and dnmbR <= TM[-1]), 'Requested date is outside the time range in the file'    
  DD = abs(TM-dnmbR)
  itime = np.argmin(DD)
  dnmb_found = TM[itime]
  yF,mF,dF,hF = mtime.datevec(dnmb_found, round_hrs=True)[:4]
  print(f'Found rec {itime}: {yF}/{mF}/{dF}:{hF:02d}') 
  assert(DD[itime] < 1.e-3),f'Could not find requested date  {dnmbR}'
  A2d = dset[varnm].isel(time=itime).data

  return A2d, TM, dset, dflera
  
# Read data:
print(f'Date: {date1}')
A1, TMp, _, dflera = read_ERA5_field(varnm, dnmb1)
print(f'Date: {date2}')
A2, _, _, _ = read_ERA5_field(varnm, dnmb2)
dltA = A1-A2

clrmp = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBu_r')
clrmp.set_bad(color=[0.2,0.2,0.2])
rmin = -10.
rmax = 10.

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
img = ax1.pcolormesh(dltA, cmap=clrmp, vmin=rmin, vmax=rmax)
#CS = ax1.contour(xR,yR,dltT,tscntrs, linestyles='solid', linewidths=1, colors=[(0., 0., 0.)])
#ax1.clabel(CS, tslabels,inline=1, fontsize=10)
sttl = f'Diff {varnm} {date1}-{date2} \n {dflera}' 
ax1.set_title(sttl)

# extend: min, max, both
ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')

ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


btx = 'era5_plot_diff2rcrds.py'
bottom_text(btx, fsz=6, pos=[0.08, 0.03])



