"""
  Plot time series of positive/ negative snow components from GFSv17 forecasts
  POsitive: snowfall rate
  all others - considered "negative"

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
from mpl_toolkits.basemap import Basemap, cm
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
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob

expt = 'rt13_upd01_stream3'
init_date = 20250104
init_hr = 0

# hs_h - grid cell mean (!) snow thickness, m
# snow_ai - snowfall rate cm/day
# dsnow_h - snow formation (cm/day)
# snoice_h - snow-ice formation (cm/day)
# melts_h  - top snow melt (cm/day)
parser = argparse.ArgumentParser()
parser.add_argument("--expt", help="expt name, e.g. rt13_upd01_stream3", type=str)
parser.add_argument("--init", help=f"init date, default={init_date}", type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--fhrS", help=f"forecast hour, Start avrg: 6,12,...,240", type=int, required=True)
parser.add_argument("--fhrE", help=f"forecast hour, End avrg: 6,12,...,240", type=int)
args = parser.parse_args()

expt = args.expt if args.expt else expt
init_date = args.init if args.init else init_date
init_hr = args.ihr if args.ihr else init_hr
fhrS = args.fhrS if args.fhrS else None
fhrE = args.fhrE if args.fhrE else fhrS
TLON = TLAT = LMSK = None

dlt_hr = 6  # delta hours between saved/avrg output 
dlt_day = dlt_hr/24.
HRFCST = np.arange(fhrS,fhrE+1,dlt_hr).astype(int)

pthoutp = f"/work/Dmitry.Dukhovskoy/GFSv17/{expt}/gfs.{init_date}/{init_hr:02d}"



sinfo = pthoutp

irec = 0
VFall_sn = []
VMelt_sn = []
Vol_sn = []
for hrf in HRFCST:
  flinp = f"gfs.ice.t00z.6hr_avg.f{hrf:03d}.nc"
  dflice = os.path.join(pthoutp,flinp)

  print(f"Reading {dflice}")
  # hs_h - grid cell mean (!) snow thickness, m
  # snow_ai - snowfall rate cm/day
  # dsnow_h - snow thickn change (cm/day)
  # snoice_h - snow-ice formation (cm/day)
  # melts_h  - top snow melt (cm/day)
  with xarray.open_dataset(dflice) as ds:
    Aice = ds['aice_h'].isel(time=0).squeeze().data 
    fsnow = ds['snow_ai_h'].isel(time=0).squeeze().data*0.01   # cm/day --> m/day
    melts = ds['melts_h'].isel(time=0).squeeze().data*0.01     # cm/day --> m/day
    hsnow = ds['hs_h'].isel(time=0).squeeze().data             # m, grid cell mean
    if TLON is None:
      TLON = ds['TLON'].data
      TLAT = ds['TLAT'].data
      LMSK = ds['tmask'].data

      DX, DY = mmom6.dx_dy(TLON, TLAT)
      Acell = DX*DY
      Lantrc = np.where(TLAT>-50., 0, LMSK)

  # Saved snowfall rate is weighted by ice area to give mean
  # grid cell mean rate
  # Convert m3(snow)/m2(cell)*sec --> m3(snow)/m2(ice)*sec 
  #fsnow = np.divide(fsnow, Aice, out=np.zeros_like(fsnow), where=Aice > 0)
  
  # Volume of snowfall over Antarctica
  fsnow = np.where(Lantrc==0, np.nan, fsnow)
  vfs = np.nansum(fsnow * Acell * dlt_day)  # m3/day --> m3
  VFall_sn.append(vfs)

  # Volume of melted snow, m3
  melts = np.where(Lantrc==0, np.nan, melts)
  vmelt = np.nansum(melts * Aice * Acell * dlt_day)  # m3, ice area only
  VMelt_sn.append(vmelt)

  # Volume of snow on ice, m3:
  hsnow = np.where(Lantrc==0, np.nan, hsnow)  # grid cell mean snow thickness, m
  vsn   = np.nansum(hsnow * Acell)    # m3
  Vol_sn.append(vsn)

  irec += 1


VFall_sn = np.array(VFall_sn) * 1e-9   # m3 --> km3
VMelt_sn = np.array(VMelt_sn) * 1e-9   # m3 --> km3 
Vol_sn   = np.array(Vol_sn) * 1.e-9    # m3 --> km3

fday = HRFCST/24.
time_days = np.arange(np.ceil(max(fday))+1)

sttl = f'Vol snowfall & top snow melt, km3,  {expt} init:{init_date}/{init_hr} fcast:{fhrS}-{fhrE}'

plt.ion()

clrfsn = [0., 0.5, 0.9]
clrmlt = [0.9, 0.3, 0]
clrvsn = [0.8, 0., 0.6]

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])

ln1, = ax1.plot(fday, VFall_sn, '.-', linewidth=2, color=clrfsn, label='snowfall')
ln2, = ax1.plot(fday, VMelt_sn, '.-', linewidth=2, color=clrmlt, label='snow melt')

# Secondary axis (right y-axis), sharing the same x-axis
ax12 = ax1.twinx()
ln3, = ax12.plot(fday, Vol_sn, '.-', linewidth=2, color=clrvsn, label='vol snow')
ax12.set_ylabel('Snow Volume on ice, km3')
ax12.yaxis.label.set_color(clrvsn)
ax12.tick_params(axis='y', colors=clrvsn)

ax1.grid('on')
ax1.set_xticks(time_days)
ax1.set_xlim([0, time_days[-1]+0.2])
ax1.set_xlabel('F/cast days')
ax1.set_ylabel('Snow Volume, km3')
ax1.set_title(sttl)

# Legend
# Combine legends from both axes:
#lns = [ln1, ln2, ln3]
#labels = [l.get_label() for l in lns]
#ax1.legend(lns, labels, loc='upper left')

ax2 = plt.axes([0.8, 0.35, 0.1, 0.1])
lgd = plt.legend(handles=[ln1,ln2,ln3], loc='upper right')
ax2.axis('off')

ax3 = fig1.add_axes([0.1, 0.32, 0.8, 0.06])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')

#btx = 'plot_snowfall_ant.py'
btx = 'plot_snow_components.py'
bottom_text(btx, pos=[0.1, 0.28])


