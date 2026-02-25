"""
  Use snowfall rate from GFSv17 forecasts
  Output from CICE
  snowfall rate - averaged over output time period (6hr)

  in CICE, snowfall rate dumped to history is converted from kg/m2*s -->
  cm/day in liquid water equivalent  !!!

  in define_hist_field:
        call define_hist_field(n_snow_ai,"snow_ai","cm/day",tstr2D, tcstr, &
             "snowfall rate",                                             &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,           &
             ns1, f_snow_ai)

    rhofresh = 1000. 
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
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob

expt = 'rt13_upd01_stream3'
init_date = 20250104
init_hr = 0
val_mean = 'cell' # plot grid-cell mean values - given in CICE output, or ice area mean

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
parser.add_argument(
    "--mean", 
    help="cell mean (i.e. *ai) or ice area mean (/ai), default={val_mean}",
    choices=['cell','ice'],
    type=str
)

#parser.add_argument(
#    "--varnm",
#    help="Variable to plot",
#    choices=['snow_ai'],
#    required=True,
#    type=str
#)
args = parser.parse_args()

expt = args.expt if args.expt else expt
init_date = args.init if args.init else init_date
init_hr = args.ihr if args.ihr else init_hr
fhrS = args.fhrS if args.fhrS else None
fhrE = args.fhrE if args.fhrE else fhrS
val_mean = args.mean if args.mean else val_mean
#varnm = args.varnm if args.varnm else None
varnm = 'snow_ai'
rho_snow = 300.
TLON = TLAT = LMSK = None

dlt_hr = 6  # delta hours between saved/avrg output 
HRFCST = np.arange(fhrS,fhrE+1,dlt_hr).astype(int)

pthoutp = f"/work/Dmitry.Dukhovskoy/GFSv17/{expt}/gfs.{init_date}/{init_hr:02d}"


if varnm == 'dsnow_h':
  varnc = varnm
  units = 'cm/day'
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  sinfo = 'snow formation, cm/day\n'
elif varnm == 'snow_ai':
  varnc = 'snow_ai_h'
  units = 'cm/day'
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 1.
  if val_mean == 'ice':
    sinfo = 'snowfall rate, cm/day, ice area mean\n'
  else:
    sinfo = 'snowfall rate, cm/day, wegithed by ice area (grid cell mean)\n'

  sinfo = sinfo + f'Converted from snowfall in liquid equiv. water cm to snow cm, using rho_snow={rho_snow:.0f}\n'

sinfo = sinfo + pthoutp

irec = 0
AIsum = None
AAsum = None
for hrf in HRFCST:
  flinp = f"gfs.ice.t00z.6hr_avg.f{hrf:03d}.nc"
  dflice = os.path.join(pthoutp,flinp)

  print(f"Reading {varnm} from {dflice}")
  with xarray.open_dataset(dflice) as ds:
    Aice = ds['aice_h'].isel(time=0).squeeze().data 
    A2d = ds[varnc].isel(time=0).squeeze().data
    if TLON is None:
      TLON = ds['TLON'].data
      TLAT = ds['TLAT'].data
      LMSK = ds['tmask'].data

  # Saved snowfall rate is weighted by ice area to give mean
  # grid cell mean rate
  # Convert m3(snow)/m2(cell)*day --> m3(snow)/m2(ice)*day (??)
  if val_mean == 'ice':
    A2d = np.divide(A2d, Aice, out=np.zeros_like(A2d), where=Aice > 0)

  if varnm == 'snow_ai':
    # convert from liq. equiv. water --> snow depth
    A2d = A2d*1000./rho_snow

  if AIsum is None:
    AIsum = Aice.copy()
    AAsum = A2d.copy()
  else:
    AIsum = AIsum + Aice
    AAsum = AAsum + A2d

  irec += 1

if irec > 1:
  Aice = AIsum / float(irec)
  A2d = AAsum / float(irec)

A2d[LMSK==0] = np.nan
jdim, idim = A2d.shape


clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.6,0.6,0.6]

if val_mean == 'cell':
  sttl = f'{varnm} {units}, {val_mean}-mean, GFSv17 {expt} init:{init_date}/{init_hr} fcast:{fhrS}-{fhrE}'
else:
  sttl = f'{varnm}/ai {units}, {val_mean}-mean, GFSv17 {expt} init:{init_date}/{init_hr} fcast:{fhrS}-{fhrE}'

plt.ion()

m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
#lons, lats = m.makegrid(idim, jdim) # get lat/lons of ny by nx evenly spaced grid.
#x, y = m(lons, lats) # compute map proj coordinates.
xh, yh = m(TLON,TLAT) # GFS coords

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])
m.drawcoastlines()

# draw parallels.
parallels = np.arange(-80,-10,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,359.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
# Contour ice edge:
CS = ax1.contour(xh, yh, Aice, [0.15], linestyles='solid', colors=[cntr_clr], linewidths=1)

ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.03, 0.8, 0.06])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')

btx = 'plot_snowfall_ant.py'
#btx = 'plot_snow_ant.py'
bottom_text(btx, pos=[0.2, 0.01])


