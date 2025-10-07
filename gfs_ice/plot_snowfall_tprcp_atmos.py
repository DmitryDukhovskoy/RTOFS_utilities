"""
  Derive snowfall rate from GFSv17 forecasts
  output from atmospheric model
  use tprcp and cpofp to get instant. snow fall

  From Ruiy Sun:
  In the GFSv17 surface history files there is a field called tprcp. 
  Tprcp is the instantaneous total precipitation (LS + convection). 
  There is another field called cpofp, which is the ratio : 
  SR(i,j) = (pptsnow + pptgraul + pptice) / (RAINNCV(i,j)+R1). 

  You can use the two to derive the time step snow + ice + graupel  
  at the location where no convective precipitation happens. 

  float cpofp(time, grid_yt, grid_xt) ;
    cpofp:_FillValue = 9.99e+20f ;
    cpofp:cell_methods = "time: point" ;
    cpofp:long_name = "Percent frozen precipitation" ;
    cpofp:missing_value = 9.99e+20f ;
    cpofp:output_file = "sfc" ;
    cpofp:units = "fraction" ;


  float tprcp(time, grid_yt, grid_xt) ;
    tprcp:_FillValue = 9.99e+20f ;
    tprcp:cell_methods = "time: point" ;
    tprcp:long_name = "total time-step precipitation" ;
    tprcp:missing_value = 9.99e+20f ;
    tprcp:output_file = "sfc" ;
    tprcp:units = "m" ;

  dt - from ufs.model_configure
  dt = 150 sec

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
varnm = 'tprcp'
#rho_snow = 300.
pltfld = 'mean'  # mean snowfall rate or cumulative fields to plot

parser = argparse.ArgumentParser()
parser.add_argument("--expt", help="expt name, e.g. rt13_upd01_stream3", type=str)
parser.add_argument("--init", help=f"init date, default={init_date}", type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--fhrS", help=f"forecast hour, Start avrg: 6,12,...,240", type=int, required=True)
parser.add_argument("--fhrE", help=f"forecast hour, End avrg: 6,12,...,240", type=int)
parser.add_argument("--pfld", help=f"Plot mean or accumulated field, default={pltfld}",
                    choices=['mean','cumul'],
                    type=str)
args = parser.parse_args()

expt = args.expt if args.expt else expt
init_date = args.init if args.init else init_date
init_hr = args.ihr if args.ihr else init_hr
fhrS = args.fhrS if args.fhrS else None
fhrE = args.fhrE if args.fhrE else fhrS
pltfld = args.pfld if args.pfld else pltfld
TLON = TLAT = LMSK = None

dlt_hr = 3  # delta hours between saved/avrg output 
dTsec = 150.   # atm. time step, preceip are instantaneous, "time-step precip"
HRFCST = np.arange(fhrS,fhrE+1,dlt_hr).astype(int)

pthoutp = f"/work/Dmitry.Dukhovskoy/GFSv17/{expt}/gfs.{init_date}/{init_hr:02d}/atmos"

varnc = varnm
irec = 0
AAsum = None
for hrf in HRFCST:
  flinp = f"gfs.t00z.sfcf{hrf:03d}.nc"
  dflice = os.path.join(pthoutp,flinp)

  print(f"Reading {varnm} from {dflice}")
  with xarray.open_dataset(dflice) as ds:
    A2d = ds[varnc].isel(time=0).squeeze().data
    PercFrz = ds['cpofp'].isel(time=0).squeeze().data
    #Rho_snow = ds['rhonewsn'].isel(time=0).squeeze()

    if TLON is None:
      TLON = ds['lon'].data
      TLAT = ds['lat'].data
      LMSK = ds['land'].data.squeeze()   # sea-land-ice mask (0-sea, 1-land, 2-ice)

  # Convert snow weight (kg/m2) --> cm of snow accumulated over a time step
  Fsnow = (A2d * PercFrz / dTsec) * 100.  # cm/sec

  if AAsum is None:
    #AIsum = Aice.copy()
    AAsum = Fsnow.copy()
  else:
    #AIsum = AIsum + Aice
    AAsum = AAsum + Fsnow

  irec += 1


# Convert AAsum --> cm/day - mean snowfall rate 
A2d = AAsum * 86400./float(irec)

if pltfld == 'cumul':
  A2d = AAsum.copy() * dTsec    # cumul snowfall in cm over all hourly output


#if irec > 1:
#  Aice = AIsum / float(irec)
#  A2d = AAsum / float(irec)

# Mask land:
#A2d[LMSK==1] = np.nan
jdim, idim = A2d.shape


if pltfld == 'mean':
  varplt = 'snowfall'
  units = 'cm/day'
  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  sinfo = 'Derived: (tprcp * cpofp)/dt_atmos *100 cm), tprcp - total time-step precip \n'
elif pltfld == 'cumul':
  varplt = 'cumul snowfall'
  units = 'cm'
  clrmp = mclrmps.colormap_warm()
  rmin = 0.
  rmax = 10.
  sinfo = 'Derived: (tprcp * cpofp)/dt_atmos *100 cm), tprcp - total time-step precip \n'

sinfo = sinfo + pthoutp

clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.6,0.6,0.6]

#sttl = f'{varnm}/ai {units}, GFSv17 {expt} init:{init_date}/{init_hr} fcast:{fhrS}-{fhrE}'
sttl = f'{varplt} {units}, GFSv17 {expt} init:{init_date}/{init_hr} fcast:{fhrS}-{fhrE}'

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
#CS = ax1.contour(xh, yh, Aice, [0.15], linestyles='solid', colors=[cntr_clr], linewidths=1)

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

btx = 'plot_snowfall_atmos.py'
#btx = 'plot_snow_ant.py'
bottom_text(btx, pos=[0.2, 0.01])


