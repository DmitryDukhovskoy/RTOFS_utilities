"""
  Plot accumulated snow melt
  Antarctica Region

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
varnm = 'melts_h'  
TLON = TLAT = LMSK = None

dlt_hr = 6  # delta hours between saved/avrg output 
dlt_day = float(dlt_hr)/24.
HRFCST = np.arange(fhrS,fhrE+1,dlt_hr).astype(int)

pthoutp = f"/work/Dmitry.Dukhovskoy/GFSv17/{expt}/gfs.{init_date}/{init_hr:02d}"


varnc = varnm
units = 'cm/day'
clrmp = mclrmps.colormap_haline2(start_clr=[1,1,1])
rmin = 0.
rmax = 10.
sinfo = 'cumul. snow melt, cm\n'

sinfo = sinfo + pthoutp

iedge = 0.15
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

  # Compute accumulated snow, cm
  # Need to convert grid cell mean cm/day to ice mean cm/day -- or not ???
  # Convert grid cell mean to ice mean, i.e. cm/(day * m2 of ice)
  # = A2d/Aice where Aice>0
  #A2d = np.divide(A2d, Aice, out=np.zeros_like(A2d), where=Aice > 0)
  

  if AIsum is None:
    AIsum = Aice.copy()
    AAsum = A2d * dlt_day  # cm/ (day * m2_ice) --> cm / m2_ice
  else:
    AIsum = AIsum + Aice
    AAsum = AAsum + A2d * dlt_day

  irec += 1

if irec > 1:
  Aice = AIsum / float(irec)

A2d  = AAsum.copy()
A2d[LMSK==0] = np.nan
jdim, idim = A2d.shape


clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.5,0.5,0.5]

sttl = f'Cumul.snow melt, cm, GFSv17 {expt} init:{init_date}/{init_hr} fcast:{fhrS}-{fhrE} hr'

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
CS = ax1.contour(xh, yh, Aice, [iedge], linestyles='solid', colors=[cntr_clr], linewidths=1)

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


btx = 'plot_accumsnow_ant.py'
bottom_text(btx, pos=[0.2, 0.01])


