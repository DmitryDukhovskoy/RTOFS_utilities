"""
  Plot monthly or N-month average PIOMAS fields
  interpoalted to NEP10k grid

"""
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
from yaml import safe_load
import argparse
from pathlib import Path
import sys
                      
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from MyPython.mod_time import datevec, datenum
from MyPython.mod_mom6 import read_mom6grid
from MyPython.mod_colormaps import colormap_ice_thkn, colormap_conc
from MyPython.mod_misc1 import convert_nptime_to_datenum
from MyPython.mod_utils_fig import bottom_text


parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="year to plot: 1993, ..., 2020", type=int, required=True)
parser.add_argument("--varnm", help="field to plot", 
                    choices=['ithkn','iarea','iconc'], type=str, required=True)
parser.add_argument(
    "--mo",
    help="List of months to plot N-month average, or 1 month (e.g., 4 5 6 - AMJ mean)",
    type=int,
    nargs="+",             # <-- allows one or more integers and will generate a list
    required=True
)
args = parser.parse_args()


YR     = args.yr
MONTHS = args.mo
varnm  = args.varnm
if varnm == 'iconc':
  varnm = 'iarea'


fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

# MOM6 NEP topo/grid:
run_name   = 'seasonal_fcst_daily'
pthtopo    = gridfls['MOM6_NEP'][run_name]['pthgrid']
fgrid      = gridfls['MOM6_NEP'][run_name]['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"][run_name]["ftopo"]
outdir     = gridfls['MOM6_NEP'][run_name]['pthoutp']
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
# Hgrid lon. lat:
hlon, hlat  = read_mom6grid(dfgrid_mom, grdpnt='hgrid')
HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

# PIOMAS fields on MOM6 NEP grid, relax fields:
pthsis  = gridfls['MOM6_NEP'][run_name]['pthsis']
#pthdata = '/work/Dmitry.Dukhovskoy/data/PIOMAS_ice'
#flthck  = f'piomas_heff{YR0}_v21.nc'
#varthck = 'heff'
#flconc  = f'piomas_area{YR0}_v21.nc'
#varconc = 'area'


# Read saved relax. fields:
YR1 = YR
YR2 = YR + 1
flout = f'PIOMASv21_ithkn_iconc_{YR1}_{YR2}_monthly.nc'
diclim = os.path.join(pthsis, flout)
print(f'Reading relax fields from {diclim}')
ds_rlx = xarray.open_dataset(diclim)
Time = ds_rlx['time'].data
TM = convert_nptime_to_datenum(Time)


A2d = np.zeros_like(HH, dtype=float)
icnt = 0
for MM in MONTHS:
  print(f"Reading {varnm} PIOMAS month={MM}")
  dnmb0 = datenum([YR,MM,15,12])
  D = abs(TM-dnmb0)
  itime = np.argmin(D)
  dv0 = datevec(TM[itime])
  assert dv0[0]==YR, f'Requested YR={YR}, year in rlx file={dv0[0]}'
  assert dv0[1]==MM, f'Requested month={MM}, month in rlx file={dv0[1]}'

  FLD = ds_rlx[varnm].isel(time=itime).values
  assert FLD.shape == A2d.shape
  A2d += FLD
  icnt += 1

A2d /= icnt
A2d[HH >= 0] = np.nan

match varnm:
  case('ithkn'):
    clrmp = colormap_ice_thkn()
    rmin = 0.
    rmax = 4.
  case('iarea'):
    clrmp = colormap_conc()
    rmin = 0.
    rmax = 1.

clrmp.set_bad(color=[0.2, 0.2, 0.2])


from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
          projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

if len(MONTHS) > 1:
  sttlS = f'PIOMASv2.1 {varnm} {YR} avrg:{MONTHS[0]}-{MONTHS[-1]} \n{flout}'
else:
  sttlS = f'PIOMASv2.1 {varnm} {YR}/{MONTHS[0]} \n{flout}'


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

img = ax1.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)

# Contour Thick ice
cntr_clr = [0.95, 0.95, 0.95]
hcntrs = [1,2,3,4,5]

CS = ax1.contour(xR, yR, A2d, hcntrs, linestyles='solid', colors=[cntr_clr], linewidths=1)
#ax1.clabel(CS, inline=1, fontsize=12)

ax1.set_title(sttlS)

ax3 = fig1.add_axes([0.2, 0.08, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ticks = np.linspace(rmin, rmax, 11)
clb.set_ticks(ticks)
clb.ax.set_xticklabels([f"{t:.1f}" for t in ticks], fontsize=12)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_piomas_ithkn_iconc_NEP10k.py'
bottom_text(btx, pos=[0.2, 0.03])




