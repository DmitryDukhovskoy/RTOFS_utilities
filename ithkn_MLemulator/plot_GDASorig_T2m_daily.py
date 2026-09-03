"""
  Plot daily mean  GDAS atm surface temperature
  Original fields
  For checking interpolation 

  Usa 6-hr fields downloaded from HPSS
  not processed
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray as xr
from yaml import safe_load
from mpl_toolkits.basemap import Basemap, cm
import argparse

# Append custom module paths
PPTHN = None
if 'PPTHN' not in locals() or PPTHN is None:
  cwd = os.getcwd()
  parts = cwd.split(os.sep)
  if 'python' in parts:
    idx = parts.index('python')
    PPTHN = os.sep + os.path.join(*parts[:idx + 1])
  else:
    raise RuntimeError("Directory 'python' not found in current working directory path.")

sys.path.extend([
    os.path.join(PPTHN, 'MyPython', 'hycom_utils'),
    os.path.join(PPTHN, 'MyPython', 'draw_map'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])


from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_sis2_relax as msisrlx
import mod_regmom as mrmom 

init_date = 20250701
init_hr = 0

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", 
                    choices=['north','south','global'],
                    type=str, required=True)
parser.add_argument("--rdate",  
             help="GDAS date YYYYMMDD to calc daily average and interpolate onto mesh025", 
             type=int,
             required=True)
parser.add_argument("--hr", help="use GDAS hour(s) for daily avrg., list, default=[0,6,12,18]",
                    type=int,
                    nargs="+")

args     = parser.parse_args()
regn     = args.regn if args.regn else None
gdas_date = args.rdate
HRS      = args.hr if args.hr is not None else [0,6,12,18]

def read_gdas(dflgdas, k2c=True, flip_north=True):
  assert os.path.isfile(dflgdas), f"Not found: {dflgdas}"
  with xr.open_dataset(dflgdas) as ds:
    A2d = ds["tmp2m"].isel(time=0).values

  if k2c:
    A2d = A2d -273.15   # K --> Celsius

  if flip_north:
    # Flip array to have N. at the top:
    A2d = np.flipud(A2d)

  return A2d


def get_gdas_coord(dflgdas, flip_north=True):
  assert os.path.isfile(dflgdas), f"Not found: {dflgdas}"
  with xr.open_dataset(dflgdas) as ds:
    LON = ds["lon"].values
    LAT = ds["lat"].values

  if flip_north:
    # Flip array to have N. at the top:
    LON = np.flipud(LON)
    LAT = np.flipud(LAT)

  return LON, LAT


def T2m_daily_mean(pthgdas, HRS):
  """
    Derive daily mean T2m from GDAS
    using HRS hours
  """
  T2m_mean = None
  for ik, hrz in enumerate(HRS):
    flname = f"gdas.t{hrz:02d}z.sfc.f000.nc"
    print(f"Reading {flname}")
    dflgdas = os.path.join(pthgdas, flname)
    T2m = read_gdas(dflgdas)

    if T2m_mean is None:
      T2m_mean = T2m.copy()
    else:
      T2m_mean += T2m

  T2m_mean /= len(HRS)

  return T2m_mean

dnmb0 = mtime.rdate2datenum(gdas_date)
YR, MM, DD = mtime.datevec(dnmb0)[:3]
# GDAS data
pthgdas = f"/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/GFSv17/gdas.{gdas_date}"

hrz = HRS[0]
flname = f"gdas.t{hrz:02d}z.sfc.f000.nc"
print(f"Reading {flname}")
dflgdas = os.path.join(pthgdas, flname)
LON, LAT = get_gdas_coord(dflgdas)

# Derive daily mean:
T2m_day = T2m_daily_mean(pthgdas, HRS)

print(f"Plotting {gdas_date}")

clrmp = mclrmps.colormap_temperature_coldwarm()
rmin = -10.
rmax = 20.

clrmp = mclrmps.colormap_difference_negpos()
rmin = -20.
rmax = 20.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

#AP[HH >= 0] = np.nan   # land
#AP[np.isnan(AP) & (HH < 0)] = -1.  # ocean

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.15, 0.8, 0.8])
    
if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)

xL, yL = m(LON, LAT)

m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
m.drawcoastlines()

img = m.pcolormesh(xL, yL, T2m_day, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(f"Air T2m, GDAS orig, {YR}/{MM:02d}/{DD:02d}")

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='both')
ticks = np.linspace(rmin, rmax, 11)
clb.set_ticks(ticks)
clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_GDASorig_T2m_daily.py'
bottom_text(btx, pos=[0.02,0.02], fsz=8)


