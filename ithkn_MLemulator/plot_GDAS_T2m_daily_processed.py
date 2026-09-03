"""
  Plot daily mean  GDAS atm surface temperature
  For checking interpolation 

  Use processed (daily averaged) fields downloaded from HPSS
  see: 
  For deriving several daily means, with fetching GDAS from HPSS, prior to interpolation:
  use sbatch derive_dailyT2m_Ndays.sh --sdate 20250701 --edate 20250731 --dt 3
  This will save daily felds at 3 day intervals for 2025/07

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


parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", 
                    choices=['north','south','global'],
                    type=str, required=True)
parser.add_argument("--rdate",  
             help="GDAS date YYYYMMDD to calc daily average and interpolate onto mesh025", 
             type=int,
             required=True)

args     = parser.parse_args()
regn     = args.regn if args.regn else None
gdas_date = args.rdate

fyaml = 'paths_ML.yaml'
with open(fyaml) as ff:
  pths_ml = safe_load(ff)


pthdata = pths_ml["GDAS"]["pthdata"]  # root dir for processed data
pthdaily = pths_ml["GDAS"]["pthdaily"]  # daily GDAS fields

# Get gmapi 4 GDAS grid points for interpolation
pthdata = pths_ml["GDAS"]["pthdata"]
pthdump = os.path.join(pthdata,'gmapi_mapping2mesh025')
fgmapi = f'GDASatm_reggrid_to_mesh025_gmapi_1440x1080_{regn}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)
print(f'Loading gmapi --> {dfgmapi}')

# Get GDAS coord, flipped N
with xr.open_dataset(dfgmapi) as dgmapi:
  LON = dgmapi['longit'].values
  LAT = dgmapi['latit'].values


dnmb0 = mtime.rdate2datenum(gdas_date)
YR, MM, DD = mtime.datevec(dnmb0)[:3]
# GDAS data

# Read daily mean:
dflt2m = os.path.join(pthdaily, f"GDAS_flipN_T2m_{YR}{MM:02d}{DD:02d}.npy")
assert os.path.isfile(dflt2m), f"File missing: {dflt2m}"

print(f"Loading T2m: {dflt2m}")
T2m_day = np.load(dflt2m)

print(f"Plotting {gdas_date}")

clrmp = mclrmps.colormap_temperature_coldwarm()
rmin = -10.
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

btx = 'plot_GDAS_T2m_daily_processed.py'
bottom_text(btx, pos=[0.02,0.02], fsz=8)


