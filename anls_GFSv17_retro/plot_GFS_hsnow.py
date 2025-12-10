"""
  Plot snow depth fields from
  GFSv17 retro

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
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6

init_date = 20250604
init_hr = 6
regn = 'south'
fld_avrg = "cell"
runnm = 'retrov17_01' 
strnm = '4'
hrS = 0      # 1st forecast, hr, =0 - initial state
hrE = 384    # last forecast, hr
dltHR = 6    # output time freq, hrs


# hs_h - grid cell mean (!) snow thickness, m
# snow_ai - snowfall rate cm/day (in liquid water equivalent !)
# dsnow_h - snow formation (cm/day) - can be > or < 0
# snoice_h - snow-ice formation (cm/day)
# melts_h  - top snow melt (cm/day)
parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--avrg", help=f"plot grid cell or ice area mean: cell or ice, default{fld_avrg}", type=str)
parser.add_argument("--runnm", help=f"Name of the run, default={runnm}", type=str)
parser.add_argument("--strnm", help=f"Stream, default={strnm}", choices=['1a','4'], type=str)
parser.add_argument("--init", help=f"init date, default {init_date}", type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fhr", help=f"f/cast hr to plot: {hrS}:{hrE} default={hrE}", type=int)
args = parser.parse_args()

regn      = args.regn if args.regn else None
run_name  = args.runnm if args.runnm else runnm
str_name  = args.strnm if args.strnm else strnm
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fcst_hr   = args.fhr if args.fhr is not None else hrE
fld_avrg  = args.avrg if args.avrg else fld_avrg
  
assert hrS <= fcst_hr <= hrE, f"Requested f/cast hour {fcst_hr} is outside time range: {hrS}/{hrE}"

# Get date:
plot_init = fcst_hr == 0  # initial conditions
 
syst_info = os.uname() 
machine = syst_info.nodename
  
if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
else:
  print("Unknown machine:", machine)

run_full = f"{run_name}_stream{str_name}"
if init_date is not None and init_hr is not None:
  RUNS = [f"{init_date}{init_hr:02d}"]
else:
  RUNS = mgfscice.gfs_retro_runs(run_full, node_nm, model='ice')

fyaml = 'gfs17_paths.yaml'
with open(fyaml) as ff:
  pths_gfs = safe_load(ff)

# Get MOM6 grid
pthgrid = pths_gfs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

# Output to plot:
# Assumed all runs have same forecast duration
hrs_outp = np.array([x for x in range(hrS, hrE+1, dltHR)])
dhr = np.abs(hrs_outp - fcst_hr)
iplot = np.argmin(dhr)
fcst_hr = hrs_outp[iplot]  # correct f/cast hour if needed

dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb

# Date to plot:
dnmb0 = dnmbI + fcst_hr/24
yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
YR,MM,DD = mtime.datevec(dnmb0)[:3]
nsec0 = hr0*3600

comroot = pths_gfs[node_nm]["CICE6"]["comroot"].format(
  run_name=run_name,
  stream=str_name,
  init=init_date,
  ihr=init_hr
  )

pthout_cice = os.path.join(pths_gfs[node_nm]["CICE6"]["pthcice"].format(
  comroot=comroot
))
print(f"cice dir: {pthout_cice}")


if plot_init:
  flinp = f"gfs.t{init_hr:02d}z.ic.nc"
else:
  flinp = f"gfs.t{init_hr:02d}z.{dltHR}hr_avg.f{fcst_hr:03d}.nc"

dflice = os.path.join(pthout_cice,flinp)

print(f"Processing {YR}/{MM}/{DD}:{fcst_hr:02d} run {init_date}:{init_hr:02d}...")
print(f"Processing {dflice}")
with xarray.open_dataset(dflice) as dcice:
  A2d  = dcice['hs_h'].data.squeeze()
  Aice = dcice['aice_h'].squeeze().data
  LMSK = dcice['tmask'].data

if fld_avrg == 'ice':
  # Plot hsnow avrg over ice area:
  A2d = np.divide(A2d, Aice, out=np.zeros_like(A2d), where=Aice > 0)
A2d[LMSK==0] = np.nan
Aice[LMSK==0] = np.nan

units = 'm'
clrmp = mclrmps.colormap_temp()
rmin = 0.
rmax = 0.5


clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]

runname_date =f"{run_name}_stream{str_name} {init_date}{init_hr:02d}"

if fld_avrg == 'cell':
  sttl = f'{runname_date} hsnow m3/m2_cell, {YR}/{MM:02d}/{DD:02d}'
  sinfo = 'grid cell mean snow thickness, (m3 per m2 of grid cell)\n'
else:
  sttl = f'{runname_date} hsnow m3/m2_ice, {YR}/{MM:02d}/{DD:02d}'
  sinfo = 'mean snow thickness, (m3 per m2 of ice area)\n'
if plot_init:
  sttl = sttl + ' INIT'
sinfo = sinfo + dflice

plt.ion()


if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
  xl1 = -8.5e6
  xl2 = -1.e6
  yl1 = xl1
  yl2 = xl2

xh, yh = m(hlon, hlat) # CICE6 coordinates


print("Plotting ...")

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])
m.drawcoastlines()

# draw parallels.
if regn == 'south':
  parallels = np.arange(-80,-10,10.)
else:
  parallels = np.arange(40,89,10.)

m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,359.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

# Contour ice edge
CS = ax1.contour(xh, yh, Aice, [0.15], linestyles='solid', colors=[cntr_clr], linewidths=1)

if regn == 'south':
  ax1.set_xlim([xl1, xl2])
  ax1.set_ylim([yl1, yl2])
  ax1.invert_yaxis()
  ax1.invert_xaxis()

ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
if rmin < 0:
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
else:
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

btx = 'plot_GFS_hsnow.py'
bottom_text(btx, pos=[0.2, 0.01])


