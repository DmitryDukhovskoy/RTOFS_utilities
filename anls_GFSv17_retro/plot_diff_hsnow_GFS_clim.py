"""
  Plot difference of snow depth fields from sensitivity tests with
  atm.-forced UFS (datm UFS) vs NASA SSM/I daily clim
  interpolated to mesh025 

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
parser.add_argument("--punit", help="plot units default=m", choices=['cm','m'], type=str)
args = parser.parse_args()

regn      = args.regn if args.regn else None
run_name  = args.runnm if args.runnm else runnm
str_name  = args.strnm if args.strnm else strnm
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fcst_hr   = args.fhr if args.fhr is not None else hrE
fld_avrg  = args.avrg if args.avrg else fld_avrg
plot_units = args.punit if args.punit else 'm'

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


def units_conversion(units, plot_units):
  if units == 'm':
    if plot_units == 'm':
      cff = 1.
    else:
      cff = 0.01
  elif units == 'cm':
    if plot_units == 'm':
      cff = 100.
    else:
      cff = 1.
  else:
    raise Exception(f"input units are not recognized: {units}")

  return cff

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
YR,MM,DD = mtime.datevec(dnmb0)[:3]

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
  units = dcice['hs_h'].attrs.get('units')
  cff_units = units_conversion(units, plot_units)
  A2d  = dcice['hs_h'].data.squeeze()*cff_units
  Aice = dcice['aice_h'].squeeze().data
  LMSK = dcice['tmask'].data

if fld_avrg == 'ice':
  # Plot hsnow avrg over ice area:
  A2d = np.divide(A2d, Aice, out=np.zeros_like(A2d), where=Aice > 0)
A2d[LMSK==0] = np.nan
Aice[LMSK==0] = np.nan

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

# Read interpolated snow depths:
# Snow depth climatology, Interpolated fields mesh025:
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthsnow = os.path.join(pthdata,'snow_nasa','daily_clim')
flhsn = f"SSMI_hsnow_mesh025_1440x1080_dailyclim_{MM:02d}_south.nc" 
dflhsn = os.path.join(pthsnow,flhsn)
print(f"Reading interpolated hsnow {dflhsn}")
with xarray.open_dataset(dflhsn) as ds_snow:
  units = ds_snow['snow_depth'].attrs.get('units')
  cff_units = units_conversion(units, plot_units)
  HSi = cff_units* ds_snow['snow_depth'].isel(time=DD-1).data.squeeze()
  LON = ds_snow['lon'].data
  LAT = ds_snow['lat'].data

dHS = A2d - HSi

units = plot_units
clrmp = mclrmps.colormap_uv()
if plot_units == 'cm':
  rmin = -50
  rmax = 50
else:
  rmin = -0.5
  rmax = 0.5

clrmp.set_bad(color=[0.2, 0.2, 0.2])

runname_date =f"{run_name}_stream{str_name} {init_date}{init_hr:02d}"
if fld_avrg == 'cell':
  sttl = f'diff hsnow m3/m2_cell, GFS vs SSM/I daily clim\n {runname_date} hsnow m3/m2_cell, {YR}/{MM:02d}/{DD:02d}'
  sinfo = f'difference GFS-climatology grid cell mean snow thickness, {cff_units}*(m3 per m2 of grid cell\n'
else:
  sttl = f'diff hsnow m3/m2_ice, GFS vs SSM/I daily clim\n {runname_date} hsnow m3/m2_ice, {YR}/{MM:02d}/{DD:02d}'
  sinfo = f'difference GFS-climatology snow thickness, {cff_units}*(m3 per m2 of ice\n'
if plot_init:
  sttl = sttl + ' INIT'
sinfo = sinfo + dflice + '\n' + dflhsn


plt.ion()

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
  xl1 = -8.5e6
  xl2 = -1.e6
  yl1 = xl1
  yl2 = xl2

xh, yh = m(hlon, hlat) # CICE6 coordinates

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

img = ax1.pcolormesh(xh, yh, dHS, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

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

btx = 'plot_diff_hsnow_GFS_clim.py'
bottom_text(btx, pos=[0.2, 0.01])


