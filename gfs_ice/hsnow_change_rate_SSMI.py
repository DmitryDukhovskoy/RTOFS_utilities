"""
  Calc snow vol change rate from 
  atm.-forced UFS (datm UFS) vs NASA SSM/I daily clim

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
import mod_gfs_cice_anls as mgfscice
importlib.reload(mgfscice)

expt = 'ufs_datm_mx025_v02'
init_date = 20250103
init_hr = 0
regn = 'south'
fdays = 14

# hs_h - grid cell mean (!) snow thickness, m
# snow_ai - snowfall rate cm/day (in liquid water equivalent !)
# dsnow_h - snow formation (cm/day) - can be > or < 0
# snoice_h - snow-ice formation (cm/day)
# melts_h  - top snow melt (cm/day)
parser = argparse.ArgumentParser()
parser.add_argument("--init", help=f"init date, default={init_date}", type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--fdays",help=f"Number of forecast days, default={fdays}", type=int)
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument(
    "--enmb",
    help="List of experiment numbers (e.g., 1 3 9 12)",
    type=int,
    nargs="+",             # <-- allows one or more integers and will generate a list
    required=True
)
args = parser.parse_args()

init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
regn      = args.regn if args.regn else regn
fdays     = args.fdays if args.fdays is not None else fdays
ENMBS = args.enmb if args.enmb else None

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

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)


# Get date:
dnmbS = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
#dnmb0 = dnmbS + fday-1                              # day to plot

#yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
YRS, MMS, DDS, hr0 = mtime.datevec(dnmbS, round_hrs=True)[:4]
nsec0 = hr0*3600
dnmbE = dnmbS + fdays-1 
YRE, MME, DDE  = mtime.datevec(dnmbE, round_hrs=True)[:3]


# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdim, idim = HH.shape

if regn == 'south':
  RMsk = np.where(HH>=0, 0, 1)
  RMsk = np.where(hlat > -60., 0, RMsk)

DX, DY = mmom6.dx_dy(hlon,hlat)
Acell = DX*DY

# Create an array of day numbers with 0hr = init cond, 12 hr - daily means
RECS = [int(dnmbS)] + [x + 0.5 for x in range(int(dnmbS), int(dnmbE) + 1)]
RECS = np.array(RECS)
RECScl = RECS.copy()
RECScl[0] = RECScl[0]-0.5

nexpts = len(ENMBS)
nrecs  = RECS.shape[0]
SVOL = np.zeros((nrecs,nexpts))
iens = -1
for enmb in ENMBS:
  iens += 1
  irec = -1
  pthout_cice = pths_ufs[node_nm]["MOM6"]["pthcice"].format(enmb=enmb)

  for nn in range(nrecs):
    dnmb = RECS[nn]
    YR,MM,DD,hr = mtime.datevec(dnmb, round_hrs=True)[:4]

    if hr == 0:
      # Initial state
      nsec0 = 0
      flcice = f"iceh_ic.{YR}-{MM:02d}-{DD:02d}-{nsec0:05d}.nc"
      dflcice = os.path.join(pthout_cice,flcice)

      if not os.path.isfile(dflcice):
        print(f"Initial state file is missing, proceed without it ...")
        irec += 1
        SVOL[irec, iens] = np.nan
        continue
    else:
      print(f"Processing {YR}/{MM}/{DD} expt{enmb:02d}...")
      flcice = f"iceh.{YR}-{MM:02d}-{DD:02d}.nc"
      dflcice = os.path.join(pthout_cice,flcice)

    print(f"Processing {dflcice}")
    with xarray.open_dataset(dflcice) as dcice:
      aice  = dcice['aice_d'].data.squeeze()
      hsnow = dcice['hs_d'].data.squeeze()     # grid cell-mean snow thickness, m

    snow_vol = hsnow * Acell * RMsk * 1e-9  # km3
    SVOL[nn, iens] = np.nansum(snow_vol)

# Read daily clim snow depth:
pthdata = pths_ufs[node_nm]['MOM6']['pthdata']
pthsnow = os.path.join(pthdata,'snow_nasa/daily_clim')

SVOL_clim = np.zeros((nrecs))
for nn in range(nrecs):
  dnmb = RECScl[nn]
  YR,MM,DD = mtime.datevec(dnmb, round_hrs=True)[:3]
  flsnow = f"SSMI_hsnow_mesh025_1440x1080_dailyclim_{MM:02d}_south.nc"
  dflsnow = os.path.join(pthsnow, flsnow)
  print(f"Reading hsnow clim {dflsnow}")
  with xarray.open_dataset(dflsnow) as ds_snow:
    hsnow_clim = ds_snow['snow_depth'].isel(time=DD-1).data.squeeze()

  snow_vol = hsnow_clim * Acell * RMsk * 1e-9  # km3
  SVOL_clim[nn] = np.nansum(snow_vol)

dlt_snow  = np.diff(SVOL, axis=0)
dlt_snow[0] = dlt_snow[0] * 2.     # half-day change
dlt_sclim = np.diff(SVOL_clim) 

# Line colors:
CLRS = mgfscice.sens_tests_colors()
clr_clim = [0.5, 0.5, 0.5]

XT  = RECS - np.floor(RECS[0])
XTC = RECScl - np.floor(RECS[0])
xticks = np.arange(np.floor(XT[0]),np.ceil(XT[-1]))
sttl = f"Snow vol SSM/I daily clim and datmUFS expts, region={regn}, \n"+\
       f"{YR}/{MMS:02d}/{DDS:02d}-{YR}/{MME:02d}/{DDE:02d}"

plt.ion()

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.55, 0.7, 0.4])
LNS = []
ln1, = ax1.plot(XTC,SVOL_clim, 'o-', linewidth=2, color=clr_clim, label='SSMI clim')
LNS.append(ln1)

for iens in range(nexpts):
  enmb = ENMBS[iens]
  svol0 = SVOL[:,iens]
  clr0  = CLRS[iens,:]
  lbl_expt = mgfscice.sens_tests_info(enmb)
  line_lbl  = f"{enmb:02d}: {lbl_expt}"
  ln1, = ax1.plot(XT,svol0, 'o-', linewidth=2, color=clr0, label=line_lbl)
  LNS.append(ln1)

yl1 = 0
yl2 = np.nanmax([np.nanmax(SVOL_clim),np.nanmax(SVOL)]) * 1.05
ax1.set_xticks(xticks)
ax1.set_ylim(yl1, yl2)
ax1.grid('on')
#ax1.set_xlabel('Forecast days')
ax1.set_ylabel('Snow Vol, km3')
ax1.set_title(sttl)

# Plot snow vol change rate km3/day
ax2 = plt.axes([0.08, 0.1, 0.7, 0.4])
ax2.plot(XTC[:-1], dlt_sclim, 'o-', linewidth=2, color=clr_clim)
for iens in range(nexpts):
  enmb = ENMBS[iens]
  dsvol = dlt_snow[:,iens]
  clr0  = CLRS[iens,:]
  ax2.plot(XT[:-1], dsvol, 'o-', linewidth=2, color=clr0)

ax2.set_xticks(xticks)
ax2.grid('on')
ax2.set_xlabel('Forecast days')
ax2.set_ylabel('km3/day')
ax2.set_title('d/dt snow vol')

ax3 = plt.axes([0.782, 0.1, 0.21, 0.4])
lgd = plt.legend(handles=LNS, loc='upper left', fontsize=8)
ax3.axis('off')


btx = 'hsnow_change_rate_SSMI.py'
bottom_text(btx, pos=[0.1, 0.02])


