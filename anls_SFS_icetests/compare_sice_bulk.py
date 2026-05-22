"""
  Compare ice bulk salinity
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib  
import xarray
import matplotlib.colors as colors 
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

init_date = 20240701
init_hr = 0    # nominal hr, actual: -6 hrs for IAU, and -3 FHROT (f/cast hr rotation)
regn = 'north'

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--init", help=f"init date", choices=[20240701, 20250101], default=init_date, type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--dend", help="End date to plot YYYYMMDD or provide --ndays", type=int)
parser.add_argument("--ndays", help=f"Optional: N days to show from init, will override dend", type=int)
parser.add_argument(
    "--enmb",
    help="List of experiment numbers (e.g., 1 3 9 12)",
    type=int,
    nargs="+",             # <-- allows one or more integers and will generate a list
    required=True
)
args = parser.parse_args()

regn      = args.regn if args.regn else regn
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
ENMBS     = args.enmb if args.enmb else None
end_date  = args.init if args.dend else None
ndays     = args.ndays if args.ndays else None

# Get date:
plot_init = False  # initial conditions

if end_date is None and ndays is None:
  raise RuntimeError("Both end_day and ndays are None, one of them has to be provided")

dnmbS = int(mtime.rdate2datenum(init_date*100))  # init. day nmb
YRS, MMS, DDS = mtime.datevec(dnmbS)[:3]

if ndays is not None:
  # ndays overrides end_date
  dnmbE = dnmbS + ndays
else:
  dnmbE = mtime.rdate2datenum(end_date*100)
YRE, MME, DDE = mtime.datevec(dnmbE)[:3]

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
    
fyaml = 'paths_sfs.yaml'
with open(fyaml) as ff:
  pths_sfs = safe_load(ff)
    
pthgrid = pths_sfs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
    
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

# Analysis region:
RMsk = (HH < 0) 
if regn == 'north':
  RMsk &= (hlat > 50)
elif regn == 'south':
  RMsk &= (hlat < -50)

# Daily mean output time is assumed except for the initial field
if plot_init:
  RECS = [dnmbS] + [x + 0.5 for x in range(dnmbS, dnmbE + 1)]
else:
  RECS = [x + 0.5 for x in range(dnmbS, dnmbE + 1)]

RECS   = np.array(RECS)
nrecs  = RECS.shape[0]
nexpts = len(ENMBS)
PRCT = np.zeros((nexpts, nrecs, 5))

iens = -1
for enmb in ENMBS:
  iens += 1
  irec = -1
  pthout_cice = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/sfs_C192mx025_cice_test/expt{enmb:02d}/cice6"

  for nn in range(nrecs):
    irec += 1
    dnmb = RECS[nn]
    YR, MM, DD, hr = mtime.datevec(dnmb, round_hrs=True)[:4]

    if hr == 0:
      # Initial state
      nsec0 = 0
      flcice = f"iceh_ic.{YR}-{MM:02d}-{DD:02d}-{nsec0:05d}.nc"
      dflcice = os.path.join(pthout_cice,flcice)

      if not os.path.isfile(dflcice) or not plot_init:
        print(f"Initial state file is missing or not requested plt_init, proceed without it ...")
        RMSE[irec, iens] = np.nan
        continue
    else:
      print(f"Processing {YR}/{MM}/{DD} expt{enmb:02d}...")
      flcice = f"iceh.{YR}-{MM:02d}-{DD:02d}.nc"
      dflcice = os.path.join(pthout_cice,flcice)

    print(f"Processing {dflcice}")
    with xarray.open_dataset(dflcice) as dcice:
      AI = dcice['aice_d'].values.squeeze()
      FRZ = dcice['sice_d'].values.squeeze()

    FRZ[AI < 0.01] = np.nan
    FRZ[~RMsk] = np.nan

    if np.all(np.isnan(FRZ)):
      frz10 = frz25 = frz50 = frz75 = frz90 = np.nan
    else:
      frz10 = np.nanpercentile(FRZ, 10)
      frz25 = np.nanpercentile(FRZ, 25)
      frz50 = np.nanpercentile(FRZ, 50)
      frz75 = np.nanpercentile(FRZ, 75)
      frz90 = np.nanpercentile(FRZ, 90)

    PRCT[iens,irec,:] = frz10, frz25, frz50, frz75, frz90


# Line colors:
import mod_gfs_cice_anls as mgfscice
importlib.reload(mgfscice)
CLRS = mgfscice.sens_tests_colors()

# Prepare stats for boxplots:
stats = []
positions = []
gap = 1.0         # space between time groups
width = 0.6       # spacing within group
for irec in range(nrecs):
  for iens in range(nexpts):
    p10, p25, p50, p75, p90 = PRCT[iens, irec, :]

    stats.append({
        'med': p50,
        'q1': p25,
        'q3': p75,
        'whislo': p10,
        'whishi': p90,
        'fliers': []
    })

    # centered grouping
    pos = irec * (nexpts + gap) + iens
    positions.append(pos)


colors = plt.cm.tab10(np.linspace(0, 1, nexpts))


plt.ion()

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.4, 0.85, 0.5])

box = ax1.bxp(
    stats,
    positions=positions,
    widths=width,
    patch_artist=True,
    showfliers=False
)

# Color ensembles
for iens, patch in enumerate(box['boxes']):
  iens = iens % nexpts
  patch.set_facecolor(colors[iens])
  patch.set_edgecolor('black')

# whiskers, medians styling
for median in box['medians']:
  median.set_color([1.,0.8,0])
  median.set_linewidth(1.5)

# X-axis: group centers
group_centers = [
    irec * (nexpts + gap) + (nexpts - 1) / 2
    for irec in range(nrecs)
]

ax1.grid(True, axis='y', alpha=0.3)
ax1.set_xticks(group_centers)
ax1.set_xticklabels([f"{irec}" for irec in range(nrecs)])

ax1.set_xlabel('Forecast days')
ax1.set_title(f'S ice bulk, SFS init {init_date}, regn={regn}')

# Legend
ax2 = plt.axes([0.08, 0.15, 0.5, 0.15])
handles = [
    plt.Line2D([0], [0], color=colors[i], lw=6)
    for i in range(nexpts)
]

leg_labels = []
for enmb in ENMBS:
  lbl = mgfscice.sfs_tests_info(enmb)
  leg_labels.append(f"{lbl}")

ax2.legend(handles, leg_labels, loc='upper left')
ax2.axis('off')

btx = 'compare_sice_bulk.py'
bottom_text(btx, pos=[0.1,0.1])

