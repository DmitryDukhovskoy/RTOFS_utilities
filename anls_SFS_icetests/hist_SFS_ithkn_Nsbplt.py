"""
  Plot histograms of ithkn daily fields from
  SFS runs
  N subpots showing N days

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray as xr
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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6

init_date = 20250701
init_hr = 0
regn = 'north'

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--init", help=f"init date", choices=[20240701, 20250701], default=init_date, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fday", help=f"forecast day to plot: 0, 1, 2, ... , =0 - init. cond.", 
                    nargs="+",
                    type=int, 
                    required=True)
parser.add_argument("--enmb", help="experiment number: 0, 1, 2, ...", required=True, type=int)
parser.add_argument("--nrow", help=f"N of rows for subplots, default=N dates", type=int)
parser.add_argument("--ncol", help=f"N of columns for subplots, default=1", type=int)
args = parser.parse_args()

regn      = args.regn if args.regn else regn
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
FDAYS     = args.fday
enmb      = args.enmb
ncol      = args.ncol if args.ncol is not None else 1
nrow      = args.nrow if args.nrow is not None else len(args.fday)

if min(FDAYS) < 0:
  raise RuntimeError(f"fday should be > 0")
  

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

# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xr.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape


# Create mask of analyzed region:
Mask = HH < -20.
if regn == 'north':
  lat0 = 70
  Mask &= hlat > lat0
elif regn == 'south':
  lat0 = -60
  Mask &= hlat < lat0
else:
  raise ValueError(f"Unknown region: {regn}")

pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/sfs_C192mx025_cice_test/expt{enmb:02d}/cice6"

ITbins = np.arange(0,5,0.5)
ITbins[-1] = 100. 

def read_icefld(dflice, Mask):
  with xr.open_dataset(dflice) as dcice:
    for varnm in ['hi_h', 'hi_d']:
      if varnm in dcice:
        A2d = dcice[varnm].values.squeeze()
        break
    else:
      print(list(dcice.data_vars))
      raise KeyError(f"No aice_h or aice_d in {dflice}")

  A2d[~Mask] = np.nan

  return A2d

def plot_hist(ax1,hist_frac, bin_x0, edges, sttl):
  plt.sca(ax1)
  
  widths = 0.98 * np.diff(edges).copy()
  # Make the open-ended last bin visually the same width
  widths[-1] = widths[-2]
  edges[-1] = edges[-2] + 2*(bin_x0[-1] - edges[-2])

  plt.bar(bin_x0, hist_frac, width=widths, align='center', color=(0.2, 0.5, 0.8))

  # Ticks exactly at bin edges
  ax1.set_xticks(edges)

  ax1.set_xlim(-0.01, edges[-1]+0.01)
  ax1.set_ylim(0, 0.8)

  ax1.set_axisbelow(True)

  ax1.grid(
      True,
      linestyle='--',
      linewidth=0.7,
      color='gray',
      alpha=0.5
  )

  # Tick-label font sizes
  ax1.tick_params(axis='x', labelsize=12)
  ax1.tick_params(axis='y', labelsize=12)
  #ax1.set_xlabel('ice thickness, m')
  ax1.set_title(sttl)

  return


import mod_gfs_cice_anls as mgfscice
importlib.reload(mgfscice)
expt_name = mgfscice.sfs_tests_info(enmb)


print("Plotting ...")
plt.ion()
fig1 = plt.figure(1,figsize=(12, 9))
fig1.clf()  # Clear the figure
axes = fig1.subplots(nrows=nrow, ncols=ncol, squeeze=False)

fig1.subplots_adjust(
    left=0.08,
    right=0.95,
    top=0.97,
    bottom=0.11,
    wspace=0.1,
    hspace=0.2
)

for iplt, fday in enumerate(FDAYS):
  # Get date:
  plot_init = fday == 0  # initial conditions
  
  dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
  YRI, MMI, DDI, hrI = mtime.datevec(dnmbI, round_hrs=True)[:4]

  if plot_init:
    dnmb0 = dnmbI
  else:
    dnmb0 = dnmbI + fday-1                              # day to plot

  yr0, mm0, dd0, hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
  nsec0 = int(hr0 * 3600)
  YR, MM, DD = mtime.datevec(dnmb0)[:3]
  sttl = f"ithkn SFS expt{enmb:02d} FDAY={fday:02d} {YR}/{MM:02d}/{DD:02d}"

  irow = iplt // ncol
  icol = iplt % ncol
  ax1 = axes[irow, icol]

  if plot_init:
    flinp = f"iceh_ic.{yr0}-{mm0:02d}-{dd0:02d}-{nsec0:05d}.nc"
  else:
    flinp = f"iceh.{yr0}-{mm0:02d}-{dd0:02d}.nc"

  dflice = os.path.join(pthoutp,flinp)
  print(f"Processing {YR}/{MM}/{DD} {hr0:02d}:00, SFS init {YRI}/{MMI:02d}/{DDI:02d} {hrI:02d}:00\n{dflice}")

  A2d = read_icefld(dflice, Mask)
  hithkn = A2d[np.isfinite(A2d)]
  hist, edges = np.histogram(hithkn, bins=ITbins) 
  hist_frac = hist / hist.sum()

  bin_centers = 0.5 * (edges[:-1] + edges[1:])
  #Eliminate the infinity bin for plotting:
  bin_centers[-1] = edges[-2] + (edges[-2] - bin_centers[-2])
  
  plot_hist(ax1, hist_frac, bin_centers, edges, sttl)
 

sinfo = f'daily average ice thickness from CICE6\n'
sinfo = sinfo + f"SFS expt{enmb:02d}: {expt_name}) init {init_date}\n"
sinfo = sinfo + dflice

#Plot colrbar
# Force Matplotlib to finalize axes positions
fig1.canvas.draw()

# Get FINAL positions of the subplot axes
positions = [
    ax.get_position()
    for ax in axes.flat
    if ax.get_visible()
]

left   = min(p.x0 for p in positions)
right  = max(p.x1 for p in positions)
bottom = min(p.y0 for p in positions)

pbtm = bottom - 0.07


ax4 = fig1.add_axes([0.02, pbtm, 0.8, 0.06])
ax4.text(0, 0, sinfo, fontsize=8,
         ha='left', va='bottom',
         transform=ax4.transAxes)
ax4.axis('off')

pos_info = ax4.get_position()
bot_info = pos_info.y0
pbtm = bot_info - 0.02

btx = f' @{machine}: hist_SFS_ithkn_Nsbplt.py'
#bottom_text(btx, pos=[0.02,0.02], fsz=8)
bottom_text(btx, pos=[0.02, pbtm], fsz=8)



