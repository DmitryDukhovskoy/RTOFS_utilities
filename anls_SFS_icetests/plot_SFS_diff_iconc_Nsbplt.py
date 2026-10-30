"""
  Plot difference of iconc daily fields from
  SFS runs and reference field (GLORYS)
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
parser.add_argument("--icref", help="Iconc reference fields, default=nsidc",
                   choices=['nsidc','amsr2'],
                   default='nsidc',
                   type=str)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fday", help=f"forecast day to plot: 0, 1, 2, ... , =0 - init. cond.", 
                    nargs="+",
                    type=int, 
                    required=True)
parser.add_argument("--enmb", help="experiment number: 0, 1, 2, ...", required=True, type=int)
parser.add_argument("--ncol", help=f"N of columns for subplots, default=N dates", type=int)
parser.add_argument("--nrow", help=f"N of rows for subplots, default=1", type=int)
args = parser.parse_args()

regn      = args.regn if args.regn else regn
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
FDAYS     = args.fday
enmb      = args.enmb
iconc_ref = args.icref.lower()
ncol      = args.ncol if args.ncol is not None else len(args.fday)
nrow      = args.nrow if args.nrow is not None else 1

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

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]

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

pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/sfs_C192mx025_cice_test/expt{enmb:02d}/cice6"

def read_nsidc(regn, dnmb):
  # Error in NSIDC ice concentration fields:
  if regn == 'north':
    NSIDC_err = mtime.datenum_v2([[2024,7,12],[2025,7,27]], ref_day0=False)
  else:
    NSIDC_err = np.array([])
   
  YR, MM, DD = mtime.datevec(dnmb)[:3]

  # Interpolated NSIDC obs fields:
  pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")
  fliceout = f'NSIDC_iconc_interp_mesh025_1080x1440_{YR}{MM:02d}_{regn}.nc'
  dfliceout = os.path.join(pthnsidc,fliceout)
  print(f'Loading interpolated ice conc {dfliceout}')
  with xr.open_dataset(dfliceout) as dsint:
    AI = dsint['ice_conc'].isel(time=DD-1).squeeze()

  if int(dnmb) in NSIDC_err:
    # Error ice conc fields
    AI = AI * np.nan

  return AI

def read_amsr2(regn, dnmb):
  YR, MM, DD = mtime.datevec(dnmb)[:3]
  # AMSR2 ice conc interpolated onto mesh025
  pthamsr = os.path.join(pthdata,f"AMSR2_OSISAF_L4/{regn}/{YR}")
  fliceout = f'AMSR2_iconc_interp_mesh025_1080x1440_{YR}{MM:02d}_{regn}.nc'
  dfliceout = os.path.join(pthamsr,fliceout)
  assert os.path.isfile(dfliceout), f"File not found {dfliceout}"
  print(f'Loading interpolated ice conc {dfliceout}')
  with xr.open_dataset(dfliceout) as dsint:
    AI = dsint['ice_conc'].isel(time=DD-1).squeeze()
      
  return AI

def plot_field(ax1, m, xh, yh, AP_s, clrmp, rmin, rmax, regn, dnmb, sttl):
  YR, MM, DD = mtime.datevec(dnmb)[:3]
  plt.sca(ax1)
  if regn == 'north':
    m.drawparallels(np.arange(60, 90, 10), labels=[0,0,0,0])
  elif regn == 'south':
    m.drawparallels(np.arange(-80, -50, 10), labels=[0,0,0,0])

  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,0])
  m.drawcoastlines()

  img = ax1.pcolormesh(xh, yh, AP_s, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(sttl, fontsize=8)

  return img

clrmp = mclrmps.colormap_difference_negpos(n_white=0)
rmin = -1.
rmax = 1.

clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]

import mod_gfs_cice_anls as mgfscice
importlib.reload(mgfscice)
expt_name = mgfscice.sfs_tests_info(enmb)


if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)

xh, yh = m(hlon,hlat) # GFS coords


print("Plotting ...")
plt.ion()
fig1 = plt.figure(1,figsize=(12, 9))
fig1.clf()  # Clear the figure
axes = fig1.subplots(nrows=nrow, ncols=ncol, squeeze=False)

fig1.subplots_adjust(
    left=0.02,
    right=0.99,
    top=0.95,
    bottom=0.1,
    wspace=0.05,
    hspace=0.1
)

pthglr = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/GLORYS_ithkn_interp_UFSmesh025'

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
  rdate = YR*10000 + MM*100 + DD
  sttl = f"diff iconc SFS expt{enmb:02d} - {iconc_ref.upper()}\nFDAY={fday:02d} {YR}/{MM:02d}/{DD:02d}"

  irow = iplt // ncol
  icol = iplt % ncol
  ax1 = axes[irow, icol]

  if plot_init:
    flinp = f"iceh_ic.{yr0}-{mm0:02d}-{dd0:02d}-{nsec0:05d}.nc"
  else:
    flinp = f"iceh.{yr0}-{mm0:02d}-{dd0:02d}.nc"

  dflice = os.path.join(pthoutp,flinp)
  print(f"Processing {YR}/{MM}/{DD} {hr0:02d}:00, SFS init {YRI}/{MMI:02d}/{DDI:02d} {hrI:02d}:00\n{dflice}")

  print(f"Processing {dflice}")
  # Different ending is handled _d (daily mean) _h (N-hr means)
  # Depending on specification in ice_in, daily mean may have _d or _h
  with xr.open_dataset(dflice) as dcice:
    for varnm in ['aice_h', 'aice_d']:
      if varnm in dcice:
        A2d = dcice[varnm].data.squeeze()
        break
    else:
      raise KeyError(f"No aice_h or aice_d in {dflcice}")

  if iconc_ref == 'nsidc':
    R2d = read_nsidc(regn, dnmb0)
  elif iconc_ref == 'amsr2':
    R2d = read_amsr2(regn, dnmb0)
  else:
    raise RuntimeError(f"Unrecognized reference iconc {iconc_ref}")

  dIT = A2d - R2d
  img = plot_field(ax1, m, xh, yh, dIT, clrmp, rmin, rmax, regn, dnmb0, sttl)
 

sinfo = f'Difference of ice thickness: Model-{iconc_ref.upper()}\n'
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

# Colorbar dimensions
cb_height = 0.018
cb_gap = 0.016

pos_clrb = [
    left,
    bottom - cb_gap - cb_height,
    right - left,
    cb_height
]

# Create colorbar
ax3 = fig1.add_axes(pos_clrb)

clb = fig1.colorbar(
    img,
    cax=ax3,
    orientation='horizontal',
     extend='both'
)

#clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ticks = np.linspace(rmin, rmax, 11)
clb.set_ticks(ticks)
clb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
clb.ax.tick_params(direction='in', length=12, labelsize=12)

fig1.canvas.draw()

pos_clb = ax3.get_position()
bot_clb = pos_clb.y0
pbtm = bot_clb - 0.07


ax4 = fig1.add_axes([0.02, pbtm, 0.8, 0.06])
ax4.text(0, 0, sinfo, fontsize=8,
         ha='left', va='bottom',
         transform=ax4.transAxes)
ax4.axis('off')

pos_info = ax4.get_position()
bot_info = pos_info.y0
pbtm = bot_info - 0.02

btx = f' @{machine}: plot_SFS_diff_ithkn_Nsbplt.py'
#bottom_text(btx, pos=[0.02,0.02], fsz=8)
bottom_text(btx, pos=[0.02, pbtm], fsz=8)



