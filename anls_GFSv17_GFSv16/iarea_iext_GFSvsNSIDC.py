"""
Calculate and plot ice extent or ice area
for north/ south polar regions
from GFSv17 vs NSIDC or AMSR2
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray as xr
import datetime
import pandas as pd
from yaml import safe_load
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

fdays = 16

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", type=str, required=True)
parser.add_argument("--fdays", help=f"N of f/cast days, default={fdays}", type=int)
parser.add_argument(
  "--init", 
  help=f"Start of the f/cast YYYYMMDD",
  choices=[20240715, 20251231], 
  type=int
)
args = parser.parse_args()

regn = args.regn
fdasy = args.fdays
init_date = args.init

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

HH = np.where(HH < 1.e-20, 1, -HH)
jdm, idm = HH.shape
DX, DY = mmom6.dx_dy(hlon, hlat)
Acell = DX*DY*1e-6  # km2

# Region mask:
if regn == 'south':
  RMsk = (hlat < -55.) & (HH < 0)
elif regn == 'north':
  RMsk = (hlat > 60.) & (HH < 0)

def ice_extent_area(aice, RMsk, Acell, iedge=0.15):
  IMsk = (aice > iedge) & RMsk
  IA = np.sum(aice[IMsk] * Acell[IMsk]) #ice area
  IE = np.sum(Acell[IMsk])

  return IA, IE


def read_nsidc(YR,MM,DD):
  # Interpolated NSIDC obs fields:
  pthdata = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data'
  pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")
  fliceout = f'NSIDC_iconc_interp_mesh025_1080x1440_{YR}{MM:02d}_{regn}.nc'
  dfliceout = os.path.join(pthnsidc,fliceout)
  print(f'Loading interpolated ice conc {dfliceout}')
  with xr.open_dataset(dfliceout) as dsint:
    AI = dsint['ice_conc'].isel(time=DD-1).values.squeeze()

  return AI


# Read GFSv17 
pthgfs17 = f"/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/GFSv17/{init_date}/daily"
fgfs17   = f"gfsv17_iconc_init{init_date}_days000_016.nc"
dfgfs17  = os.path.join(pthgfs17, fgfs17)

varnm = 'aice_h'
with xr.open_dataset(dfgfs17) as ds17:
  A3d_gfs17 = ds17[varnm].values
  TM17 = pd.to_datetime(ds17['time'].data)

nrecs = len(TM17)
fdays = min([nrecs, fdays])
IAREA = []
IEXT = []
TM = []
for iday in range(fdays):
  YR = TM17[iday].year
  MM = TM17[iday].month
  DD = TM17[iday].day

  A17 = A3d_gfs17[iday,:,:].squeeze()
  iarea, iext = ice_extent_area(A17, RMsk, Acell)

  AI = read_nsidc(YR,MM,DD)
  iarea_nsidc, iext_nsidc = ice_extent_area(AI, RMsk, Acell)

  IAREA.append((iarea, iarea_nsidc))
  IEXT.append((iext, iext_nsidc))
  TM.append(mtime.datenum([YR,MM,DD]))


IAREA = np.asarray(IAREA) * 1e-6 # km2 * 1e6
IEXT  = np.asarray(IEXT) * 1e-6
TM    = np.asarray(TM)


# Plot 
XTn = (TM - TM[0]) + 1    # lead time, days
Xplt_N = XTn-0.5            # dayly avrg

xticks = np.arange(np.floor(XTn[0]),np.ceil(XTn[-1]+1))
yticks = np.arange(0.,1,0.05)
sttl1 = f"Ice Extent, km2 x 1e6 {regn} GFSv17 {init_date}"
sttl2 = f"Ice Area, km2 x 1e6 {regn} GFSv17 {init_date}"

clr_gfs = [0., 0.5, 0.8]
clr_obs = [1, 0.4, 0]

yext_upper = 1.1 * np.max(IEXT)
yarea_upper = 1.1 * np.max(IAREA)

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()

ax1 = plt.axes([0.12, 0.6, 0.8, 0.35])
ln1, = ax1.plot(Xplt_N, IEXT[:,0], 'o-', linewidth=2, color=clr_gfs, label='GFSv17')
ln2, = ax1.plot(Xplt_N, IEXT[:,1], 'o-', linewidth=2, color=clr_obs, label='NSIDC')
ax1.set_xticks(xticks)
ax1.set_ylim(bottom=0, top=yext_upper)
ax1.grid('on')
ax1.set_ylabel('Ice Extent x 1e6, km2')
ax1.set_xlabel('Forecast days')
ax1.tick_params(axis='x', labelsize=12)
ax1.tick_params(axis='y', labelsize=12)
ax1.set_title(sttl1)

# Ice Area
ax2 = plt.axes([0.12, 0.16, 0.8, 0.35])
ax2.plot(Xplt_N, IAREA[:,0], 'o-', linewidth=2, color=clr_gfs)
ax2.plot(Xplt_N, IAREA[:,1], 'o-', linewidth=2, color=clr_obs)
ax2.set_xticks(xticks)
ax2.set_ylim(bottom=0, top=yarea_upper)
ax2.grid('on')
ax2.set_ylabel('Ice Area x 1e6, km2')
ax2.set_xlabel('Forecast days')
ax2.tick_params(axis='x', labelsize=12)
ax2.tick_params(axis='y', labelsize=12)
ax2.set_title(sttl2)

ax3 = plt.axes([0.75, 0.03, 0.25, 0.2])
lgd = plt.legend(handles=[ln1,ln2], loc='lower left')
ax3.axis('off')

btx = 'iarea_iext_GFSvsNSIDC.py'
bottom_text(btx, pos=[0.02,0.02])




