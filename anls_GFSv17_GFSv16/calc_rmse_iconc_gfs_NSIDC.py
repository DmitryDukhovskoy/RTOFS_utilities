"""
  RMSE of ice conc btw GFSv16, GFSv17 and NRT NSIDC

  all fields are on mesh025 grid

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
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
import pandas as pd
  
init_date = 20251231  # initialization date
init_hr = 0
fdays = 16    # N of f/cast days
prst = 0
fld_name = 'iconc'   
 
parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere to analyze", type=str, 
                    choices=['north','south'], required=True)
parser.add_argument("--init", help=f"init date", choices=[20251231, 20240715], required=True, type=int)
parser.add_argument("--fdays", help=f"N of f/cast days, default={fdays}", type=int)
parser.add_argument("--fv16", help="Show GFSv16 (1 default), no (0)",
                    choices=[0,1],
                    default=1,
                    type=int)
parser.add_argument(
    "--prst", 
    help=f"Show persitance using GFSv16 or GFSv17 IC, 0=no, 1=yes",
    type=int,
    required=True,
    choices=[0,1]
) 
args = parser.parse_args()
show_v16 = args.fv16 == 1
  
regn  = args.regn if args.regn else None
init_date = args.init if args.init else init_date
fdays = args.fdays if args.fdays else fdays
prst = args.prst if args.prst is not None else prst

track_prst = prst == 1

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

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape

DX, DY = mmom6.dx_dy(hlon, hlat)
Acell = DX*DY*1e-6  # km2

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]

RMsk = np.where(HH>=0, 0, 1)
if regn == 'south':
  RMsk = np.where(hlat > -60., 0, RMsk)
elif regn == 'north':
  RMsk = np.where(hlat < 50, 0, RMsk)

def rmse2d(AA, AI, Acell, wgt_area=True):
  """
    RMSE of 2d fields
  """
  mask = ((AA > 1.e-3) | (AI > 1.e-3)) & ~np.isnan(AA) & ~np.isnan(AI)
  assert np.any(mask), "No ice grid found"

  sqerr = (AA - AI)**2
  if wgt_area:
    rmse = np.sqrt(np.sum(Acell[mask] * sqerr[mask]) / np.sum(Acell[mask]))
  else:
    rmse = np.sqrt(np.mean(sqerr[mask]))

  return rmse

pthgfs16 = f"/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/GFSv16/{init_date}/intrp_mesh025"
pthgfs17 = f"/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/GFSv17/{init_date}/daily"
fgfs16   = f"gfsv16_{fld_name}_init{init_date}_days000_016.nc"
fgfs17   = f"gfsv17_{fld_name}_init{init_date}_days000_016.nc"
dfgfs16  = os.path.join(pthgfs16, fgfs16)
dfgfs17  = os.path.join(pthgfs17, fgfs17)

with xarray.open_dataset(dfgfs16) as ds16:
  A3d_gfs16 = ds16['ICEC_surface'].values
  TM16 = pd.to_datetime(ds16['time'].data)

varnm = 'aice_h'
with xarray.open_dataset(dfgfs17) as ds17:
  A3d_gfs17 = ds17[varnm].values
  TM17 = pd.to_datetime(ds17['time'].data)

A16p = A17p = None
RMSEp16 = []
RMSEp17 = []
RMSE16  = []
RMSE17  = []
TM = []
for iday in range(0,fdays):
  assert TM16[iday] == TM17[iday], f"Mismatched dates in GFSv16 and GFSv16 {iday}"
  YR = TM17[iday].year
  MM = TM17[iday].month
  DD = TM17[iday].day
  HR = TM17[iday].hour

  A16 = A3d_gfs16[iday,:,:].squeeze()
  A17 = A3d_gfs17[iday,:,:].squeeze()

  # Interpolated NSIDC obs fields:
  pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")
  fliceout = f'NSIDC_iconc_interp_mesh025_{jdm}x{idm}_{YR}{MM:02d}_{regn}.nc'
  dfliceout = os.path.join(pthnsidc,fliceout)
  print(f'Loading interpolated ice conc {dfliceout}')
  with xarray.open_dataset(dfliceout) as dsint:
    AI = dsint['ice_conc'].isel(time=DD-1).squeeze()

  A16 = np.where(RMsk == 0, np.nan, A16)
  A17 = np.where(RMsk == 0, np.nan, A17)
  AI = np.where(RMsk == 0, np.nan, AI)
  rmse16 = rmse2d(A16, AI, Acell)
  rmse17 = rmse2d(A17, AI, Acell)

  if track_prst:
    if A17p is None:
      A16p = A16.copy()
      A17p = A17.copy()
    rmse_p16 = rmse2d(A16p, AI, Acell)
    RMSEp16.append(rmse_p16)
    rmse_p17 = rmse2d(A17p, AI, Acell)
    RMSEp17.append(rmse_p17)
    print(f"RMSE 16/17={rmse16:.3f}/{rmse17:.3f}  RMSE_prst 16/17={rmse_p16:.3f}/{rmse_p17:.3f}")
  else:
    print(f"RMSE 16/17={rmse16:.3f}/{rmse17:.3f}")

  RMSE16.append(rmse16)
  RMSE17.append(rmse17)
  TM.append(mtime.datenum([YR,MM,DD]))

# Ice conc does not change in GFSv16 --> rmse = rmse(persistence)
RMSE16  = np.array(RMSE16)
RMSE17  = np.array(RMSE17)
RMSEp17 = np.array(RMSEp17)
TM = np.array(TM)

# Plot

XT = (TM - TM[0]) + 1    # lead time, days
Xplt = XT-0.5            # dayly avrg
xticks = np.arange(np.floor(XT[0]),np.ceil(XT[-1]+1))
yticks = np.arange(0.,1,0.05)
sttl = f"RMSE iconc NRT NSIDC and GFSv16, GFSv17 {init_date}, {regn}\n"

#clr16 = [0., 0.6, 1]
clr16 = [0., 0, 0]
clr17 = [1, 0., 0]
yl1 = 0
if track_prst:
  yl2 = max([np.max(RMSE16), np.max(RMSE17), np.max(RMSEp17)]) * 1.3
else:
  yl2 = max([np.max(RMSE16), np.max(RMSE17)]) * 1.3

# To keep yl2 same for Arc/ S. Ocean:
if yl2 > 0.5:
  yl2 = 0.8
else:
  yl2 = 0.35

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()

ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])
ln1 = []
if show_v16:
  ln1, = ax1.plot(Xplt, RMSE16, 'o-', linewidth=2, color=clr16, label='GFSv16')
ln2, = ax1.plot(Xplt, RMSE17, 'o-', linewidth=2, color=clr17, label='GFSv17')
ln3, = ax1.plot(Xplt, RMSEp17, '--', linewidth=2, color=clr17, label='persist')

ax1.set_yticks(yticks)
ax1.set_xticks(xticks)
ax1.set_ylim(yl1, yl2)
ax1.grid('on')
ax1.set_ylabel('Ice partial area')
ax1.set_xlabel('Forecast days')
ax1.set_title(sttl)

ax3 = plt.axes([0.1, 0.15, 0.6, 0.2])
if show_v16:
  LNS = [ln1, ln2, ln3]
else:
  LNS = [ln2,ln3]

lgd = plt.legend(handles=LNS, loc='upper left')
ax3.axis('off')

btx = 'calc_rmse_iconc_gfs_NSIDC.py'
bottom_text(btx, pos=[0.1,0.1])













