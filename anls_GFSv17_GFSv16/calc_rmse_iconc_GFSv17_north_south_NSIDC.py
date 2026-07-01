"""
  RMSE of ice conc btw GFSv17 and NRT NSIDC
  both Arctic and Antarctic are shown on the same plot
  E.g., show winter seasons for both regions

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
    
import importlib
import mod_utils_fig
importlib.reload(mod_utils_fig)
bottom_text = mod_utils_fig.bottom_text

import mod_time as mtime
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import pandas as pd
  
init_hr = 0
fdays = 16    # N of f/cast days
prst = 0
fld_name = 'iconc'   
 
parser = argparse.ArgumentParser()
parser.add_argument("--initn", help=f"init date: North", choices=[20251231, 20240715], required=True, type=int)
parser.add_argument("--inits", help=f"init date: South", choices=[20251231, 20240715], required=True, type=int)
parser.add_argument("--fdays", help=f"N of f/cast days, default={fdays}", type=int)
parser.add_argument("--prst", help="Show persistence for both regions (1), or not (0)", 
                    choices=[0,1], required=True, type=int)
args = parser.parse_args()
  
init_north = args.initn 
init_south = args.inits
fdays = args.fdays if args.fdays else fdays
prst = args.prst 

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

def calc_rmse(init_date, fdays, track_prst, regn, RMsk):
  pthgfs17 = f"/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/GFSv17/{init_date}/daily"
  fgfs17   = f"gfsv17_iconc_init{init_date}_days000_016.nc"
  dfgfs17  = os.path.join(pthgfs17, fgfs17)

  if regn == 'south':
    RMsk = np.where(hlat > -60., 0, RMsk)
  elif regn == 'north':
    RMsk = np.where(hlat < 50, 0, RMsk)


  varnm = 'aice_h'
  with xarray.open_dataset(dfgfs17) as ds17:
    A3d_gfs17 = ds17[varnm].values
    TM17 = pd.to_datetime(ds17['time'].data)

  A17p = None
  RMSEp17 = []
  RMSE17  = []
  TM = []
  for iday in range(0,fdays):
    YR = TM17[iday].year
    MM = TM17[iday].month
    DD = TM17[iday].day
    HR = TM17[iday].hour

    A17 = A3d_gfs17[iday,:,:].squeeze()

    # Interpolated NSIDC obs fields:
    pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")
    fliceout = f'NSIDC_iconc_interp_mesh025_{jdm}x{idm}_{YR}{MM:02d}_{regn}.nc'
    dfliceout = os.path.join(pthnsidc,fliceout)
    print(f'Loading interpolated ice conc {dfliceout}')
    with xarray.open_dataset(dfliceout) as dsint:
      AI = dsint['ice_conc'].isel(time=DD-1).squeeze()

    A17 = np.where(RMsk == 0, np.nan, A17)
    AI = np.where(RMsk == 0, np.nan, AI)
    rmse17 = rmse2d(A17, AI, Acell)

    if track_prst:
      if A17p is None:
        A17p = A17.copy()
      rmse_p17 = rmse2d(A17p, AI, Acell)
      RMSEp17.append(rmse_p17)
      print(f"RMSE 17={rmse17:.3f}  RMSE_prst 17={rmse_p17:.3f}")

    RMSE17.append(rmse17)
    TM.append(mtime.datenum([YR,MM,DD]))

  RMSE17  = np.array(RMSE17)
  RMSEp17 = np.array(RMSEp17)
  TM = np.array(TM)

  return RMSE17, RMSEp17, TM

RMSE_N, RMSEp_N, TM_N = calc_rmse(init_north, fdays, track_prst, 'north', RMsk)
RMSE_S, RMSEp_S, TM_S = calc_rmse(init_south, fdays, track_prst, 'south', RMsk)

# Plot 

XTn = (TM_N - TM_N[0]) + 1    # lead time, days
Xplt_N = XTn-0.5            # dayly avrg
XTs = (TM_S - TM_S[0]) + 1    # lead time, days
Xplt_S = XTs-0.5            # dayly avrg

xticks = np.arange(np.floor(XTn[0]),np.ceil(XTn[-1]+1))
yticks = np.arange(0.,1,0.05)
sttl = f"RMSE iconc NRT NSIDC GFSv17 north: {init_north}, south: {init_south}\n"

clrN = [0., 0.6, 1]
clrS = [1, 0.3, 0]

yl1 = 0
if track_prst:
  yl2 = max([np.max(RMSE_N), np.max(RMSE_S), np.max(RMSEp_N), np.max(RMSEp_S)]) * 1.3
else:
  yl2 = max([np.max(RMSE_N), np.max(RMSE_S)]) * 1.3

# To keep yl2 same for Arc/ S. Ocean:
if yl2 > 0.5:
  yl2 = 0.8
else:
  yl2 = 0.35

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()

ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])
ln1, = ax1.plot(Xplt_N, RMSE_N, 'o-', linewidth=2, color=clrN, label='north')
ln2, = ax1.plot(Xplt_S, RMSE_S, 'o-', linewidth=2, color=clrS, label='south')
if track_prst:
  ln3, = ax1.plot(Xplt_N, RMSEp_N, '--', linewidth=2, color=clrN, label='persist north')
  ln4, = ax1.plot(Xplt_S, RMSEp_S, '--', linewidth=2, color=clrS, label='persist south')

ax1.set_yticks(yticks)
ax1.set_xticks(xticks)
ax1.set_ylim(yl1, yl2)
ax1.grid('on')
ax1.set_ylabel('Ice partial area')
ax1.set_xlabel('Forecast days')
ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)

ax1.set_title(sttl)

ax3 = plt.axes([0.1, 0.15, 0.6, 0.2])
if track_prst:
  LNS = [ln1, ln2, ln3, ln4]
else:
  LNS = [ln1, ln2]

lgd = plt.legend(handles=LNS, loc='upper left', fontsize=14)
ax3.axis('off')

btx = 'calc_rmse_iconc_GFSv17_north_south_NSIDC.py'
bottom_text(btx, pos=[0.02,0.15])













