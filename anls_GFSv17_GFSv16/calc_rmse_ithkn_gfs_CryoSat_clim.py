"""
  RMSE of ice conc btw GFSv16, GFSv17 and 
  ice thickness climatology from monthly CryoSat AWI ice thickness fields

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
fld_name = 'ithkn'   
 
parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere to analyze", type=str, 
                    choices=['north','south'], required=True)
parser.add_argument("--init", help=f"init date, default={init_date}", type=int)
parser.add_argument("--fdays", help=f"N of f/cast days, default={fdays}", type=int)
args = parser.parse_args()
  
regn  = args.regn if args.regn else None
init_date = args.init if args.init else init_date
fdays = args.fdays if args.fdays else fdays


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
  A3d_gfs16 = ds16['ICETK_surface'].values
  TM16 = pd.to_datetime(ds16['time'].data)

varnm = 'hi_h'
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

  # CryoSat AWI ithkn climatology:
  if regn == 'north':
    pthithkn = os.path.join(pthdata,'CryoSat_arctic_ice_snow_thkn','clim')
    fclim = 'ithkn_CryoSat_arcticAWI_mnthclim_2015-2024_1080x1440.nc'

  dfliceout = os.path.join(pthithkn, fclim)
  print(f'Reading ice thkn clim {dfliceout}')
  with xarray.open_dataset(dfliceout) as dsint:
    AI = dsint['ithkn'].isel(time=MM-1).squeeze()

  A16 = np.where(RMsk == 0, np.nan, A16)
  A17 = np.where(RMsk == 0, np.nan, A17)
  AI = np.where(RMsk == 0, np.nan, AI)
  rmse16 = rmse2d(A16, AI, Acell)
  rmse17 = rmse2d(A17, AI, Acell)

  print(f"RMSE 16/17={rmse16:.3f}/{rmse17:.3f}")

  RMSE16.append(rmse16)
  RMSE17.append(rmse17)
  TM.append(mtime.datenum([YR,MM,DD]))

# Ice conc does not change in GFSv16 --> rmse = rmse(persistence)
RMSE16  = np.array(RMSE16)
RMSE17  = np.array(RMSE17)
TM = np.array(TM)

# Plot

XT = (TM - TM[0]) + 1    # lead time, days
Xplt = XT-0.5            # dayly avrg
xticks = np.arange(np.floor(XT[0]),np.ceil(XT[-1]+1))
yticks = np.arange(0.,0.8,0.05)
sttl = f"RMSE ithkn CryoSat clim and GFSv16, GFSv17 {init_date}, {regn}\n"

clr16 = [0., 0.6, 1]
clr17 = [0.9, 0.3, 0]
yl1 = 0
yl2 = max([np.max(RMSE16), np.max(RMSE17)]) * 1.1


plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()

ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])
ln1, = ax1.plot(Xplt, RMSE16, 'o-', linewidth=2, color=clr16, label='GFSv16')
ln2, = ax1.plot(Xplt, RMSE17, 'o-', linewidth=2, color=clr17, label='GFSv17')

ax1.set_yticks(yticks)
ax1.set_xticks(xticks)
ax1.set_ylim(yl1, yl2)
ax1.grid('on')
ax1.set_xlabel('Forecast days')
ax1.set_title(sttl)

ax3 = plt.axes([0.1, 0.15, 0.6, 0.2])
LNS = [ln1, ln2]
lgd = plt.legend(handles=LNS, loc='upper left')
ax3.axis('off')

btx = 'calc_rmse_iconc_gfs_NSIDC.py'
bottom_text(btx, pos=[0.1,0.1])













