"""
  Compare snow thickness in test runs with datm UFS

  Note: snowfall is in liquid water equiv. 
  in CICE, snowfall rate dumped to history is converted from kg/m2*s -->
  cm/day in liquid water equivalent  !!!

"""

NOT FINISHED

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
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob

expt = 'ufs_datm_mx025_v02'
init_date = 20250103
init_hr = 0
regn = 'south'

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", type=str, required=True)
parser.add_argument(f"--dinit", help=f"init date of f/cast, default={init_date}", type=int)
parser.add_argument(f"--ihr", help=f"init hour of f/cast, default={init_hr}", type=int)
parser.add_argument("--enmb1", help="Experiment 1 = 1, ..., ", type=int, required=True)
parser.add_argument("--enmb2", help="Experiment 2 = 1, ...", type=int, required=True)
parser.add_argument("--dfcst", help="Forecast day to plot = 1, ..., 14", type=int, required=True)
args = parser.parse_args()
  
regn = args.regn if args.regn else None
init_date = args.dinit if args.dinit else init_date
init_hr = args.ihr if args.ihr else init_hr
dfcst = args.dfcst if args.dfcst is not None else None
enmb1 = args.enmb1 if args.enmb1 else None
enmb2 = args.enmb2 if args.enmb2 else None



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

pthoutp = pths_ufs[node_nm]["MOM6"]["pthcice"].format(enmb=enmb)

sinfo = pthoutp

varnc = varnm
units = 'cm'
clrmp = mclrmps.colormap_temp()
rmin = 0.
rmax = 20.
sinfo = 'grid cell mean snow thickness\n'

if regn == 'south':
  RMsk = np.where(HH>=0, 0, 1)
  RMsk = np.where(hlat > -60., 0, RMsk)

dnmbS = int(mtime.datenum([YR,MMS,DDS]))
dnmbE = int(mtime.datenum([YR,MME,DDE]))

dlt_day = 86400.
rho_snow = 330.
lwe2snow = 1000./rho_snow # convert m of water equivalent snow fall to m of snow
irec = 0
VFall_sn = []
for dnmb in range(dnmbS,dnmbE+1):
  YR,MM,DD = mtime.datevec(dnmb)[:3]
  print(f"Processing {YR}/{MM}/{DD} ...")
  flcice = f"iceh.{YR}-{MM:02d}-{DD:02d}.nc"
  dflcice = os.path.join(pthoutp,flcice)
  print(f"Processing {dflcice}")
  with xarray.open_dataset(dflcice) as dcice:
    Aice = dcice['aice_d'].data.squeeze()
    fsnow = dcice['snow_ai_d'].squeeze().data*0.01*lwe2snow   # cm/day of water  --> m/day of snow, cell mean
    Acell = dcice['tarea'].data.squeeze()

  fsnow = np.where(RMsk==0, np.nan, fsnow)
  # Convert m3(snow)/m2(cell)*sec --> m3(snow)/m2(ice)*sec 
  fsnow_ice = np.divide(fsnow, Aice, out=np.zeros_like(fsnow), where=Aice > 0)
  # Volume of snowfall over Antarctica
  # to compare different experiments assume ice conc = 1 
  vfs = np.nansum(fsnow_ice * Acell * dlt_day)  # m3/day --> m3 assuming that ice fraction = 1 everywhere
  VFall_sn.append(vfs)


  irec += 1


VFall_sn = np.array(VFall_sn) * 1e-9   # m3 --> km3
XT = np.arange(1, dnmbE-dnmbS+2)
sttl = f"IntegrSnowfall, km3 datm UFS expt{enmb:02d}, {YR}/{MMS:02d}/{DDS:02d}-{YR}/{MME:02d}/{DDE:02d}"
btx  = 'integrated_snowfall_datmUFS.py'

print('Plotting ...')

plt.ion()

clrfsn = [0., 0.5, 0.9]
clrmlt = [0.9, 0.3, 0]
clrvsn = [0.8, 0., 0.6]

clr1 = [0., 0.5, 1]
clr2 = [0.3, 0.9, 0]

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])

ln1, = ax1.plot(XT, VFall_sn, '.-', linewidth=2, color=clrfsn, label='snowfall')

ax1.set_xticks(XT)
ax1.grid('on')
ax1.set_xlabel('Forecast days')
ax1.set_title(sttl)

bottom_text(btx, pos=[0.02, 0.4])


