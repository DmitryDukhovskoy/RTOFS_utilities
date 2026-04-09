"""
  Plot MOM6 surfT on a tri-polar grid
  Show lon/ lat

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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6

init_date = 20250714
init_hr = 0    # nominal hr, actual: -6 hrs for IAU, and -3 FHROT (f/cast hr rotation)

parser = argparse.ArgumentParser()
parser.add_argument("--init", help=f"init date: YYYYMMDD", required=True, type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--fday", help=f"forecast day to plot: 1,...,14, =0 - init. cond.", type=int, required=True)
parser.add_argument("--enmb", help="experiment number: 1, 2, ...", type=int, required=True)
args = parser.parse_args()
  
enmb      = args.enmb
init_date = args.init
init_hr   = args.ihr if args.ihr else init_hr
fday      = args.fday


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


# Get date:

dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
dnmb0 = dnmbI + fday-1                              # day to plot

yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
YR,MM,DD = mtime.datevec(dnmb0)[:3]
nsec0 = hr0*3600

RMsk = np.where(HH>=0, 0, 1)

pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/ufs_datm_mx025/expt{enmb:02d}/mom6"
flinp = f"ocean_{yr0}_{mm0:02d}_{dd0:02d}.nc"
dflmom = os.path.join(pthoutp,flinp)

print(f"Processing {YR}/{MM}/{DD} expt{enmb:02d} {dflmom}...")
ilevel=1
with xarray.open_dataset(dflmom) as dmom:
  AA = dmom['potT'].isel(time=0, zl=ilevel).data.squeeze()    # grid cell mean ice thickness

AA = np.where(RMsk == 0, np.nan, AA)
clrmp = mclrmps.colormap_temp(skip_dark=0.1)
rmin = -2
rmax = 32
clrmp.set_bad(color=[0.1, 0.1, 0.1])

hlon1 = hlon.copy()
hlon1[hlon1<-5] = np.nan

hlon2 = np.where(hlon > 60, hlon-360., hlon)
hlon2[hlon2>=5] = np.nan

plt.ion()
fgnmb = 1
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
img  = ax1.pcolormesh(AA, cmap = clrmp, vmin=rmin, vmax=rmax)

clr_grid = (0.2,0.2,0.2)
lon_levels1 = np.arange(0, 180, 10)
lon_levels2 = np.arange(-180, 0, 10)
ax1.contour(hlon1, levels=lon_levels1, 
            linestyles='solid', linewidths=0.5, colors=[clr_grid])
ax1.contour(hlon2, levels=lon_levels2, 
            linestyles='solid', linewidths=0.5, colors=[clr_grid])

# Latitude grid lines every 10 degrees
lat_levels = np.arange(-80, 90, 10)
ax1.contour(hlat, levels=lat_levels, 
            linestyles='solid', linewidths=0.5, colors=[clr_grid])

ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xticklabels([]) 
ax1.set_yticklabels([])  

ax1.axis('scaled')

ax1.set_title(f'GFSv17 MOM6 daily T {YR}/{MM:02d}/{DD:02d} expt{enmb:02d} layer {ilevel:02d} \nmesh025 on a tripolar grid')

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='both')
ticks = np.linspace(rmin, rmax, 11)
clb.set_ticks(ticks)
clb.ax.set_xticklabels([f"{t:.0f}" for t in ticks], fontsize=12)
clb.ax.tick_params(direction='in', length=12)


btx = 'plot_surfTmom6_grid.py'
bottom_text(btx)




