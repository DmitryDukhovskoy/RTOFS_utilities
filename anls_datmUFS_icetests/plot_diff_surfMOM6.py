"""
  Plot difference of ocean surface fields from sensitivity tests with
  atm.-forced UFS (datm UFS) between 2 experiments

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


expt = 'ufs_datm_mx025_v02'
init_date = 20250103
init_hr = 0
regn = 'south'

parser = argparse.ArgumentParser()
parser.add_argument("--fldnm", help="Field name to analyze", 
                    choices=['sst','sss'], required=True, type=str)
parser.add_argument("--enmb", help="2 experiment numbers to compare",
                    type=int, nargs="+", required=True)
parser.add_argument("--init", help=f"init date, default={init_date}", type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--fday", help=f"forecast day to plot: 1,...,14", 
                    type=int, required=True)
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
args = parser.parse_args()

ENMBS     = args.enmb if args.enmb else None
init_date = args.init if args.init else init_date
init_hr   = args.ihr  if args.ihr else init_hr
fday      = args.fday if args.fday > 0 else None
regn      = args.regn if args.regn else regn
fld_name  = args.fldnm 

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


dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
dnmb0 = dnmbI + fday-1                              # day to plot

yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
YR,MM,DD = mtime.datevec(dnmb0)[:3]
nsec0 = hr0*3600

clrmp = mclrmps.colormap_uv()
if fld_name == 'sst':
  varnm = 'potT'
  units = 'dgrC'
  rmin = -1
  rmax = 1
elif fld_name == 'sss':
  varnm = 'salt'
  units = 'psu'
  rmin = -1
  rmax = 1

clrmp.set_bad(color=[0.2, 0.2, 0.2])

flinp = f"ocean_{yr0}_{mm0:02d}_{dd0:02d}.nc"

import mod_gfs_cice_anls as mgfscice
Alist  = []
AIlist = []
LBL    = []
for enmb in ENMBS:
  pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/ufs_datm_mx025/expt{enmb:02d}/mom6"
  dflmom = os.path.join(pthoutp,flinp)

  print(f"Processing {YR}/{MM}/{DD} expt{enmb:02d} {dflmom}...")
  with xarray.open_dataset(dflmom) as dcice:
    A2d  = dcice[varnm].isel(zl=0).data.squeeze()

  Alist.append(A2d)
  line_lbl = mgfscice.sens_tests_info(enmb)
  LBL.append(line_lbl)

  # Get sea ice edge:
  pthout_cice = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/ufs_datm_mx025/expt{enmb:02d}/cice6"
  flcice = f"iceh.{YR}-{MM:02d}-{DD:02d}.nc"
  dflcice = os.path.join(pthout_cice,flcice)
  print(f"Processing {dflcice}")
  with xarray.open_dataset(dflcice) as dcice:
    Aice = dcice['aice_d'].data.squeeze()

  AIlist.append(Aice)

dFld = Alist[1] - Alist[0]



sttl = f'diff {fld_name}, datmUFS expt{ENMBS[1]:02d} vs expt{ENMBS[0]:02d}\n'
sttl = sttl + f'expt{ENMBS[1]:02d}={LBL[1]}, expt{ENMBS[0]:02d}={LBL[0]}\n'
sttl = sttl + f'init:{init_date}/{init_hr} fcast:{YR}/{MM:02d}/{DD:02d}'

sinfo = ''
#sinfo = 'difference expt{ENMBS[1]:02d} ({LBL[1]}) - expt{ENMBS[0]:02d} ({LBL[0]}))\n'
sinfo = sinfo + dflmom

plt.ion()

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)

meridians = np.arange(-360,359.,45.)

xh, yh = m(hlon, hlat) # MOM6 coords


fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])
m.drawcoastlines()

m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)

img = ax1.pcolormesh(xh, yh, dFld, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

# Show ice contours:
#CLRS = mgfscice.sens_tests_colors()
CLRS = [[0., 1., 0.4],
        [1, 0., 0.]]

LNS = []
for icc, enmb in enumerate(ENMBS):
  line_lbl = mgfscice.sens_tests_info(enmb)
  clr = CLRS[icc]
  ax1.contour(xh, yh, AIlist[icc], [0.15], linestyles='solid', colors=[clr], linewidths=1)
  # Create legend proxy
  ln1 = matplotlib.lines.Line2D(
      [], [], color=clr,
      linestyle='solid',
      linewidth=1,
      label=line_lbl
  )

  LNS.append(ln1)


ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
if rmin < 0:
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
else:
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')

ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.03, 0.8, 0.06])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')

ax4 = plt.axes([0.63, 0.02, 0.37, 0.07])
lgd = plt.legend(handles=LNS, loc='upper left')
ax4.axis('off')

btx = 'plot_diff_surfMOM6.py'
bottom_text(btx, pos=[0.02, 0.01])


