"""
  Check snow depth at a grid point

  From a debug run with saved variables in ice thickness categories

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
import mod_misc1 as mmisc

init_date = 20250701
init_hr = 0

# 61 - debug run
enmb = 61

parser = argparse.ArgumentParser()
parser.add_argument("--init", help=f"init date", choices=[20240701, 20250701], default=init_date, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fday1", help=f"forecast day to plot: 0, 1, 2, ... , =0 - init. cond.", type=int, required=True)
parser.add_argument("--fday2", help=f"forecast day to plot: 0, 1, 2, ... , =0 - init. cond.", type=int, required=True)
parser.add_argument("--enmb", help="experiment number: 0, 1, 2, ...", default=enmb, type=int)
args = parser.parse_args()

init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fday1     = args.fday1
fday2     = args.fday2
enmb      = args.enmb

def read_cice(fday, HH, init_date, init_hr):
# Get date:
  plot_init = fday == 0  # initial conditions
  
  dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
  YRI,MMI,DDI,hrI = mtime.datevec(dnmbI, round_hrs=True)[:4]

  if plot_init:
    dnmb0 = dnmbI
  else:
    dnmb0 = dnmbI + fday-1                           # day to plot

  yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
  nsec0 = int(hr0 * 3600)
  YR,MM,DD = mtime.datevec(dnmb0)[:3]

  pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/sfs_C192mx025_cice_test/expt{enmb:02d}/cice6"
  if plot_init:
    flinp = f"iceh_ic.{yr0}-{mm0:02d}-{dd0:02d}-{nsec0:05d}.nc"
  else:
    flinp = f"iceh.{yr0}-{mm0:02d}-{dd0:02d}.nc"
  varnm = 'hs_d'

  dflice = os.path.join(pthoutp,flinp)

  print(f"Reading {YR}/{MM}/{DD}, SFS init {YRI}/{MMI:02d}/{DDI:02d}\n  {dflice}")

  with xarray.open_dataset(dflice) as dcice:
    vsno_m2c = dcice['hs_d'].data.squeeze()   # grid cell mean snow vol, m3/m2_cell
    aicen = dcice['aicen_d'].values.squeeze()  # ice area in cats
    vicen = dcice['vicen_d'].values.squeeze()  # ice vol m2_ice in cats
    vsnon = dcice['vsnon_d'].values.squeeze()  # snow vol m2_ice in cats
    snfrn = dcice['snowfracn_d'].values.squeeze() # snow frac m2_cell, in cats

  vsno_m2c[HH >= 0] = np.nan

  return vsno_m2c, aicen, vicen, vsnon, snfrn


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

vsno1_m2c, aicen1, vicen1, vsnon1, snfrn1 = read_cice(fday1, HH, init_date, init_hr)
vsno2_m2c, aicen2, vicen2, vsnon2, snfrn2 = read_cice(fday2, HH, init_date, init_hr)

clrmp = mclrmps.colormap_temp()
rmin = 0.
rmax = 0.1

clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]

plt.ion()

print("Plotting ...")

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])

ax1.pcolormesh(vsno2_m2c, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

ax1.set_aspect('equal', adjustable='box')
ax1.set_ylim(800, 1080)
ax1.set_xlim(100, 670)

# Grid point analysis
ii0 = 482
jj0 = 1062

aice1 = np.sum(aicen1[:,jj0,ii0])
aice2 = np.sum(aicen2[:,jj0,ii0])
hicen1 = vicen1[:,jj0,ii0] / aicen1[:,jj0,ii0]
hicen2 = vicen2[:,jj0,ii0] / aicen2[:,jj0,ii0]
vsno1 = np.sum(vsnon1[:,jj0,ii0] * aicen1[:,jj0,ii0])
vsno2 = np.sum(vsnon2[:,jj0,ii0] * aicen2[:,jj0,ii0])


str = f"vsnon j={jj0} i={ii0}:"
mmisc.print_cols(vsnon1[:,jj0,ii0], vsnon2[:,jj0,ii0], str=str)

str = f"aicen j={jj0} i={ii0}:"
mmisc.print_cols(aicen1[:,jj0,ii0], aicen2[:,jj0,ii0], str=str)

str = f"vicen j={jj0} i={ii0}:"
mmisc.print_cols(vicen1[:,jj0,ii0], vicen2[:,jj0,ii0], str=str)

str = f"hice j={jj0} i={ii0}:"
mmisc.print_cols(hicen1, hicen2, str=str)





