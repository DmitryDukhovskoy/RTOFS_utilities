"""
  Check restart CICE6
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
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)
  
rest_date = 20250103
rest_hr   = 0
hunits    = 'cm'
  
parser = argparse.ArgumentParser()
parser.add_argument("--rdate", help=f"restart date input file, default={rest_date}", type=int)
parser.add_argument("--rhr", help=f"input file, restart hour = 0, ..., 23, default={rest_hr}", type=int)
parser.add_argument("--rdate_out", help="output file, restart date if different from input", type=int)
parser.add_argument("--rhr_out", help="output file, restart hour if date is different from input", type=int)
args = parser.parse_args()

rest_date     = args.rdate if args.rdate else rest_date
rest_hr       = args.rhr if args.rhr else rest_hr
rest_date_out = args.rdate_out if args.rdate_out else rest_date
rest_hr_out   = args.rhr_out if args.rhr_out else rest_hr

change_rest_time = (rest_date != rest_date_out) or (rest_hr != rest_hr_out)

# Get date numbers:
# Input restart file
dnmbR = mtime.rdate2datenum(rest_date*100+rest_hr)  # restart day nmb
yrR,mmR,ddR,hrR = mtime.datevec(dnmbR, round_hrs=True)[:4]
nsecR = hrR*3600

# Dates of the output fields in the new restart:
dnmbN = mtime.rdate2datenum(rest_date_out*100+rest_hr_out)
yrN, mmN, ddN, hrN = mtime.datevec(dnmbN, round_hrs=True)[:4]
nsecN = hrN*3600 

syst_info = os.uname()
machine = syst_info.nodename

if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
elif 'an' in machine:
  print("Running on PPAN node:", machine)
  node_nm = "ppan"
else:
  print("Unknown machine:", machine)

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

pthrest = os.path.join(pths_ufs[node_nm]["MOM6"]["pthrest"],'new')
flrst = f"cice_model.res.{yrN}{mmN:02d}{ddN:02d}.{hrN:02d}.iconc.nc"
#flrst = 'cice_model.res.20250103.000000.nc'
dflrst = os.path.join(pthrest,flrst)

print(f"Reading {dflrst}")

# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

LON, LAT = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')


ds_rst = xarray.open_dataset(dflrst)
aicen = ds_rst['aicen'].data
vicen = ds_rst['vicen'].data
vsnon = ds_rst['vsnon'].data
ncat = vsnon.shape[0]

# Total ice vol:
vice = np.sum(vicen*aicen, axis=0).squeeze()

# Aggreg iconc:
aice = np.sum(aicen, axis=0).squeeze()

# mean ice thickness over ice area:
hice = np.divide(vice, aice, out=np.zeros_like(aice), where=aice != 0)

# Check hice(n) as it is caclulated in icepack_therm_vertical.F90
# hice(n) = vice(n) / aice(n) 
for k in range(1,ncat):
  aice_n = aicen[k-1,:].squeeze()
  vice_n = vicen[k-1,:].squeeze()
  hice_n = np.divide(vice_n, aice_n, out=np.zeros_like(aice), where=aice_n != 0)
  jmin, imin = np.unravel_index(hice_n.argmin(), hice_n.shape)
  jmax, imax = np.unravel_index(hice_n.argmax(), hice_n.shape)
  print(f"Cat {k}, j={jmin}, i={imin}, min hice(n): {np.nanmin(hice_n)}, "+\
        f"aice(n): {aice_n[jmin,imin]}, vice(n): {vice_n[jmin,imin]}")
  print(f"         j={jmax}, i={imax}, max hice(n): {np.nanmax(hice_n)}, "+\
        f"aice(n): {aice_n[jmax,imax]}, vice(n): {vice_n[jmax,imax]}")

f_plt = False
if f_plt:
  clrmp = mclrmps.colormap_conc()
  rmin = 0.
  rmax = 5.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  hlon = LON
  hlat = LAT

  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
  #lons, lats = m.makegrid(idim, jdim) # get lat/lons of ny by nx evenly spaced grid.
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
  xl1 = -8.e6
  xl2 = -1.2e6
  yl1 = xl1
  yl2 = xl2

  xh, yh = m(hlon,hlat) # GFS coords

  plt.ion()
  fig1 = plt.figure(1, figsize=(8,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
  img1 = ax1.pcolormesh(xh,yh,aice, cmap=clrmp, vmin=rmin, vmax=rmax)

  ax1.contour(xh,yh,HH,[0], linestyles='solid', colors=[(0.,0.,0.)], linewidths=1)

  ax1.set_xlim([xl1, xl2])
  ax1.set_ylim([yl1, yl2])
  ax1.invert_yaxis()
  ax1.invert_xaxis()

  # Plot pnt:
  x0 = hlon[j0,i0]
  y0 = hlat[j0,i0]
  xh0 = xh[j0,i0]
  yh0 = yh[j0,i0]
  ax1.plot(xh0,yh0,'o')

  # Colorbars
  ax3 = fig1.add_axes([0.2, 0.05, 0.6, 0.02])
  clb = plt.colorbar(img1, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)






