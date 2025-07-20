"""
  Plot NSIDC ice conc climatology 
  calculated in calc_NSIDC_iconc_clim.py
  Plot on 12 subplots
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
import argparse

PPTHN = '/home/Dmitry.Dukhovskoy/python'
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="Start year used for averaging", type=int, required=True)
parser.add_argument("--YRE", help="End year of averaging", type=int, required=True)
parser.add_argument("--MMS", help="Cal. month start to plot, default 1", type=int)
parser.add_argument("--MME", help="Cal. month end to plot, default 12", type=int)
parser.add_argument("--ncol", help="N of columns for subplots", type=int)
parser.add_argument("--nrow", help="N of rows for subplots", type=int)
args = parser.parse_args()

# Default values: 
expt_nmb  = 2
trlx = 12   # max relax. time scale, hrs
expt_name = f'NEPbgc_nudged_hindcast{expt_nmb:02d}'
hcst_time = 3 # f/csat time interval, months
hcst_interv = np.array([x for x in range(1,12+hcst_time,hcst_time)], dtype=int)

YRS  = args.YRS if args.YRS else None
YRE  = args.YRE if args.YRE else YRS
MMS  = args.MMS if args.MMS else 1
MME  = args.MME if args.MME else 12
nmnths = MME-MMS+1
ncol = args.ncol if args.ncol else None
nrow = args.nrow if args.nrow else nmnths

if not ncol:
  ncol = min([nmnths,4])
  nrow = nmnths // ncol

assert ncol*nrow == nmnths, f'Specified columns/rows ({ncol}/{nrow}) do not match N months {nmnths}'

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

run_name   = 'seasonal_fcst_daily'
pthtopo    = gridfls['MOM6_NEP'][run_name]['pthgrid']
fgrid      = gridfls['MOM6_NEP'][run_name]['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"][run_name]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

#if varnm == 'iconc':
clrmp = mclrmps.colormap_conc()
clrmp.set_bad(color=[0.2, 0.2, 0.2])
rmin = 0.
rmax = 1.

ci0=0.15
cntr_clr = [1, 0.4, 0]
#cntr_clr = [0.2, 0.7, 1.0]
cntr_nsidc = [0.,0.2,0.6]
cntr_irlx = [1,0,1]

def plot_field(ax1, fig1, m, xR, yR, A2d, clrmp, rmin, rmax, plt_clrb, sttl=[]):
  fig1.sca(ax1)
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(sttl)

  # extend: min, max, both
  if plt_clrb:
    ax2 = fig1.add_axes([0.9,0.1,0.013,0.8])
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')

    ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
    ax2.set_yticklabels(ax2.get_yticks())
    ticklabs = clb.ax.get_yticklabels()
    #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
    clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
    clb.ax.tick_params(direction='in', length=12)

  return ax1

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

plt.ion()
fig1 = plt.figure(1,figsize=(15, 10))
fig1.clf()  # Clear the figure
axes = fig1.subplots(nrows=nrow, ncols=ncol)

fig1.subplots_adjust(
    left=0.05,
    right=0.85,  # More right-side room
    top=0.95,
    bottom=0.1,
    wspace=0.05,
    hspace=0.1
)

iplt = 0
for MMA in range(MMS,MME+1):
  print(f"Processing clim {YRS}-{YRE} {MMA} NSIDC IceConc ...")
  imo = MMA-1      # current month in the archive output

  pthclim = '/work/Dmitry.Dukhovskoy/data/NRT_NOAA_NSIDC_seaconc/clim'
  flnm = f'NSIDC_NRT_interpNEP_iconc_clim_{YRS}_{YRE}.nc'

  dcice = os.path.join(pthclim,flnm)
  print(f'Reading {dcice}')

  ds = xarray.open_dataset(dcice)
  A2d = ds['ice_conc'].isel(months=imo).data.squeeze()

  irow = iplt // ncol
  icol = iplt % ncol
  iplt += 1

  sttl = f'NSIDC IceConc avrg {YRS}-{YRE} {MMA:02d}'

  ax1 = axes[irow, icol]
  if iplt == 1:
    plt_clrb = True
  else:
    plt_clrb = False
  ax1 = plot_field(ax1, fig1, m, xR, yR, A2d, clrmp,rmin,rmax,plt_clrb,sttl=sttl)

btx = 'plot_NSIDC_clim_Nsbplts.py'
bottom_text(btx)



