"""
  Plot monthly ice thkn from NEP BGC hindcasts with ocean GLORYS and ice PIOMAS nudging
  Specify months (calendar numbering!) to average statistics by seasons

  use --help for more information on keywargs

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import pickle
from copy import copy
import matplotlib.colors as colors
from yaml import safe_load
import argparse
import pickle

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
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_sis2_relax as msisrlx
importlib.reload(mutob)
importlib.reload(manseas)

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="Calendar year to start stat. averaging", type=int, required=True)
parser.add_argument("--YRE", help="Calendar year to end averaging of statistics", type=int)
parser.add_argument("--MMS", help="Cal. month to start averaging of statistics", type=int, required=True)
parser.add_argument("--MME", help="Cal. month to end averaging of statistics", type=int)
args = parser.parse_args()

f_cntrobs = True   # Plot observation-derived ice edge

# experiments: 2 - daily OB seasonal forecasts, 3 - same as 2 but with sea ice relaxation
# Default values: 
expt_nmb  = 2
trlx = 12   # max relax. time scale, hrs
expt_name = f'NEPbgc_nudged_hindcast{expt_nmb:02d}'
hcst_time = 3 # f/csat time interval, months
hcst_interv = np.array([x for x in range(1,12+hcst_time,hcst_time)], dtype=int)

YRS = args.YRS if args.YRS else None
YRE = args.YRE if args.YRE else YRS
MMS = args.MMS if args.MMS else None
MME = args.MME if args.MME else MMS

YAVRG = [x for x in range(YRS,YRE+1)]
run_info = f'{expt_name}  ice conc: {min(YAVRG)}-{max(YAVRG)} months={MMS}-{MME}'

fyaml = 'bgc_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP']['hindcast']['pthgrid']
fgrid      = pthseas['MOM6_NEP']['hindcast']['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"]['hindcast']["ftopo"]
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

TM = []
icc = 0
Nrec = 0
for YRA in YAVRG:
  for MMA in range(MMS,MME+1):
    print(f"Processing {YRA}/{MMA} IceThkn ...")

    # Find init date for given month, assuming hcst_time (n months) f/cast interval
    kint = np.searchsorted(hcst_interv, MMA, side='right') - 1
    assert(hcst_interv[kint] <= MMA < hcst_interv[kint+1]), f'Wrong time bin {kint} for {MMA}'
    MINIT = hcst_interv[kint]
    imo = MMA-MINIT      # current month in the archive output

    pthhcst = (
         pthseas['MOM6_NEP']['hindcast']['pthoutp'].format(
           hindcast_name=expt_name, YR=YRA, MM=MINIT, DD=1
      )
    )
    dnmb0 = mtime.datenum([YRA, MMA, 15])

    dcice = os.path.join(pthhcst,f'ice_month.nc')
    print(f'Reading {dcice}')

    ds = xarray.open_dataset(dcice)
    C2d = ds['siconc'].isel(time=imo).data.squeeze()
    H2d = ds['sithick'].isel(time=imo).data.squeeze()   # ice thicknes m, need m3/m2 
    A2d = H2d*C2d      # m3/m2 - grid cell mean thickness

    if icc == 0:
      AMN = A2d
    else:
      AMN = AMN + A2d

    icc += 1

if icc > 1:
  AMN  = AMN.squeeze()/icc

clrmp = mclrmps.colormap_ice_thkn()
clrmp.set_bad(color=[0.2, 0.2, 0.2])
rmin = 0.
rmax = 4.


def plot_field(fgnmb, m, xR, yR, A2d, clrmp, rmin, rmax, sttl=[]):
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(sttl)

  cntr_clr = [0.95, 0.95, 0.95]
  hcntrs = [6,8,10,12,14,16,18,20]

  CS = ax1.contour(xR, yR, A2d, hcntrs, linestyles='solid', colors=[cntr_clr], linewidths=1)
  ax1.clabel(CS, inline=1, fontsize=10)
    
  # extend: min, max, both
  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  if fgnmb>1:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')
  else:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')

  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'plot_fcst_ithkn_mnthly.py'
  bottom_text(btx, pos=[0.2, 0.01])

  return ax1

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

if MMS==MME and YRS==YRE:
  sttl = f'{expt_name} IceThkn {YRS} M={MMS}, Trlx_max={trlx}hrs'
else:
  sttl = (f'{expt_name} IceThkn Trlx_max={trlx}hrs\n' + \
         f'avrg {YRS}-{YRE} M={MMS}-{MME}')

plt.ion()
fgnmb=1
ax1 = plot_field(1, m, xR, yR, A2d, clrmp,rmin,rmax,sttl=sttl)


