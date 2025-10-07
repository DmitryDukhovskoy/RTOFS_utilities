"""
  Plot monthly ice thickness from seas f/cast experiments
  Specify months (calendar numbering!) to average statistics by seasons

  use only for NEP, as ARC does not need orthonormal projection

  use --help for more information on keywargs

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
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_sis2_relax as msisrlx
importlib.reload(mutob)
importlib.reload(manseas)

parser = argparse.ArgumentParser()
parser.add_argument("--expt", help="f/cast experiment number: 1, 2, 3, ..5, 11, 12,", type=int)
parser.add_argument("--yr", help="year to plot, default 2001 for NEP and 1995 for ARC", type=int)
#parser.add_argument("--regn", help="NEP or ARC", type=str, required=True)
parser.add_argument("--varnm", help="iconc or ithkn", type=str, required=True)
parser.add_argument("--mms", help="Calendar month to start averaging", type=int, required=True)
parser.add_argument("--mme", help="Calendar month to end averaging, default=MMS", type=int)
args = parser.parse_args()

# Test runs were performed for only 1 year
#regn = args.regn if args.regn else None
YRS = args.yr if args.yr else None
MMS = args.mms if args.mms else None
MME = args.mme if args.mme else MMS
varnm = args.varnm if args.varnm else None
expt_nmb = args.expt if args.expt else None

regn = 'NEP'

if YRS is None:
  if regn == 'NEP':
    YRS = 2001
  else:
    YRS = 1995

expt_name = f'{regn}phys_irlxtest_-expt{expt_nmb:02d}'

mstart = 1  # current test runs all started on Jan 1, 2001

# relax hours:
RLXH = [0,1,24,120,360]

if regn == 'NEP':
  pthtest = f'/archive/Dmitry.Dukhovskoy/fre/NEP/test_ice_relax/NEPphys_expt{expt_nmb:02d}/{YRS}-{mstart:02d}'
  pthrlx  = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax'
elif regn == 'ARC':
  pthtest = f'/archive/Dmitry.Dukhovskoy/fre/ARC12/test_ice_relax/ARCphys_expt{expt_nmb:02d}/{YRS}-{mstart:02d}'
  pthrlx  = '/work/Dmitry.Dukhovskoy/ARC12/irlx'


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP']['seasonal_daily']['pthgrid']
fgrid      = pthseas['MOM6_NEP']['seasonal_daily']['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"]['seasonal_daily']["ftopo"]
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
YRA = YRS
for MMA in range(MMS,MME+1):
  print(f"Processing {YRA}/{MMA} IceConc ...")

  dnmb0 = mtime.datenum([YRA, MMA, 15])

  dcice = os.path.join(pthtest,f'ice_month.nc')
  print(f'Reading {dcice}')

  # Note assumed start month = 1, if not - update imo 
  imo = MMA-1

  with xarray.open_dataset(dcice) as ds:
    C2d = ds['siconc'].isel(time=imo).data.squeeze()
    H2d = ds['sithick'].isel(time=imo).data.squeeze()   # ice thicknes m, need m3/m2 

  if varnm == 'iconc':
    A2d = C2d
  else:  
    A2d = H2d*C2d      # m3/m2 - grid cell mean thickness

  if icc == 0:
    AMN = A2d
  else:
    AMN = AMN + A2d

  icc += 1

if icc > 1:
  AMN  = AMN.squeeze()/icc

if varnm == 'iconc':
  clrmp = mclrmps.colormap_conc()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 1.
else:
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
  hcntrs = [4,6,8,10,12,14,16,18,20]

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

  btx = 'plot_fcst_iconc_mnthly.py'
  bottom_text(btx, pos=[0.2, 0.01])

  return ax1

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

if expt_nmb < 10:
  rlx_time = RLXH[expt_nmb-1]
else:
  enmb = expt_nmb // 10
  rlx_time = RLXH[enmb-1]

if MME == MMS:
  sttl = f'{expt_name} rlx={rlx_time}hr {varnm} {YRS}/{MMS:02d}'
else:
  sttl = f'{expt_name} rlx={rlx_time}hr {varnm} {YRS}/{MMS:02d}-{YRS}/{MME:02d}'


plt.ion()
fgnmb=1
ax1 = plot_field(1, m, xR, yR, A2d, clrmp,rmin,rmax,sttl=sttl)


