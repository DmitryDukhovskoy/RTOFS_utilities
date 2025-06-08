"""
  Plot monthly ice thickness from NEP BGC spinup and hindcasts
  Specify months (calendar numbering!) to average statistics by seasons

  Usage: plot_BGChndcst_iconc_mnthly.py --MMI=4 --YRS=1993 --YRE=1993 --MMS=10 --MME=10
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
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_sis2_relax as msisrlx
importlib.reload(mutob)
importlib.reload(manseas)

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="Calendar (! not init.) year to start stat. averaging", type=int)
parser.add_argument("--YRE", help="Calendar year to end averaging of statistics", type=int)
parser.add_argument("--MMS", help="Calendar month to start averaging of statistics", type=int)
parser.add_argument("--MME", help="Calendar month to end averaging of statistics", type=int)
args = parser.parse_args()

f_cntrobs = True   # Plot observation-derived ice edge

# experiments: 2 - daily OB seasonal forecasts, 3 - same as 2 but with sea ice relaxation
# Default values: 
expt   = 'spinup_bgc'   # spinup runs for NEP BGC hindcasts
# Default Averaging time period:
MMI   = 1    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YRS = 1993   # Start: f/cast init. year to use for monthly averaging
YRE = YRS   # End
MMS = 10
MME = MMS
expt_nmb = 3


if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
else:
  YRE=YRS
if args.MMS:
  MMS = args.MMS
if args.MME:
  MME = args.MME
else:
  MME = MMS

YAVRG = [x for x in range(YRS,YRE+1)]

expt_name = f'NEPbgc_nudged_spinup'
run_info = f'{expt_name} init MM={MMI}, ice conc: {min(YAVRG)}-{max(YAVRG)}'


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
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
for YRA in (YAVRG):
  pthoutp = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YR=YRA)
  for MMA in range(MMS,MME+1):
    print(f"Processing {YRA}/{MMA} IceConc ...")

    dnmb0 = mtime.datenum([YRA, MMA, 15])
    #YRI = manseas.yr_init_fcst_from_datenum(dnmb0, MMI)  # init year
    YRI = YRA

    dcice = os.path.join(pthoutp,f'ice_month.nc')
    print(f'Reading {dcice}')

    #MMF = manseas.mofcst_from_mocalend(YRI,MMI,MMA) # forecast month #
    MMF = MMA
    imo = MMF-1

    ds = xarray.open_dataset(dcice)
    C2d = ds['siconc'].isel(time=imo).data.squeeze()
    H2d = ds['sithick'].isel(time=imo).data.squeeze()

    A2d = H2d*C2d
    if icc == 0:
      AMN = A2d
    else:
      AMN = AMN + A2d

    icc += 1

if icc > 1:
  AMN  = AMN.squeeze()/icc

# Mask out deep region:
#AMN1 = np.where(HH<-500, 1.e3, AMN1)
#AMN2 = np.where(HH<-500, 1.e3, AMN2)

#if varnm == 'iconc':
#clrmp = mclrmps.colormap_conc()
clrmp = mclrmps.colormap_ice_thkn()
clrmp.set_bad(color=[0.2, 0.2, 0.2])
rmin = 0.
rmax = 4.

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

mo_fcast = msisrlx.cal_mo_to_fcast(MMI,MMS)
mE_fcast = msisrlx.cal_mo_to_fcast(MMI,MME)
sttl = (f'{expt_name} IceThkn {YRS} mo={MMS}')

plt.ion()
fgnmb=1
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(sttl)

cntr_clr = [0.2, 0.7, 1.0]
if fgnmb==1 : 
  clevel=0.
else:
  clevel=0.15

#CS = ax1.contour(xR, yR, A2d, [clevel], linestyles='solid', colors=[cntr_clr], linewidths=1)
#ax1.clabel(CS, inline=1, fontsize=10)
  
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

btx = 'plot_BGChndcst_ithkn_mnthly.py'
bottom_text(btx, pos=[0.2, 0.01])



