"""
  Plot monthly ice conce from NEP PHYS only spinup and hindcasts
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
importlib.reload(msisrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="Calendar year to start stat. averaging", type=int, required=True)
parser.add_argument("--YRE", help="Calendar year to end averaging of statistics", type=int)
parser.add_argument("--MMS", help="Calendar month to start averaging of statistics", type=int)
parser.add_argument("--MME", help="Calendar month to end averaging of statistics", type=int)
parser.add_argument("--varnm", help="iconc or ithkn", type=str, required=True)
args = parser.parse_args()

f_cntrobs = True   # Plot observation-derived ice edge

# experiments: 2 - daily OB seasonal forecasts, 3 - same as 2 but with sea ice relaxation
# Default values: 
expt   = 'hindcast_phys'   # spinup runs for NEP BGC hindcasts
# Default Averaging time period:
MMI   = 1    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YRS = 1993   # Start: f/cast init. year to use for monthly averaging
YRE = YRS   # End
MMS = 10
MME = MMS
expt_nmb = 3

YRS = args.YRS if args.YRS else None
YRE = args.YRE if args.YRE else YRS
varnm = args.varnm if args.varnm else None
if args.MMS:
  MMS = args.MMS
if args.MME:
  MME = args.MME
else:
  MME = MMS
 
# No ice relaxation, old hindcast with GLORYS nudging only 
expt = "hindcast"

YAVRG = [x for x in range(YRS,YRE+1)]

expt_name = 'NEPphys_nudged_' + expt
f_hcast = True
run_info = f'{expt_name} init MM={MMI}, ice conc: {min(YAVRG)}-{max(YAVRG)}'


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP']['seasonal_fcst']['pthgrid']
fgrid      = pthseas['MOM6_NEP']['seasonal_fcst']['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"]['seasonal_fcst']["ftopo"]
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

DX, DY = mmom6.dx_dy(hlon, hlat)
Acell  = DX*DY    # m2

TM = []
icc = 0
Nrec = 0
for YRA in (YAVRG):
  for MMA in range(MMS,MME+1):
    print(f"Processing {YRA}/{MMA} IceConc ...")
  
    dnmb0 = mtime.datenum([YRA, MMA, 15])
    #YRI = manseas.yr_init_fcst_from_datenum(dnmb0, MMI)  # init year
    YRI = YRA
    MMF = MMA

    if f_hcast:
    # assumed file naming for 3-month segment runs:
    # YYYYMM01.ice_month.nc
      if MMA < 4:
        MINIT=1
      elif MMA > 3 and MMA < 7:
        MINIT=4
      elif MMA > 6 and MMA < 10:
        MINIT=7
      elif MMA > 9:
        MINIT=10
    pthoutp = f"/archive/Dmitry.Dukhovskoy/fre/NEP/2024/NEP_physics_202404_nudging-15d/"+\
            f"gfdl.ncrc5-intel22-repro/history/{YRA}-{MINIT:02d}"

    dcice = os.path.join(pthoutp,'ice_month.nc')
    imo = MMF-MINIT  # consequtive months in the file from init. month

    print(f'Reading {dcice}')

    #MMF = manseas.mofcst_from_mocalend(YRI,MMI,MMA) # forecast month #

    ds = xarray.open_dataset(dcice)
    C2d = ds['siconc'].isel(time=imo).data.squeeze()
    H2d = ds['sithick'].isel(time=imo).data.squeeze()
    #V2d = ds['sivol'].isel(time=imo).data.squeeze()
    

    if varnm == 'iconc':
      A2d = C2d.copy()
    else:
      A2d = H2d*C2d  # Possible error in this approach for cell ice thickness estimate
      #A2d = msisrlx.cell_ithkn_sis2(dcice, imo)
      #A2d = H2d.copy()

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



# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

mo_fcast = msisrlx.cal_mo_to_fcast(MMI,MMS)
mE_fcast = msisrlx.cal_mo_to_fcast(MMI,MME)
sttl = (f'{expt_name} {varnm} {YRS} mo={MMS}')

plt.ion()
fgnmb=1
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)

if varnm == 'ithkn':
  cntr_clr = [0.95, 0.95, 1.0]
  clevel = [x for x in range(5,30,5)]
  CS = ax1.contour(xR, yR, A2d, clevel, linestyles='solid', colors=[cntr_clr], linewidths=1)
  ax1.clabel(CS, inline=1, fontsize=12)


if varnm == 'iconc' and f_cntrobs:
  YR0 = YRA
  MM0 = MMA
  ci0 = 0.15
  #cntr_clr = [1, 0.4, 0]
  #ax1.contour(xR,yR,A2d,[ci0],linestyles='solid', colors=[cntr_clr], linewidths=2)
  # Use interpolated fields
  pthnsidc = f'/work/Dmitry.Dukhovskoy/data/NRT_NOAA_NSIDC_seaconc/{YR0}_mnth'
  fliceout = f'NSIDC_iconc_mnth_interpNEP816x342_{YR0}.nc'
  dfliceout = os.path.join(pthnsidc,fliceout)
  dset_nsidc = xarray.open_dataset(dfliceout)
  imo = MM0-1
  ICnrt = dset_nsidc['ice_conc'].isel(time=imo).data

  cntr_nsidc = [1.,0.3,0.]
  ax1.contour(xR,yR,ICnrt,[ci0],linestyles='solid', colors=[cntr_nsidc], linewidths=2)
  

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

btx = 'plot_PHYShndcst_iconc_ithkn_mnthly.py'
bottom_text(btx, pos=[0.2, 0.01])



