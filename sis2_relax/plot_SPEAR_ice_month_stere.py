"""
  Plot monthly SPEAR ice fields
 
  Extract and subsample for NEP domain using bash script:
  ./subset_spear_ice.sh 2010 2010 1 1

  Note that year 2010 and month 1 designate the initialization time
  SPEAR monthly fields have 12 f/cast months in each file


  usage: plot_SPEAR_ice_month_stere.py --varnm={ithkn,iarea or iconc} --YRI=1993 --MMI=1 --mo=8

  Only 12 months of the f/cast are saved in the SPEAR files
  mo - calendar month, depending on the init month, it may be at the end/start of the f/cast

"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
from yaml import safe_load
import argparse

PPTHN = '/Users/ddmitry/python'
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc
from mod_utils_fig import bottom_text
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--YRI", help="init year of SPEAR f/cast: 1993, ..., 2020", type=int)
parser.add_argument("--MMI", help="init month of SPEAR f/cast: 1, ..., 12", type=int)
parser.add_argument("--mo", help="cal. month to plot: 1,..., 12, ...", type=int)
parser.add_argument("--ensmb", help="ensemble number, 1,..., 15", type=int)
parser.add_argument("--varnm", help="field to plot: ithkn or iarea", type=str)
args = parser.parse_args()

f_cntrobs = True   # Plot observation-derived ice edge
# Years in the relax file also used in the rlx file name:
YRI = 1993  # init yr
MMI = 1     # init month
MM0 = MMI   # month to plot
ifld = 'iarea'  # ithkn, iarea
ens_nmb = 1  # SPEAR ensemble #

if args.YRI:
  YRI = args.YRI
if args.MMI:
  MMI = args.MMI
  MM0 = MMI
if args.mo:
  MM0 = args.mo
if args.varnm:
  ifld = args.varnm
  if ifld=='iconc' or ifld=='iarea':
    ifld = 'siconc'
  elif ifld=='ithkn' or ifld=='ithk':
    ifld = 'sithick'
if args.ensmb:
  ens_nmb = args.ensmb

varnm = ifld
pthdata = f'/work/Dmitry.Dukhovskoy/tmp/spear_subset/{YRI}/ens{ens_nmb:02d}'

# Read saved relax. fields:
flout = f'NEP_spear_{YRI}{MMI:02d}.{ifld}.nc'
dfspear = os.path.join(pthdata,flout)
ds_spear = xarray.open_dataset(dfspear)
# Assumed that rec #1 = init month
#Time = ds_spear['time'].data
#TM = mmisc.convert_nptime_to_datenum(Time)
#dnmb0 = mtime.datenum([YR0,MM0,15,12])
mcal = np.arange(MMI,MMI+12)
mcal = np.where(mcal>12, mcal-12, mcal)
D = abs(mcal-MM0)
itime = np.argmin(D)

YR0=YRI
if MM0 < MMI:
  YR0=YRI+1

A2d = ds_spear[varnm].isel(time=itime).data
LON = ds_spear['GEOLON'].data
LAT = ds_spear['GEOLAT'].data

match ifld:
  case('sithick'):
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 4.
  case('siconc'):
    clrmp = mclrmps.colormap_conc()
    rmin = 0.
    rmax = 1.

clrmp.set_bad(color=[0.2, 0.2, 0.2])


def plot_ice(fgnmb, xR, yR, A2d, clrmp, rmin, rmax, sttl, xTst=-1, yTst=-1):
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = ax1.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  #  img = ax1.pcolormesh(RLXHR, cmap=clrmp)

  if xTst >=0 and yTst >= 0:
    ax1.plot(xTst,yTst,'o')

  ax1.set_title(sttl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  # extend: min, max, both
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'plot_SPEAR_ice_month_stere.py'
  bottom_text(btx, pos=[0.2, 0.01])

  return ax1

plt.ion()

#A2d = np.where(HH>=0, np.nan, A2d)

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
#m = Basemap(width=5000*1.e3,height=5000*1.e3, resolution='l',\
#            projection='stere', lat_ts=50, lat_0=62, lon_0=-165)
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
          projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(LON, LAT)

fgnmb=1
mo_fcast = msisrlx.cal_mo_to_fcast(MMI,MM0)
sttlS = f'SPEAR init: {YRI}/{MMI:02d}-e{ens_nmb:02d} Lead: {mo_fcast}, {ifld} {YR0}/{MM0}'
ax1 = plot_ice(fgnmb, xR, yR, A2d, clrmp, rmin, rmax, sttlS)

if f_cntrobs:
  ci0=0.15
  cntr_clr = [1, 0.4, 0]
  ax1.contour(xR,yR,A2d,[ci0],linestyles='solid', colors=[cntr_clr], linewidths=2)
  # Use interpolated fields
  pthnsidc = f'/work/Dmitry.Dukhovskoy/data/NRT_NOAA_NSIDC_seaconc/{YR0}_mnth'
  fliceout = f'NSIDC_iconc_mnth_interpNEP816x342_{YR0}.nc'
  dfliceout = os.path.join(pthnsidc,fliceout)
  dset_nsidc = xarray.open_dataset(dfliceout)
  imo = MM0-1
  ICnrt = dset_nsidc['ice_conc'].isel(time=imo).data

  # Hgrid lon. lat:
  fyaml = 'paths_seasfcst.yaml'
  with open(fyaml) as ff:
    pthseas = safe_load(ff)

  expt = 'seasonal_daily'
  pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
  fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
  ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
  hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
  hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
  dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
  dfgrid_mom = os.path.join(pthtopo, fgrid)
  hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

  xRm, yRm = m(hlon,hlat)

  cntr_nsidc = [0.,0.2,0.6]
  ax1.contour(xRm,yRm,ICnrt,[ci0],linestyles='solid', colors=[cntr_nsidc], linewidths=1)
  
  ax3 = plt.axes([0.02, 0.12, 0.1, 0.1])
  dx = 0.15
  ax3.plot([0,dx],[0.1,0.1],'-',color=cntr_clr)
  ax3.text(2*dx, 0.1, 'SPEAR', va='center')
  ax3.plot([0,dx],[0.2,0.2],'-',color=cntr_nsidc)
  ax3.text(2*dx,0.2, 'NSIDC NRT', va='center')
  ax3.set_xlim([0,0.7])
  ax3.set_ylim([0,0.3])
  ax3.axis('off')







