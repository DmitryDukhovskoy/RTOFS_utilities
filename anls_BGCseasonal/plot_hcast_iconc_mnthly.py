"""
  Plot monthly ice conce from NEP BGC hindcasts with ocean GLORYS and ice PIOMAS nudging 
  Specify months (calendar numbering!) to average statistics by seasons

  Usage: plot_fcst_iconc_mnthly.py --MMI=4 --YRS=1993 --YRE=1993 --MMS=10 --MME=10
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
importlib.reload(mutob)
importlib.reload(manseas)

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="Calendar year to start stat. averaging", type=int, required=True)
parser.add_argument("--YRE", help="Calendar year to end averaging of statistics", type=int)
parser.add_argument("--MMS", help="Cal. month to start averaging of statistics", type=int, required=True)
parser.add_argument("--MME", help="Cal. month to end averaging of statistics", type=int)
args = parser.parse_args()

f_cntrobs = True   # Plot observation-derived ice edge
f_cntrrlx = True   # Plot ice edge from the ice concentration target fields

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
    print(f"Processing {YRA}/{MMA} IceConc ...")

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
    A2d = ds['siconc'].isel(time=imo).data.squeeze()

    if icc == 0:
      AMN = A2d.copy()
    else:
      AMN = AMN + A2d

    icc += 1

if icc > 1:
  AMN  = AMN.squeeze()/icc

A2d = AMN.copy()

#if varnm == 'iconc':
clrmp = mclrmps.colormap_conc()
clrmp.set_bad(color=[0.2, 0.2, 0.2])
rmin = 0.
rmax = 1.


def plot_field(fgnmb, m, xR, yR, A2d, clrmp, rmin, rmax, sttl=[]):
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

  btx = 'plot_hcst_iconc_mnthly.py'
  bottom_text(btx, pos=[0.2, 0.01])

  return ax1

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)


if MMS==MME and YRS==YRE:
  sttl = f'{expt_name} IceConc {YRS} M={MMS}, Trlx_max={trlx}hrs'
else:
  sttl = (f'{expt_name} IceConc Trlx_max={trlx}hrs\n' + \
         f'avrg {YRS}-{YRE} M={MMS}-{MME}')

plt.ion()
fgnmb=1
ax1 = plot_field(1, m, xR, yR, A2d, clrmp,rmin,rmax,sttl=sttl)

if f_cntrobs:
  # Not averaged, for 1 year/month
  # Change if averaged contour is needed
  ci0=0.15
  cntr_clr = [1, 0.4, 0]
  ax1.contour(xR,yR,A2d,[ci0],linestyles='solid', colors=[cntr_clr], linewidths=1.2)
  # Use interpolated fields
  kcc = 0
  for YR0 in YAVRG:
    for MM0 in range(MMS,MME+1):
      pthnsidc = f'/work/Dmitry.Dukhovskoy/data/NRT_NOAA_NSIDC_seaconc/{YR0}_mnth'
      fliceout = f'NSIDC_iconc_mnth_interpNEP816x342_{YR0}.nc'
      dfliceout = os.path.join(pthnsidc,fliceout)
      dset_nsidc = xarray.open_dataset(dfliceout)
      imo = MM0-1
      Itmp = dset_nsidc['ice_conc'].isel(time=imo).values

      if kcc == 0:
        Isum = Itmp.copy()
      else:
        Isum = Isum + Itmp
      kcc += 1
  
  if kcc == 1:
    ICnrt = Isum
  else:
    ICnrt = Isum / float(kcc)

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
  ax1.contour(xRm,yRm,ICnrt,[ci0],linestyles='solid', colors=[cntr_nsidc], linewidths=1.5)

if f_cntrrlx:
  # Target fields:
  # Read saved relax. fields:
  ifld = 'iarea'
  pthsis = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax/'
  kcc = 0
  for YR0 in YAVRG:
    YR1 = YR0
    YR2 = YR1+1
    flout = f'PIOMASv21_ithkn_iconc_{YR1}_{YR2}_monthly.nc'
    diclim = os.path.join(pthsis, flout)
    print(f'Reading relax fields from {diclim}')
    ds_rlx = xarray.open_dataset(diclim)
    Time = ds_rlx['time'].data
    TM = mmisc.convert_nptime_to_datenum(Time)
    dnmb0 = mtime.datenum([YR0,MM0,15,12])
    D = abs(TM-dnmb0)
    itime = np.argmin(D)
    dv0 = mtime.datevec(TM[itime])
    assert dv0[0]==YR0, f'Requested YR={YR0}, year in rlx file={dv0[0]}'
    assert dv0[1]==MM0, f'Requested month={MM0}, month in rlx file={dv0[1]}'

    Itmp = ds_rlx[ifld].isel(time=itime).data
    if kcc == 1:
      Isum = Itmp
    else:
      Isum = Isum + Itmp
    kcc += 1

  if kcc == 1:
    ICrlx = Isum
  else:
    ICrlx = Isum / float(kcc)

  # Use grid from the hindcast - same grid
  cntr_irlx = [1,0,1]
  ax1.contour(xR,yR,ICrlx,[ci0],linestyles='solid', colors=[cntr_irlx], linewidths=1.5)


ax3 = plt.axes([0.02, 0.12, 0.1, 0.1])
dx = 0.15
dy = 0.08
yy0 = 0.1
ax3.plot([0,dx],[yy0,yy0],'-',color=cntr_clr)
ax3.text(2*dx, 0.1, 'F/cast', va='center')
if f_cntrobs:
  yy0 = yy0 + dy
  ax3.plot([0,dx],[yy0,yy0],'-',color=cntr_nsidc)
  ax3.text(2*dx,yy0, 'NSIDC NRT', va='center')
if f_cntrrlx:
  yy0 = yy0 + dy
  ax3.plot([0,dx],[yy0,yy0],'-',color=cntr_irlx)
  ax3.text(2*dx,yy0, 'PIOMAS target', va='center')
  
ax3.set_xlim([0,0.7])
ax3.set_ylim([0,0.3])
ax3.axis('off')


