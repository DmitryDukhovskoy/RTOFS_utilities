"""
  Create N x M subplots to show N*M months

  Plot monthly ice conc/thickness from NEP BGC hindcasts with ocean GLORYS and ice PIOMAS nudging 
  Specify months (calendar numbering!) to average statistics by seasons

  Usage: plot_fcst_iconc_mnthly.py --MMI=4 --YRS=1993 --YRE=1993 --MFS=10 --MFE=10
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

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="Calendar year to start stat. averaging", type=int, required=True)
parser.add_argument("--YRE", help="Calendar year to end averaging of statistics", type=int)
parser.add_argument("--MMI", help="Initialization month=1,4,7,10", type=int, required=True)
parser.add_argument("--ens", help="Ensemble number=1,...,10", type=int, required=True)
parser.add_argument("--MFS", help="Forecast month start to plot: 1,...,12", type=int, required=True)
parser.add_argument("--MFE", help="Forecast month end to plot: 1,...,12", type=int)
parser.add_argument("--ncol", help="N of columns for subplots", type=int)
parser.add_argument("--nrow", help="N of rows for subplots", type=int)
parser.add_argument("--varnm", help="iconc or ithkn", type=str, required=True)
args = parser.parse_args()

f_cntrobs = True   # Plot observation-derived ice edge, for iconc only
f_cntrrlx = True   # Plot ice edge from the ice concentration target fields

# Default values: 
expt_nmb  = 1
trlx = 24   # max relax. time scale, hrs
expt_name = f'NEPbgc_fcst_dailyOB{expt_nmb:02d}'  

YRS     = args.YRS if args.YRS else None
YRE     = args.YRE if args.YRE else YRS
MMI     = args.MMI if args.MMI else None
ens_nmb = args.ens if args.ens else None
MFS     = args.MFS if args.MFS else None       # Forecast (i.e. lead time) month
MFE     = args.MFE if args.MFE else MFS
varnm   = args.varnm if args.varnm else None
nmnths  = MFE-MFS+1
ncol    = args.ncol if args.ncol else None
nrow    = args.nrow if args.nrow else nmnths

if not ncol:
  ncol = min([nmnths,4])
  nrow = nmnths // ncol

assert ncol*nrow == nmnths, f'Specified columns/rows ({ncol}/{nrow}) do not match N months {nmnths}'

YAVRG = [x for x in range(YRS,YRE+1)]
run_info = f'{expt_name}  ice conc: {min(YAVRG)}-{max(YAVRG)} '

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

if varnm == 'iconc':
  clrmp = mclrmps.colormap_conc()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 1.

  ci0=0.15
  hcntrs = [cio]
  cntr_clr = [1, 0.4, 0]

elif varnm == 'ithkn':
  clrmp = mclrmps.colormap_ice_thkn()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 4.

  hcntrs = [x for x in range(1,10)]
  cntr_clr = [.95, 0.95, 0.95]

  
#cntr_clr = [1, 0.4, 0]
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
    #ax2 = fig1.add_axes([0.93,0.1,0.012,0.8])
    ax2 = fig1.add_axes([0.9,0.1,0.013,0.8])
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')

    ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
    ax2.set_yticklabels(ax2.get_yticks())
    ticklabs = clb.ax.get_yticklabels()
    #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
    clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
    clb.ax.tick_params(direction='in', length=12)

  return ax1

def contour_nsidc(ax1,YAVRG,MM0,xR,yR, cntr_nsidc):
  ci0=0.15
  # Use interpolated fields
  kcc = 0
  for YR0 in YAVRG:
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
  dfgrid_mom = os.path.join(pthtopo, fgrid)
  hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

  xRm, yRm = m(hlon,hlat)

  ax1.contour(xRm,yRm,ICnrt,[ci0],linestyles='solid', colors=[cntr_nsidc], linewidths=1.5)

  return ax1

def contour_piomas_irlx(ax1,xR,yR,YAVRG,MM0,cntr_irlx,HH):
  # Target fields:
  # Read saved relax. fields:
  # For BGC f/casts: 5yr average fields are used fpr target fields
  ci0=0.15
  ifld = 'iarea'
  pthsis = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax/'
  kcc = 0
  for YR0 in YAVRG:
    YR1 = YR0
    YR2 = YR1+1
    flout = f'PIOMASv21_ithkn_iconc_{YR1}_{YR2}_avrg5yr.nc'
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
    kcc += 1
    if kcc == 1:
      Isum = Itmp
    else:
      Isum = Isum + Itmp

  if kcc == 1:
    ICrlx = Isum
  else:
    ICrlx = Isum / float(kcc)

  ICrlx = np.where(HH>=0, np.nan, ICrlx)
  # Use grid from the hindcast - same grid
  ax1.contour(xR,yR,ICrlx,[ci0],linestyles='solid', colors=[cntr_irlx], linewidths=1.5)

  return ax1

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

plt.ion()
fig1 = plt.figure(1,figsize=(15, 10))
fig1.clf()  # Clear the figure
axes = fig1.subplots(nrows=nrow, ncols=ncol, squeeze=False)

# Shift subplots to the left
# and up for colorbar and text
# also keep subplots close to each other : wspace, hspace
fig1.subplots_adjust(
    left=0.05,
    right=0.85,  # More right-side room
    top=0.95,
    bottom=0.1,
    wspace=0.05,
    hspace=0.1
)

iplt = 0
# Loop over forecast months
for MFA in range(MFS,MFE+1):
  icc = 0
  for YRA in YAVRG:
    print(f"Processing {YRA}/{MFA} IceConc ...")

    # Calendar yr/months:
    TCAL = manseas.yrmo_seasonal_fcst(YRS,MMI)

    # Find init date for given month, assuming hcst_time (n months) f/cast interval
    #kint = np.searchsorted(hcst_interv, MFA, side='right') - 1        
    imo = MFA-1      # current month in the archive output

    pthfcst = '/archive/Dmitry.Dukhovskoy/fre/NEP/forecast_bgc/NEPbgc_fcst_dailyOB01/' +\
              f'{YRA}-{MMI:02d}-e{ens_nmb:02d}/history'

    dnmb0 = mtime.datenum([YRA, MFA, 15])

    dcice = os.path.join(pthfcst,f'ice_month.nc')
    print(f'Reading {dcice}')

    ds = xarray.open_dataset(dcice)
    Cice = ds['siconc'].isel(time=imo).data.squeeze()
    if varnm == 'ithkn':
      Hice = ds['sithick'].isel(time=imo).data.squeeze()
      A2d = Cice*Hice
    else:
      A2d = Cice

    if icc == 0:
      AMN = A2d.copy()
    else:
      AMN = AMN + A2d

    icc += 1

  if icc > 1:
    AMN  = AMN.squeeze()/icc

  A2d = AMN.copy()

  irow = iplt // ncol 
  icol = iplt % ncol  
  iplt += 1

  mcal = TCAL[imo,1]
  ycal = TCAL[imo,0]
  if YRS==YRE:
    sttl = f'{varnm} {ycal}/{mcal:02d} {YRS}-{MMI:02d}-e{ens_nmb:02d} '
  else:
    sttl = (f'{varnm} avrg {YRS}-{YRE}/{mcal:02d}, {MMI:02d}-e{ens_nmb:02d}')

  ax1 = axes[irow, icol]
  if iplt == 1:
    plt_clrb = True
  else:
    plt_clrb = False
  ax1 = plot_field(ax1, fig1, m, xR, yR, A2d, clrmp,rmin,rmax,plt_clrb,sttl=sttl)

  if varnm == 'iconc':
    ax1.contour(xR,yR,A2d,hcntrs,linestyles='solid', colors=[cntr_clr], linewidths=1.2)
  elif varnm == 'ithkn':
    CS1 = ax1.contour(xR,yR,A2d,hcntrs,linestyles='solid', colors=[cntr_clr], linewidths=1.2)
    ax1.clabel(CS1, inline=1, fontsize=10)

  if f_cntrobs and varnm == 'iconc':
    ax1 = contour_nsidc(ax1,YAVRG,MFA,xR,yR,cntr_nsidc)
  if f_cntrrlx and varnm == 'iconc':
    ax1 = contour_piomas_irlx(ax1,xR,yR,YAVRG,MFA,cntr_irlx,HH)

#plt.tight_layout()

if varnm == 'iconc':
  ax3 = plt.axes([0.01, 0.01, 0.08, 0.1])
  dx = 0.15
  dy = 0.05
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
    ax3.text(2*dx,yy0, 'PIOMAS trgt', va='center')
    
  ax3.set_xlim([0,0.9])
  ax3.set_ylim([0,0.4])
  ax3.axis('off')

sinfo = f'Monthly {varnm} from BGC seas. forecasts \n'
sinfo = sinfo + f'init: {YRS}/{MMI} ens={ens_nmb:02d}\n'
sinfo = sinfo + dcice

ax4 = plt.axes([0.3, 0.045, 0.6,0.03])
ax4.text(0, 0, sinfo, fontsize=8)
ax4.axis('off')

btx = 'plot_fcast_ice_mnthly_Nsubplots.py'
bottom_text(btx, pos=[0.2, 0.01])



