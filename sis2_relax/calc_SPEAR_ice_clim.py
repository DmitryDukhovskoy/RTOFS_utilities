"""
  Compute monthly SPEAR ice climatologies for ice conc and thickn. 
 
  Extract and subsample for NEP domain using bash script:
  ./subset_spear_ice.sh 2010 2010 1 1

  Note that year 2010 and month 1 designate the initialization time
  SPEAR monthly fields have 12 f/cast months in each file


  usage: 
  calc_SPEAR_ice_clim.py --varnm iconc --YRS 2010 --YRE 2020 --MMI=1 --ensmb 1

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
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc
from mod_utils_fig import bottom_text
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="start year for deriving climat.: 1993, ..., 2022", type=int)
parser.add_argument("--YRE", help="end year for deriving climat.: 1993, ..., 2022", type=int)
parser.add_argument("--MMI", help="init month of SPEAR f/cast: 1, ..., 12", type=int)
parser.add_argument("--ensmb", help="ensemble number, 1,..., 15", type=int)
parser.add_argument("--varnm", help="field: ithkn or iarea", type=str)
args = parser.parse_args()

f_save = True
# Years in the relax file also used in the rlx file name:
YRS = 2010  # init yr
YRE = 2015
MMI = 1     # init month
ifld = 'iarea'  # ithkn, iarea
ens_nmb = 1  # SPEAR ensemble #

if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI
if args.varnm:
  ifld = args.varnm
  if ifld=='iconc' or ifld=='iarea':
    ifld = 'siconc'
  elif ifld=='ithkn' or ifld=='ithk':
    ifld = 'sithick'
if args.ensmb:
  ens_nmb = args.ensmb

varnm = ifld

if not f_save:
  print(f'WARNING: fields wont be saved, f_save flag is off !!!\n')

# Read SPEAR fields:
# Note that SPEAR thickness is not "volume/unit area" (like in CICE output or PIOMAS)
# and needs to be multiplied by partial area
icc = 0
for YRI in range(YRS,YRE+1):
  print(f'Processing {YRI}')
  pthdata = f'/work/Dmitry.Dukhovskoy/tmp/spear_subset/{YRI}/ens{ens_nmb:02d}'
  flthkn = f'NEP_spear_{YRI}{MMI:02d}.sithick.nc'
  dfthkn = os.path.join(pthdata,flthkn)
  dspear_thkn = xarray.open_dataset(dfthkn)
  flconc = f'NEP_spear_{YRI}{MMI:02d}.siconc.nc'
  dfconc = os.path.join(pthdata,flconc)
  dspear_conc = xarray.open_dataset(dfconc)
  # Assumed that rec #1 = init month
  #Time = ds_spear['time'].data
  #TM = mmisc.convert_nptime_to_datenum(Time)
  #dnmb0 = mtime.datenum([YR0,MM0,15,12])
  mcal = np.arange(MMI,MMI+12)
  mcal = np.where(mcal>12, mcal-12, mcal)
  #D = abs(mcal-MM0)
  #itime = np.argmin(D)

  if icc == 0:
    xh = dspear_thkn['xh'].data
    yh = dspear_thkn['yh'].data
    xdim = len(xh)
    ydim = len(yh)
    ASUM = np.zeros((12,ydim,xdim))
    
  # Process by f/casts months
  for itime in range(12):
    H2d = dspear_thkn['sithick'].isel(time=itime).data
    C2d = dspear_conc['siconc'].isel(time=itime).data
    if ifld == 'siconc':
      A2d = C2d.copy()
    elif ifld == 'sithick':
      A2d = H2d*C2d         # ice m ---> m3/m2 

    hmax = np.nanmax(A2d)
    hmin = np.nanmin(A2d)
    print(f'N={itime+1} M={mcal[itime]} {ifld} min/max: {hmin:.2f}/{hmax:.2f}')

    ASUM[itime,:,:] = ASUM[itime,:,:] + A2d
    if icc == 0:
      LON = dspear_thkn['GEOLON'].data
      LAT = dspear_thkn['GEOLAT'].data

  icc += 1

ASUM = ASUM / icc

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


def plot_ice(fgnmb, m, xR, yR, A2d, clrmp, rmin, rmax, sttl):
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = ax1.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  #  img = ax1.pcolormesh(RLXHR, cmap=clrmp)

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

  btx = 'calc_SPEAR_ice_clim.py'
  bottom_text(btx, pos=[0.2, 0.01])


if f_save:
  pthpkl = '/work/Dmitry.Dukhovskoy/anls_output/spear_ice'
  floutp = f'spear_{varnm}_clim_{YRS}_{YRE}_MI{MMI:02d}.pkl'
  dflout = os.path.join(pthpkl,floutp)
  print(f'Dumping climtology --> {dflout}')
  with open(dflout,'wb') as fid:
    pickle.dump([ASUM,LON,LAT],fid) 

# Check
f_check = False

if f_check:
  plt.ion()

  MM0 = 12  # month to check
  D = abs(mcal-MM0)
  it0 = np.argmin(D)
  A2d = ASUM[it0,:,:].squeeze()
  #A2d = np.where(HH>=0, np.nan, A2d)

  # Stereographic Map projection:
  from mpl_toolkits.basemap import Basemap, cm
  #m = Basemap(width=5000*1.e3,height=5000*1.e3, resolution='l',\
  #            projection='stere', lat_ts=50, lat_0=62, lon_0=-165)
  m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

  xR, yR = m(LON, LAT)

  fgnmb=1
  sttlS = f'SPEAR clim init: {YRI}/{MMI:02d}-e{ens_nmb:02d}, {ifld} {YRS}-{YRE} {MM0}'
  plot_ice(fgnmb, m, xR, yR, A2d, clrmp, rmin, rmax, sttlS)

