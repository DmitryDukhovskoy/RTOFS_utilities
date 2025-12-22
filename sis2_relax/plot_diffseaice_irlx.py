"""

  Plot difference of daily ice fields (ice_day.nc) or monthly ice_month 
  from test experiments with ice relaxation irlx
  vs PIOMAS monthly
 
  May be:
  Use day 15 for a given month to compare with the monthly PIOMAS 
  since this should be the closest 


 some experiments were run with no ice ridging:
 no irlx: expt 11 - no ridging, 01 - with ridging, 12 - with ridging different advection schemes
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
importlib.reload(msisrlx)

#default values:
nexp0 = 1
regn  = 'NEP'
fday  = 0    # =0 - monthly data, >0 - plot daily, day=fday
imnth = 1    # =1 interp monthly PIOMAS rlx filds to daily for accuracy

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="year to plot, default NEP=2001, ARC=1993", type=int)
parser.add_argument("--regn", help=f"Model domain: NEP or ARC, default = {regn}", type=str)
parser.add_argument("--moS", help="month to plot: 1,..., 12 or to start averaging", type=int, required=True)
parser.add_argument("--moE", help="month to end averaging, default = moS", type=int)
parser.add_argument("--fday", help=f">0: use daily outp day=fday, =0: monthly default={fday}", type=int)
parser.add_argument("--imnth", help=f">0: interp PIOMAS rlx filds to daily for accuracy, default {imnth}", \
                   type=int)
parser.add_argument("--varnm", help="field to plot: ithkn or iconc", type=str, required=True)
parser.add_argument("--nexp0", help=f"reference test nmb, default={nexp0}", type=int)
parser.add_argument("--nexp", help="test run nmb to compare: 2,3,4,5,32,...", type=int, required=True)
args = parser.parse_args()

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
varnm  = args.varnm if args.varnm else None
yr_plt = args.yr if args.yr else None
moS    = args.moS if args.moS else None
moE    = args.moE if args.moE else moS
regn_name = args.regn if args.regn else regn
nexpC  = args.nexp if args.nexp is not None else nexp0
nexpT  = args.nexp if args.nexp is not None else None
use_mnth = not (args.fday and args.fday > 0)
intrpm = args.imnth if args.imnth is not None else 1
intrpm = intrpm if intrpm >= 0 else 1


if yr_plt is None:
  if regn_name == 'NEP':
    yr_plt = 2001
  elif regn_name == 'ARC':
    yr_plt = 1995


interp_mnthly = intrpm >0   # for more accurate comparison, do time interpolation of PIOMAS 
                            # to get mnthly mean values, similar to how it is done in SIS2
                            # when deriving iconc ithkn for day=d0 from PIOMAS target fields

# If interpolating PIOMAS rlx into daily, then do not use fday:
if interp_mnthly:
  fday = False

mstart = 1  # current test runs all started on Jan 1, 2001

# relax hours:
# assumed experiment numbering is x1 - no relaxation, x2 - 1hr run, x3 - 24hr, etc:
RLXH = [0,1,24,120,360]

if regn_name == 'NEP':
  pthtest = f'/archive/Dmitry.Dukhovskoy/fre/NEP/test_ice_relax/NEPphys_expt{nexpC:02d}/{yr_plt}-{mstart:02d}'
  pthrlx  = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax' 
elif regn_name == 'ARC':
  pthtest = f'/archive/Dmitry.Dukhovskoy/fre/ARC12/test_ice_relax/ARCphys_expt{nexpT:02d}/{yr_plt}-{mstart:02d}'
  pthrlx  = '/work/Dmitry.Dukhovskoy/ARC12/irlx'
  
f_rlxfld2 = True

flnm_out = 'ice_month.nc'
prfx = ''

if use_mnth:
  print(f"Difference {varnm} Test {nexpT:02d} using monthly ice_month MM={moS}:{moE}")
else:
  print(f"Difference {varnm} Test {nexpT:02d} using day 15 from ice_daily MM={moS}:{moE}")


def read_sis2_testrun(dnmb0, pthtest, prfx, varnm, use_mnth):
  dv0  = mtime.datevec(dnmb0)
  YR0, MM0, DD0 = dv0[:3]
  dref = dnmb0 - mtime.datenum([1993,1,1])

  print(f'{pthtest} Plot date: {YR0}/{MM0}/{DD0}')

  if use_mnth:
    floutp = 'ice_month.nc'
  else:
    floutp = 'ice_daily.nc'

  if len(prfx) > 0:
    flice_name = f'{prfx}.{floutp}'
  else:
    flice_name  = floutp

  dfsis2 = os.path.join(pthtest, flice_name)

  print(f'Reading {dfsis2}')

  dset   = xarray.open_dataset(dfsis2, decode_times=False)

  TIME = dset['time'].data
  DTM = np.abs(TIME-dref)
  itime = np.argmin(DTM)
  assert DTM[itime] < 15, f'Given date: {YR0}/{MM0}/{DD0} - Check dates in the arch file {dfsis2}'

  HIce = dset['sithick'].isel(time=itime).data
  CIce = dset['siconc'].isel(time=itime).data
  if varnm == 'iconc':
    A2d = CIce
  elif varnm == 'ithkn':
    A2d = CIce*HIce

  return A2d

def read_relax_piomas(dnmb0, diclim, varnm):
  """
    Read target PIOMAS fields, monthly values
  """
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]

  match varnm:
    case('ithkn'):
      ifld = 'ithkn'
    case('iconc'):
      ifld = 'iarea'

  print(f'Reading relax fields from {diclim}')
  #ds_rlx = xarray.open_dataset(diclim)
  with xarray.open_dataset(diclim) as ds_rlx:
    Time = ds_rlx['time'].data
    TM = mmisc.convert_nptime_to_datenum(Time)
    #dnmb0 = mtime.datenum([YR0,MM0,15,12])
    D = abs(TM-dnmb0)
    itime = np.argmin(D)
    dv0 = mtime.datevec(TM[itime])
    assert dv0[0]==YR0, f'Requested YR={YR0}, year in rlx file={dv0[0]}'
    assert dv0[1]==MM0, f'Requested month={MM0}, month in rlx file={dv0[1]}'

    A2dS = ds_rlx[ifld].isel(time=itime).data

  return A2dS


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

if regn_name == 'NEP':
  expt       = 'seasonal_daily'
  pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
  fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
  ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
  hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
  hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
  dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
  dfgrid_mom = os.path.join(pthtopo, fgrid)

  # Hgrid lon. lat:
  hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
else:
  pthtopo = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'

  dflarc  = os.path.join(pthtopo,'ocean_hgrid.nc')
  hlon, hlat = mmom6.read_mom6grid(dflarc, grdpnt='hgrid')

  dtopo = os.path.join(pthtopo,'ocean_topog.nc')
  ds_topo = xarray.open_dataset(dtopo)
  HH = -ds_topo['depth'].data
  jdm, idm = HH.shape

  assert HH[300,200] < 0., f'Check sign of topography, ocean pnts should be < 0'


# Read ice fields from the test experiments:
icc = 0
for mo_plt in range(moS,moE+1):
  dnmb0 = mtime.datenum([yr_plt,mo_plt,15])
  Am = read_sis2_testrun(dnmb0, pthtest, prfx, varnm, use_mnth)

  icc += 1
  if icc == 1:
    A1 = Am.copy()
  else:
    A1 = A1 + Am

if icc > 1:
  A1 = A1/icc


# PIOMAS target Relax fields:
# Read saved relax. fields:
kcc = 0
for mo_plt in range(moS,moE+1):
  dnmb0 = mtime.datenum([yr_plt,mo_plt,15])

  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  YR1 = YR0
  YR2 = YR0+1
  if regn_name == 'NEP':
    flrlx = f'PIOMASv21_ithkn_iconc_{YR1}_{YR2}_monthly.nc'
  else:
    flrlx = f'PIOMASv21_ARC12_ithkn_iconc_{YR1}_{YR2}_monthly.nc'
  diclim = os.path.join(pthrlx, flrlx)

  if interp_mnthly:
    Am = msisrlx.mnthly_PIOMAS_linear_daily(diclim,dnmb0,varnm)
  else:
    Am = read_relax_piomas(dnmb0, diclim, varnm)

  kcc += 1
  if kcc == 1:
    A2 = Am.copy()
  else:
    A2 = A2 + Am

if kcc > 1:
  A2 = A2/kcc

dI = A1-A2

# -------------------
#
# Plot ice fields
#
# -------------------
print(" Plotting ...")
clrmp_dlt = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBu_r')
clrmp_dlt.set_bad(color=[0.2,0.2,0.2])
if varnm == 'iconc':
  clrmp = mclrmps.colormap_conc()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 1.
  dmin = -1
  dmax = 1
elif varnm == 'ithkn':
  clrmp = mclrmps.colormap_ice_thkn()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 4.
  dmin = -1.5
  dmax = 1.5

rindx = int(nexpT % 10)
rlx_time = RLXH[rindx-1]

if moE == moS:
  sttl = f"Rlx={rlx_time} hrs, Diff {varnm} test{nexpT:02d}-PIOMAS_rlx avrg: {yr_plt}/{mo_plt:02d}"
else:
  sttl = f"Rlx={rlx_time} hrs, Diff {varnm} test{nexpT:02d}-PIOMAS_rlx avrg: {yr_plt} {moS:02d}-{moE:02d}"


# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
if regn_name == 'NEP':
  m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
              projection='stere', lat_ts=60, lat_0=65, lon_0=-175)
  xR, yR = m(hlon, hlat)
else:
  lon0 = 180.
  lat0 = 70.
  res  = 'l'
  m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
  m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
  xR, yR = m(hlon,hlat)

if use_mnth:
  sinfo = f'test run: {pthtest}/ice_month.nc\n'
else:
  sinfo = f'test run: {pthtest}/ice_daily.nc\n'
sinfo = sinfo + f'rlx fields: {diclim}'

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))


f_pnt = False
if f_pnt:
  xL1 = -164.223
  yL1 = 71.231
  xp, yp = m(xL1,yL1)

  xL2 = -174.69
  yL2 = 62.19
  xp2, yp2 = m(xL2,yL2)

img = m.pcolormesh(xR, yR, dI, cmap=clrmp_dlt, vmin=dmin, vmax=dmax)

f_cntr = True
if f_cntr:
  cntr_clr = [0.95, 0.95, 0.95]
  hcntrs = [2,4,6,8,10,12,14,16,18,20]
 
  CS = ax1.contour(xR, yR, dI, hcntrs, linestyles='solid', colors=[cntr_clr], linewidths=1)
  ax1.clabel(CS, inline=1, fontsize=10)


if regn_name == 'ARC':
  ax1.set_xlim([42e5, 95e5])  
  ax1.set_ylim([43e5, 108e5])

ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(dmin,dmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.025, 0.8, 0.05])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')

btx = 'plot_diffseaice_irlx.py'
bottom_text(btx, pos=[0.2, 0.01])


