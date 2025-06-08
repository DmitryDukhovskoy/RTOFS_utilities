"""
  Plot sea ice conc/thickness 
  from test simulations

  Plot difference of daily ice fields (icem*.nc) 
  from test experiments with relaxation and / or target fields (monthly mean)

 
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
importlib.reload(mutob)

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="year to plot: 1993, ..., 2020", type=int)
parser.add_argument("--jday", help="year day to plot: 1, ..., 366", type=int)
parser.add_argument("--mo", help="month to plot: 1,..., 12", type=int)
parser.add_argument("--day", help="day to plot: 1, ..., 31", type=int)
parser.add_argument("--varnm", help="field to plot: ithkn or iconc", type=str)
parser.add_argument("--run1", help="1st test run to compare: irlx1, irlx2, ... or rlxfld", type=str)
parser.add_argument("--run2", help="2nd test run to compare: irlx1, irlx2, ... or rlxfld", type=str)
args = parser.parse_args()

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
varnm  = 'ithkn'  # iconc or ithkn
expt1 = 'irlx3'
expt2 = 'irlx2'

jday_plt = 0

if args.varnm:
  varnm = args.varnm
if args.yr:
  yr_plt = args.yr
if args.jday:
  jday_plt = args.jday
if args.run1:
  expt1 = args.run1
if args.run2:
  expt2 = args.run2
if args.mo:
  mo_plt = args.mo
if args.day:
  day_plt = args.day

f_rlxfld1 = False
if expt1[:4] == 'irlx':
  pthtest1 = f'/work/Dmitry.Dukhovskoy/tmp/test_{expt1}'
elif expt1[:4] == 'rlxf':
  pthtest1 = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax' 
  f_rlxfld1 = True

f_rlxfld2 = False
if expt2[:4] == 'irlx':
  pthtest2 = f'/work/Dmitry.Dukhovskoy/tmp/test_{expt2}'
elif expt2[:4] == 'rlxf':
  pthtest2 = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax' 
  f_rlxfld2 = True

outfld  = 'icem'
prfx = ''  # 19930401 - time stamp used in SIS2 output in file names, note that find
           # closest archive file does not work for 19930401.icem*.nc file names
           # rename files using ./rename_archive_v0.sh 0 in the output dir
MMI = 4
ndav1 = 5
if expt1 == 'irlx3':
  ndav1 = 1
ndav2 = 5
if expt2 == 'irlx3':
  ndav2 = 1

i0=j0=0
# Test point in Fortran indices:
# make it <0 not to show
#iF0 = 224
#jF0 = 742
iF0 = 184
jF0 = 655
i0 = iF0-1 ; j0 = jF0-1

if jday_plt >0 and jday_plt <=366:
  dnmbR  = mtime.jday2dnmb(yr_plt,jday_plt)
else:
  dnmbR  = mtime.datenum([yr_plt,mo_plt,day_plt])  # day to plot
dvR     = mtime.datevec(dnmbR)

def read_test_run(dnmb0, MMI, outfld, pthtest, prfx, varnm):
  # Find init year of the f/cast:
  YRI = manseas.yr_init_fcst_from_datenum(dnmb0, MMI)
  # Find closest output:
  YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthtest, dnmb0, fld=outfld)
  dv0  = mtime.datevec(dnmb0)
  YR0, MM0, DD0 = dv0[:3]
  jday0   = int(mtime.date2jday([YR0,MM0,DD0]))

  print(f'{pthtest} Plot date: {dvR[0]}/{dvR[1]}/{dvR[2]}')

  if len(prfx) > 0:
    flice_name = f'{prfx}.icem_{YR0}_{jday0:03d}.nc'
  else:
    flice_name  = f'icem_{YR0}_{jday0:03d}.nc'
  dfsis2 = os.path.join(pthtest, flice_name)

  print(f'Reading {dfsis2}')

  dset   = xarray.open_dataset(dfsis2)

  HIce = dset['sithick'].isel(time=0).data
  CIce = dset['siconc'].isel(time=0).data
  if varnm == 'iconc':
    A2d = CIce
  elif varnm == 'ithkn':
    A2d = CIce*HIce

  return A2d

def read_relax_piomas(dnmb0, pthsis, varnm):
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  flthck  = f'piomas_heff{YR0}_v21.nc'
  varthck = 'heff'
  flconc  = f'piomas_area{YR0}_v21.nc'
  varconc = 'area'

  # Read saved relax. fields:
  YR1 = YR0
  YR2 = YR0+1
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

  match varnm:
    case('ithkn'):
      ifld = 'ithkn'
    case('iconc'):
      ifld = 'iarea'

  A2dS = ds_rlx[ifld].isel(time=itime).data

  return A2dS

# Averaging period:
dnmb_av1 = dnmbR - np.floor(ndav1/2)
#if dnmb_av1 < dnmbI: dnmb_av1=dnmbI
dnmb_av2 = dnmb_av1 + ndav1-1

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

expt       = 'seasonal_daily'
pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
ndav       = pthseas['MOM6_NEP'][expt]['ndav']  # # of days output averaged

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

# Read ice fields from the test experiments:
if f_rlxfld1:
  A1 = read_relax_piomas(dnmbR, pthtest1, varnm)
else:
  A1 = read_test_run(dnmbR,MMI,outfld,pthtest1,prfx,varnm)
if f_rlxfld2: 
  A2 = read_relax_piomas(dnmbR, pthtest2, varnm)
else:
  A2 = read_test_run(dnmbR,MMI,outfld,pthtest2,prfx,varnm)
dI = A1-A2

# -------------------
#
# Plot ice fields
#
# -------------------
clrmp_dlt = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBu_r')
clrmp_dlt.set_bad(color=[0.2,0.2,0.2])
if varnm == 'iconc':
  clrmp = mclrmps.colormap_conc()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 1.
  dmin = -0.5
  dmax = 0.5
elif varnm == 'ithkn':
  clrmp = mclrmps.colormap_ice_thkn()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 4.
  dmin = -0.5
  dmax = 0.5

dv_av1 = mtime.datevec(dnmb_av1)
yrs, mms, dds = dv_av1[:3]
dv_av2 = mtime.datevec(dnmb_av1)
yre, mme, dde = dv_av2[:3]

sttl = f"Diff {varnm} {expt1}-{expt2} avrg: {yrs}/{mms}/{dds}-{yre}/{mme}/{dde}"
if j0 >= 0 and i0 >= 0:
  sttl = sttl + f"\n Test pnt iF0/jF0 = {iF0}/{jF0} {varnm}={dI[j0,i0]:.6f}"

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

sinfo = f'run1: {pthtest1}\n'
sinfo = sinfo + f'run2: {pthtest2}'

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

img = m.pcolormesh(xR, yR, dI, cmap=clrmp_dlt, vmin=dmin, vmax=dmax)

ax1.set_title(sttl)

# Show test pnt:
if j0 >= 0 and i0 >=0:
  xTst, yTst = m(hlon[j0,i0],hlat[j0,i0])
  ax1.plot(xTst,yTst,'o')

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

btx = 'plot_diffseaice_test.py'
bottom_text(btx, pos=[0.2, 0.01])


