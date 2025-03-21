"""
  Check relax fields created from PIOMAS monthly ice thickness and concentration
  stereographic projection

  usage: check_relax_sis2_stere.py --varnm={ithkn,iarea or iconc} --yr=1993 --mo=8

  see: piomas_relaxation_yearly.py

  monthly fields
  1901 - 2010
  https://psc.apl.uw.edu/research/projects/piomas-20c/

"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
import pickle
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
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc
from mod_utils_fig import bottom_text
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="year to plot: 1993, ..., 2020", type=int)
parser.add_argument("--mo", help="month to plot: 1,..., 12", type=int)
parser.add_argument("--varnm", help="field to plot: ithkn or iarea", type=str)
args = parser.parse_args()

plot_fields = True
plot_piomas = False
# Years in the relax file also used in the rlx file name:
YR1 = 1993
YR2 = 1994 
YR0 = 1993   # year to plot

MM0 = 6      # month to plot
ifld = 'iarea'  # ithkn, iarea
file_type = 'monthly'  # monthly, daily, ... or clim
                       # for climatologies, do not need padded time - data will be recycled
                       # for monthly, daily, etc. need -dt and +dt at the beginn/end 

if args.yr:
  YR0 = args.yr
  YR1 = YR0
  YR2 = YR1+1
if args.mo:
  MM0 = args.mo
if args.varnm:
  ifld = args.varnm
  if ifld=='iconc':
    ifld = 'iarea'

if YR0 < YR1 or YR0 > YR2:
  raise Exception(f"year to plot {YR0} is outside the time window in the file: {YR1}/{YR2}")

# Test point in Fortran indices:
# make it <0 not to show
iF0 = 230
jF0 = 700
i0 = iF0-1 ; j0 = jF0-1


fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

# MOM6 NEP topo/grid:
run_name   = 'seasonal_fcst_daily'
pthtopo    = gridfls['MOM6_NEP'][run_name]['pthgrid']
fgrid      = gridfls['MOM6_NEP'][run_name]['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"][run_name]["ftopo"]
outdir     = gridfls['MOM6_NEP'][run_name]['pthoutp']
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
# Hgrid lon. lat:
hlon, hlat  = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

pthsis  = gridfls['MOM6_NEP'][run_name]['pthsis']
pthdata = '/work/Dmitry.Dukhovskoy/data/PIOMAS_ice'
flthck = 'piomas20c.heff.1901.2010.v1.0.nc'
varthck = 'sit'
flconc  = 'piomas20c.area.1901.2010.v1.0.nc'
varconc = 'sic'

dflthkn = os.path.join(pthdata, flthck)
dflconc = os.path.join(pthdata, flconc)

ds_thkn = xarray.open_dataset(dflthkn)
LAT  = ds_thkn['Latitude'].data
LON  = ds_thkn['Longitude'].data

# Read saved relax. fields:
flout = f'PIOMAS_ithkn_iconc_{YR1}_{YR2}_{file_type}.nc'
diclim = os.path.join(pthsis, flout)
ds_rlx = xarray.open_dataset(diclim)
Time = ds_rlx['time'].data
TM = mmisc.convert_nptime_to_datenum(Time)
dnmb0 = mtime.datenum([YR0,MM0,15,12])
D = abs(TM-dnmb0)
itime = np.argmin(D)
dv0 = mtime.datevec(TM[itime])
assert dv0[0]==YR0, f'Requested YR={YR0}, year in rlx file={dv0[0]}'
assert dv0[1]==MM0, f'Requested month={MM0}, month in rlx file={dv0[1]}'

A2dS = ds_rlx[ifld].isel(time=itime).data

# Read PIOMAS field:
match ifld:
  case('ithkn'):
    varnm = varthck
    dfpiomas = os.path.join(pthdata,flthck)
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 4.
  case('iarea'):
    varnm = varconc
    dfpiomas = os.path.join(pthdata,flconc)
    clrmp = mclrmps.colormap_conc()
    rmin = 0.
    rmax = 1.

clrmp.set_bad(color=[0.2, 0.2, 0.2])
A2dP = msisrlx.read_PIOMAS(YR0, MM0, dfpiomas, varnm)


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

  btx = 'check_relax_sis2_stere.py'
  bottom_text(btx, pos=[0.2, 0.01])



plt.ion()

if plot_fields: 
  A2dS = np.where(HH>=0, np.nan, A2dS)

  # Stereographic Map projection:
  from mpl_toolkits.basemap import Basemap, cm
  #m = Basemap(width=5000*1.e3,height=5000*1.e3, resolution='l',\
  #            projection='stere', lat_ts=50, lat_0=62, lon_0=-165)
  m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

  xR, yR = m(hlon, hlat)
  xRp, yRp = m(LON, LAT)

  sttlS = f'Relaxation {ifld} SIS2 from PIOMAS {YR0}/{MM0}'
  if j0 >= 0 and i0 >= 0:
    sttlS = sttlS + f"\n Test pnt iF0/jF0 = {iF0}/{jF0} {ifld}={A2dS[j0,i0]:.6f}"

  sttlP = f'{ifld} PIOMAS {YR0}/{MM0}'

  fgnmb=1
# Show test pnt:
  if j0 >= 0 and i0 >=0:
    xTst, yTst = m(hlon[j0,i0],hlat[j0,i0])
    plot_ice(fgnmb, xR, yR, A2dS, clrmp, rmin, rmax, sttlS, xTst=xTst, yTst=yTst)
  else:
    plot_ice(fgnmb, xR, yR, A2dS, clrmp, rmin, rmax, sttlS)

  if plot_piomas:
    fgnmb=2
    plot_ice(fgnmb, xRp, yRp, A2dP, clrmp, rmin, rmax, sttlP)





