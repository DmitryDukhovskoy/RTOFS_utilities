"""
  Plot SPEAR and PIOMAS monthly  ice anomalies
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
parser.add_argument("--YR", help="year to plot anomalies : 1993, ..., 2022", type=int)
parser.add_argument("--MM", help="calendar month to plot: 1, ..., 12", type=int)
parser.add_argument("--MMI", help="init month of SPEAR f/cast: 1, ..., 12", type=int)
parser.add_argument("--ensmb", help="SPEAR ensemble number, 1,..., 15", type=int)
parser.add_argument("--varnm", help="field: ithkn or iarea/iconc", type=str)
args = parser.parse_args()

YR  = 2010
MM  = 9
MMI = 1
ens_nmb = 1

if args.YR:
  YR = args.YR
if args.MM:
  MM = args.MM
if args.MMI:
  MMI = args.MMI
if args.ensmb:
  ens_nmb = args.ensmb
if args.varnm:
  ifld = args.varnm
  if ifld=='iconc' or ifld=='iarea':
    ifld = 'siconc'
    ifld_piomas = 'ice_conc_anom'
    varnm = 'iconc'
  elif ifld=='ithkn' or ifld=='ithk':
    ifld = 'sithick'
    ifld_piomas = 'ice_thkn_anom'
    varnm = 'ithkn'

# Saved climatologies:
ICLIM=[[1990,1994],[1995,1999],[2000,2004],[2005,2009],[2010, 2014],[2015,2019],[2020,2023]]
ICLIM=np.array(ICLIM)

# NEP grid:
fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

expt_name = "seasonal_daily"
pthtopo    = pthseas['MOM6_NEP'][expt_name]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt_name]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt_name]["ftopo"]
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
jdim, idim = HH.shape

iclm = np.where((ICLIM[:,0] <= YR) & (ICLIM[:,1] >= YR))[0]
assert(len(iclm)>0), f'Could not find clim. time window for {YR}'
iclm = iclm[0]
YRC1 = ICLIM[iclm,0]
YRC2 = ICLIM[iclm,1]

# Read SPEAR anomalies:
pthspear = '/work/Dmitry.Dukhovskoy/anls_output/spear_ice'
flspear = f'spear_{ifld}_monthly_anom_{YR}{MMI:02d}.nc'
dflspear = os.path.join(pthspear, flspear)

dspear = xarray.open_dataset(dflspear)
mcal = dspear['Calendar_months'].data
imo  = np.where(mcal == MM)[0][0]

Asp = dspear[f"{ifld}_anom"].isel(nmonths=imo).data

# Read PIOMAS anomalies:
pthpiomas = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax'
flpiomas = f'piomasV21_iconc_ithkn_anom_{YR}.nc'
dflpio = os.path.join(pthpiomas,flpiomas)

dspio = xarray.open_dataset(dflpio)
iMM = MM-1     # PIOMAS fields are by cal. months, 1, ..., 12
Apio = dspio[ifld_piomas].isel(nmonths=iMM).data


# Difference plot:
dA = Asp - Apio


CLRS = [[0.6, 0.02, 0.6],
        [0.2, 0.38, 1],
        [0., 0.8, 0.5],
        [0.2,1.,0.8],
        [1, 1, 1],
        [1, 0.9, 0.85],
        [1, 0.4, 0.4],
        [0.9, 0.6,0],
        [0.6, 0.2, 0]]

clrmp = mclrmps.colormap_posneg_uneven(CLRS)
clrmp.set_bad(color=[0.6,0.6,0.6])
if varnm == 'iconc':
  rmin = -1.
  rmax = 1.
  dmin = -0.5
  dmax = 0.5
else:
  rmin = -2.
  rmax = 2.
  dmin = -1.
  dmax = 1.

clrmp_dlt = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='YlGnBu_r')
clrmp_dlt.set_bad(color=[0.6,0.6,0.6])


def plot_ice(ax1, xR, yR, A2d, clrmp, rmin, rmax, sttl, btx=[]):
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

  if len(btx)>0:
    bottom_text(btx, pos=[0.2, 0.01])

  return ax1, ax2

plt.ion()

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
          projection='stere', lat_ts=60, lat_0=65, lon_0=-175)
xR, yR = m(hlon, hlat)

fgnmb=1
sttlS = f'SPEAR {varnm} anom {YR}/{MM} wrt {YRC1}-{YRC2}, init {MMI:02d}-e{ens_nmb:02d}'
sttlP = f'PIOMASv2.1 {varnm} anom {YR}/{MM} wrt {YRC1}-{YRC2}'
sttlD = 'SPEAR - PIOMAS difference'
btx = 'plot_SPEAR_PIOMAS_anom.py'

fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.05, 0.55, 0.4, 0.4])
ax1, ax12 = plot_ice(ax1, xR, yR, Asp, clrmp, rmin, rmax, sttlS)
ax2 = plt.axes([0.54, 0.55, 0.4, 0.4])
ax2, ax22 = plot_ice(ax2, xR, yR, Apio, clrmp, rmin, rmax, sttlP)
ax3 = plt.axes([0.3, 0.08, 0.4, 0.4])
ax3, ax32 = plot_ice(ax3, xR, yR, dA, clrmp_dlt, dmin, dmax, sttlD, btx=btx)







