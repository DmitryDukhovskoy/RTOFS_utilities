"""
  Compute monthly ice climatologies for ice conc
  from NEP seas f/casts
 
  Note that year 2010 and month 1 designate the initialization time
  NEP  monthly fields have 12 f/cast months in each file

  Only 12 months of the f/cast are saved in the NEP files
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
parser.add_argument("--MMI", help="init month of NEP f/cast: 1, ..., 12", type=int)
parser.add_argument("--ensmb", help="ensemble number, 1,..., 15", type=int)
parser.add_argument("--exptnmb", help="NEP seas f/cast expt=2,3 (no /with ice relax)", type=int)
parser.add_argument("--varnm", help="field: ithkn or iarea/iconc", type=str)
args = parser.parse_args()

f_save = True
# Years in the relax file also used in the rlx file name:
navrg = 5
YRS = 2010  # init yr
YRE = YRS+navrg-1
MMI = 1     # init month
ifld = 'iarea'  # ithkn, iarea
ens_nmb = 1  # NEP ensemble #
interp_NEP = True # interpolate to NEP and save on both grids (NEP & NEP), otherwise - only NEP
expt     = "seasonal_daily"
expt_nmb = 2 

if args.YRS:
  YRS = args.YRS
  YRE = YRS+navrg-1
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI
if args.ensmb:
  ens_nmb = args.ensmb
if args.exptnmb:
 expt_nmb=args.exptnmb
if args.varnm:
  ifld = args.varnm
  if ifld=='iconc' or ifld=='iarea':
    ifld = 'siconc'
  elif ifld=='ithkn' or ifld=='ithk':
    ifld = 'sithick'

varnm = ifld
runname  = f"NEPphys_frcst_dailyOB-expt{expt_nmb:02d}"

if not f_save:
  print(f'WARNING: fields wont be saved, f_save flag is off !!!\n')

# Saved climatologies:
ICLIM=[[1993,1997],[1995,1999],[2000,2004],[2005,2009],[2010, 2014],[2015,2019],[2016,2020]]
ICLIM=np.array(ICLIM)

# There are no MMI=01 in 1993, so push 1 year if YRI=1993:
if YRS==1993 and MMI==1:
  YRS=YRS+1
  YRE=YRE+1
  if YRE>2020:
    YRE=2020

if YRS==2016 and MMI>1:
  print(f'For f/casts initialized MMI>1 no 2020 runs, use 2015-2019 clim')
  assert YRS<2016, 'Change YRS<2016 for 5-yr climatology'

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)
pthoutp1 = pthseas['MOM6_NEP']['seasonal_daily']['pthoutp'].format(expt_nmb=expt_nmb)

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

# Read NEP fields:
# Note that NEP thickness is not "volume/unit area" (like in CICE output or PIOMAS)
# and needs to be multiplied by partial area
icc = 0
for YRI in range(YRS,YRE+1):
  # Process and save climatology by forecasts lead time (not calendar months!)
  print(f'Processing {YRI}')

  if icc == 0:
    ASUM = np.zeros((12,jdim,idim))

  for MMF in range(1,13):
    # Forecast months
    itime = MMF-1

    pthfcst1 = os.path.join(pthoutp1,f'{YRI}-{MMI:02d}-e01','history')
    dcice1 = os.path.join(pthfcst1,f'ice_month.nc')
    #MMF = msisrlx.cal_mo_to_fcast(MMI,MM) # forecast month #
    imo = MMF-1
    mcal = np.arange(MMI,MMI+12)
    mcal = np.where(mcal>12, mcal-12, mcal)
    mcal = mcal.astype(int)

    ds = xarray.open_dataset(dcice1)
    H2d = ds['sithick'].isel(time=itime).data
    C2d = ds['siconc'].isel(time=itime).data
    if ifld == 'siconc':
      A2d = C2d.copy()
    elif ifld == 'sithick':
      A2d = H2d*C2d         # ice m ---> m3/m2 

    hmax = np.nanmax(A2d)
    hmin = np.nanmin(A2d)
    print(f'N={itime+1} M={mcal[itime]} {ifld} min/max: {hmin:.2f}/{hmax:.2f}')

    ASUM[itime,:,:] = ASUM[itime,:,:] + A2d

  icc += 1

ASUM = ASUM / icc

kdim,jdim,idim = ASUM.shape
darr_ice = xarray.DataArray(ASUM, dims=("months","jdim","idim"), \
                  coords={"months": np.arange(kdim), \
                          "jdim": np.arange(jdim), \
                          "idim": np.arange(idim)})
darr_months = xarray.DataArray(mcal, dims=("months"), \
                     coords={"months": np.arange(kdim)})
dset = xarray.Dataset({"calend_months": darr_months, f"{varnm}": darr_ice})


if f_save:
  pthdump = pthseas['MOM6_NEP']['seasonal_daily']['pthsis2'].format(expt_nmb=expt_nmb)
  floutp = f'NEPseasfcast_{varnm}_clim_{YRS}_{YRE}_MI{MMI:02d}e{ens_nmb:02d}.nc'
  dflout = os.path.join(pthdump,floutp)
  dset['calend_months'].attrs['long_name']='Calendar months during the forecast'
  if ifld == 'siconc':
    dset[f'{varnm}'].attrs['long_name']='ice partial area'
  elif ifld == 'sithick':
    dset[f'{varnm}'].attrs['long_name']='mean cell thickness or ice volume per unit area, m3/m2'

  # Add global attributes:
  dset.attrs.update({
    "info": "Interpolated NEP {varnm} climatology to NEP SIS2 grid",
    "code": "calc_seasfcst_ice_clim.py"
  })

  print(f'Dumping {varnm} climtology --> {dflout}')
  dset.to_netcdf(dflout, format='NETCDF3_64BIT', engine='netcdf4')
 

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

  btx = 'calc_seasfcst_ice_clim.py'
  bottom_text(btx, pos=[0.2, 0.01])


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

  xR, yR = m(hlon,hlat)

  fgnmb=1
  sttlS = f'NEP clim init: {YRI}/{MMI:02d}-e{ens_nmb:02d}, {ifld} {YRS}-{YRE} {MM0}'
  plot_ice(fgnmb, m, xR, yR, A2d, clrmp, rmin, rmax, sttlS)

