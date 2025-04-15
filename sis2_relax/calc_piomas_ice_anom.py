"""
  Compute monthly PIOMAS ice anomalies for ice conc and thickn. 
  climatology computed in calc_piomas_ice_clim.py

  Save anomalies by years
 
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
args = parser.parse_args()

plot_fields = True
plot_piomas = True
# Years in the relax file also used in the rlx file name:
YRS = 2010  # init yr
YRE = 2014
MMI = 1     # init month
f_save = True

if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
#if args.varnm:
#  ifld = args.varnm
#  if ifld=='iconc' or ifld=='iarea':
#    ifld = 'siconc'
#  elif ifld=='ithkn' or ifld=='ithk':
#    ifld = 'sithick'

# 5-yr climatologies:
ICLIM=[[1990,1994],[1995,1999],[2000,2004],[2005,2009],[2010, 2014],[2015,2019],[2020,2023]]
ICLIM=np.array(ICLIM)

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
LMsk = np.where(HH<0, 1, 0)
jdm, idm = HH.shape


# Read PIOMAS relax. fields:
icc = 0
for YRI in range(YRS,YRE+1):
  # Get monthly climatology:
  iclm = np.where((ICLIM[:,0] <= YRI) & (ICLIM[:,1] >= YRI))[0]
  assert(len(iclm)>0), f'Could not find clim. time window for {YRI}'
  iclm = iclm[0]
  YRC1 = ICLIM[iclm,0]
  YRC2 = ICLIM[iclm,1]

  pthpiomas = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax'
  floutp = f'piomasV21_iconc_ithkn_clim_{YRC1}_{YRC2}.nc'
  dflout = os.path.join(pthpiomas,floutp)
  print(f'Reading climtology --> {dflout}')
  dsclim = xarray.open_dataset(dflout)
   
  Hanom = np.zeros((12,jdm,idm))
  Canom = np.zeros((12,jdm,idm))
  for MM in range(1,13):
    imo = MM-1
    print(f'Processing {YRI}/{MM}')
    # Read PIOMAS interpolated to NEP MOM6-SIS2 grid:
    fice   = f'PIOMAS_ithkn_iconc_{YRI}_{YRI+1}_monthly.nc'
    dfice  = os.path.join(pthpiomas,fice)
    ds_rlx = xarray.open_dataset(dfice)

    Time = ds_rlx['time'].data
    TM   = mmisc.convert_nptime_to_datenum(Time)
    dnmb0 = mtime.datenum([YRI,MM,15,12])
    D = abs(TM-dnmb0)
    itime = np.argmin(D)
    assert(D[itime] == 0), f"Could not find requested time {YRI}/{MM}"
    dv0 = mtime.datevec(TM[itime])
    print(f'PIOMAS field: {dv0[0]}/{dv0[1]}/{dv0[2]}')

    H2d = ds_rlx['ithkn'].isel(time=itime).data
    C2d = ds_rlx['iarea'].isel(time=itime).data
    Hclim = dsclim['ice_thikness'].isel(months=imo).data
    Cclim = dsclim['ice_conc'].isel(months=imo).data

    # Add land mask:
    C2d = np.where(HH>=0, np.nan, C2d)
    H2d = np.where(HH>=0, np.nan, H2d)

    Hanom[imo,:,:] = H2d - Hclim
    Canom[imo,:,:] = C2d - Cclim
   
  mcal = np.arange(1,13) 
  darr_hice = xarray.DataArray(Hanom, dims=("nmonths","jdim","idim"),\
                     coords={"nmonths": np.arange(12), \
                             "jdim": np.arange(jdm), \
                             "idim": np.arange(idm)})
  darr_cice = xarray.DataArray(Canom, dims=("nmonths","jdim","idim"),\
                     coords={"nmonths": np.arange(12), \
                             "jdim": np.arange(jdm), \
                             "idim": np.arange(idm)})
  darr_months = xarray.DataArray(mcal, dims=("nmonths"), \
                       coords={"nmonths": np.arange(12)})
  dset = xarray.Dataset({"Calendar_months": darr_months, \
                         "ice_thkn_anom": darr_hice, \
                         "ice_conc_anom": darr_cice})

  if f_save:
    pthpiomas = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax'
    floutp = f'piomasV21_iconc_ithkn_anom_{YRI}.nc'
    dflout = os.path.join(pthpiomas,floutp)
    print(f'Dumping anomalies --> {dflout}')
    dset.to_netcdf(dflout, format='NETCDF3_64BIT', engine='netcdf4')


# Check
f_check = False
if f_check:
  plt.ion()

  MM0 = 9  # month to check
  it0 = MM0-1
  ifld = 'siconc'

  match ifld:
    case('sithick'):
      clrmp = mclrmps.colormap_ice_thkn()
      rmin = 0.
      rmax = 4.
      A2d = Canom[it0,:,:].squeeze()
    case('siconc'):
      clrmp = mclrmps.colormap_conc()
      rmin = 0.
      rmax = 1.
      A2d = Hanom[it0,:,:].squeeze()

  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  A2d = np.where(A2d>=999, np.nan, A2d)

  # Stereographic Map projection:
  from mpl_toolkits.basemap import Basemap, cm
  #m = Basemap(width=5000*1.e3,height=5000*1.e3, resolution='l',\
  #            projection='stere', lat_ts=50, lat_0=62, lon_0=-165)
  m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

  xR, yR = m(LON, LAT)

  fgnmb=1
  sttlS = f'PIOMAS clim {ifld} {YRS}-{YRE} {MM0}'
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  #m.drawcoastlines()
  #m.drawparallels(np.arange(-90.,120.,10.))
  #m.drawmeridians(np.arange(-180.,180.,10.))

  #img = ax1.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  img = ax1.pcolormesh(A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  #  img = ax1.pcolormesh(RLXHR, cmap=clrmp)

  sttl = f'{ifld} anomaly {YRI}/{MM0}'
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

  btx = 'calc_PIOMAS_ice_clim.py'
  bottom_text(btx, pos=[0.2, 0.01])


