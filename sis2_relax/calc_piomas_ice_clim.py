"""
  Compute monthly PIOMAS ice climatologies for ice conc and thickn. 
 
  usage: 
  plot_piomas_ice_month_stere.py --YRS 2010 --YRE 2020

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
#ifld = 'iarea'  # ithkn, iarea
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
#ICLIM=[[1990,1994],[1995,1999],[2000,2004],[2005,2009],[2010, 2014],[2015,2019],[2020,2023]]
#ICLIM=np.array(ICLIM)
pthdata = '/work/Dmitry.Dukhovskoy/data/PIOMAS_ice'
varthck = 'heff'
varconc = 'area'

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
  print(f'Processing {YRI}')

  flthck  = f'piomas_heff{YRI}_v21.nc'
  flconc  = f'piomas_area{YRI}_v21.nc'
  dflthkn = os.path.join(pthdata, flthck)
  dflconc = os.path.join(pthdata, flconc)

  ds_thkn = xarray.open_dataset(dflthkn)
  ds_conc = xarray.open_dataset(dflconc)
  LAT  = ds_thkn['lat_scaler'].data
  LON  = ds_thkn['lon_scaler'].data

  if icc == 0:
   ydim, xdim = LAT.shape
   HSUM = np.zeros((12,ydim,xdim))
   CSUM = np.zeros((12,ydim,xdim))

  for imo in range(12):
    H2d   = ds_thkn[varthck].data[imo,:].squeeze()  # thikness, m
    C2d   = ds_conc[varconc].data[imo,:].squeeze()  #conc

    # Get rid off the land values along the southern boundary:
    nbnd = 4
    for ik in range(nbnd):
      H2d[ik,:] = H2d[nbnd,:]
      C2d[ik,:] = C2d[nbnd,:]

    # Get rid of nans - fill land:
    H2d = np.where(H2d>9999., np.nan, H2d)
    C2d = np.where(C2d>9999., np.nan, C2d)
    H2df = mmom6.fill_land3d(H2d, sinfo=f'thkn mo={imo+1}')
    C2df = mmom6.fill_land3d(C2d, verb=0)
    C2df  = np.where(C2df > 1., 1., C2df)
    C2df  = np.where(C2df < 0.0, 0.0, C2df)

    HSUM[imo,:,:] = HSUM[imo,:,:] + H2df
    CSUM[imo,:,:] = CSUM[imo,:,:] + C2df

  icc += 1

HSUM = HSUM / icc
CSUM = CSUM / icc

# Interpolate to NEP MOM6-SIS2 grid:
# Find gmapi - indices for bi-polar interpolation
import mod_regmom as mrmom
pthsis = gridfls['MOM6_NEP'][run_name]['pthsis']
fgmapi  = f'PIOMAS_mom6_NEP_gmapi_{jdm}x{idm}.pkl'
dfgmapi = os.path.join(pthsis, fgmapi)
print(f'Loading gmapi <-- {dfgmapi}')
with open(dfgmapi, 'rb') as fid:
  IMOM, JMOM, INDX, JNDX = pickle.load(fid)

# interpolate onto MOM6 grid
Hclim = np.zeros((12,jdm,idm))
Cclim = np.zeros((12,jdm,idm))
for imo in range(12):
  print(f'Interpolating to NEP MOM6-SIS2 grid, month={imo+1}')
  h2d  = HSUM[imo,:,:].squeeze()
  h2di = msisrlx.interp2Dfld(h2d, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
  c2d  = HSUM[imo,:,:].squeeze()
  c2di = msisrlx.interp2Dfld(c2d, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
  c2di = np.where(c2di>0.99, 1.0, c2di) # interp error 1--> 0.99, also caps max conc = 1

  Hclim[imo,:,:] = h2di.copy()
  Cclim[imo,:,:] = c2di.copy()


mcal = np.arange(1,13)

if f_save:
  kdim,jdim,idim = Hclim.shape
  darr_hice = xarray.DataArray(Hclim, dims=("months","jdim","idim"),\
                     coords={"months": np.arange(kdim), \
                             "jdim": np.arange(jdim), \
                             "idim": np.arange(idim)})
  darr_cice = xarray.DataArray(Cclim, dims=("months","jdim","idim"),\
                     coords={"months": np.arange(kdim), \
                             "jdim": np.arange(jdim), \
                             "idim": np.arange(idim)})
  darr_months = xarray.DataArray(mcal, dims=("months"), \
                       coords={"months": np.arange(kdim)})
  dset = xarray.Dataset({"Calendar_months": darr_months, "ice_conc": darr_cice, \
                         "ice_thikness": darr_hice})


  pthpiomas = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax'
  floutp = f'piomasV21_iconc_ithkn_clim_{YRS}_{YRE}.nc'
  dflout = os.path.join(pthpiomas,floutp)
  print(f'Dumping climtology --> {dflout}')
  dset.to_netcdf(dflout, format='NETCDF3_64BIT', engine='netcdf4')

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

  btx = 'calc_PIOMAS_ice_clim.py'
  bottom_text(btx, pos=[0.2, 0.01])


# Check
f_check = False
if f_check:
  plt.ion()

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

  MM0 = 9  # month to check
  it0 = MM0-1
  A2d = CSUM[it0,:,:].squeeze()
  #A2d = HSUM[it0,:,:].squeeze()
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
  plot_ice(fgnmb, xR, yR, A2d, clrmp, rmin, rmax, sttlS)

