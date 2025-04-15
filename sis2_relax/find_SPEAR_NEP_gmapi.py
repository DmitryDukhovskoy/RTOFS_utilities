"""
  Find gmapi indices: 4 vertices of SPEAR grid for each NEP grip
  for bilinear interpolation of SPEAR --> NEP 
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

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)
spear_dir = os.path.join(config['filesystem']['spear_month_ens'], 'monthly_clim')

# Check if mapping indices exist, gmapi:
dirgmapi = config['filesystem']['spear_mom_gmapi']
flgmaph  = 'spear2mom_NEP_full_gmapi_hpnt.nc'
#flgmaph  = 'spear2mom_NEP_full_gmapi_hpnt.pkl'
#flgmapu  = 'spear2mom_NEP_full_gmapi_upnt.nc'
#flgmapv  = 'spear2mom_NEP_full_gmapi_vpnt.nc'
dflgmaph = os.path.join(dirgmapi, flgmaph)
#dflgmapu = os.path.join(dirgmapi, flgmapu)
#dflgmapv = os.path.join(dirgmapi, flgmapv)
# h-point indices
if os.path.isfile(dflgmaph):
  print(f'Mapping indices hpnt already exist, {dflgmaph} ')
  raise Exception('gmapi h-pnt already exist')

YRI=2010
MMI=1
ens_nmb=1
ifld='siconc'
pthdata = f'/work/Dmitry.Dukhovskoy/tmp/spear_subset/{YRI}/ens{ens_nmb:02d}'
flout = f'NEP_spear_{YRI}{MMI:02d}.{ifld}.nc'
dfspear = os.path.join(pthdata,flout)
ds_spear = xarray.open_dataset(dfspear)
itime=2
A2d = ds_spear[ifld].isel(time=itime).data
LON = ds_spear['GEOLON'].data
LAT = ds_spear['GEOLAT'].data

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
LMsk = np.where(HH>=0, 0, 1)
JJ,II = np.where(LMsk==1)

import mod_regmom as mrmom

print('Searching gmapi for PIOMAS interpolation onto MOM6')
jS = 570
icc = -1
IMOM = II
JMOM = JJ

for kk in range(len(IMOM)):
  if kk%500 == 0:
    print(f' kk={kk} {kk/(len(IMOM))*100:.2f}% done ...')

  jj = JMOM[kk]
  ii = IMOM[kk]
  x0 = hlon[jj,ii]
  y0 = hlat[jj,ii]

  icc += 1
  ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=1.)
  ixx = np.expand_dims(ixx, axis=0)
  jxx = np.expand_dims(jxx, axis=0)

  if kk == 0:
    INDX = ixx.copy()
    JNDX = jxx.copy()
  else:
    INDX = np.append(INDX, ixx, axis=0)
    JNDX = np.append(JNDX, jxx, axis=0)

INDX = INDX.astype(int)
JNDX = JNDX.astype(int)
IMOM = IMOM.astype(int)
JMOM = JMOM.astype(int)

npnts = len(IMOM)
darri = xarray.DataArray(INDX, dims=("npnts","vertices"), \
              coords={"npnts": np.arange(npnts), \
                      "vertices": np.arange(1,5)})
darrj = xarray.DataArray(JNDX, dims=("npnts","vertices"), \
              coords={"npnts": np.arange(npnts), \
                      "vertices": np.arange(1,5)})
dimom = xarray.DataArray(IMOM, dims=("npnts"), coords={"npnts": np.arange(npnts)})
djmom = xarray.DataArray(JMOM, dims=("npnts"), coords={"npnts": np.arange(npnts)})

jdm, idm = LON.shape
darr_lon = xarray.DataArray(LON, dims=("jdim","idim"),\
                  coords={"jdim": np.arange(jdm),\
                          "idim": np.arange(idm)})
darr_lat = xarray.DataArray(LAT, dims=("jdim","idim"),\
                  coords={"jdim": np.arange(jdm),\
                          "idim": np.arange(idm)})

dset = xarray.Dataset({"indx_spear": darri, "jndx_spear": darrj,\
                       "indx_nep": dimom, "jndx_nep": djmom,\
                       "lonh_spear": darr_lon, "lath_spear": darr_lat})

print(f'Saving gmapi --> {dflgmaph}')
dset.to_netcdf(dflgmaph, format='NETCDF3_64BIT', engine='netcdf4')

#print(f'Saving gmapi --> {dflgmaph}')
#with open(dflgmaph, 'wb') as fid:
#  pickle.dump([IMOM, JMOM, INDX, JNDX], fid)


