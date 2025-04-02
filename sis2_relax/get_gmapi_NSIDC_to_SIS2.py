"""
  Derive gmapi indices for bi-linear interpolation
  of NSIDC data onto SIS2 grid
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import pickle
from copy import copy
import matplotlib.colors as colors
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
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')
import mod_swstate as msw
import conversions as gsw

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
import mod_rtofs as mrtofs
importlib.reload(mutob)
importlib.reload(manseas)

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

_, LON, LAT = manseas.avrg_cice_NSIDC(2000, 2000, 1, 1)
# get rid of ice in unneeded part of the domain
#CMobs[:,150:] = np.nan
#CMobs[:200,:] = np.nan

import mod_regmom as mrmom
jS = 570
icc = -1
IMOM = []
JMOM = []
for ii in range(idm):
  if ii%50 == 0:
    print(f' icc={icc} {ii/idm*100:.2f}% done ...')
  for jj in range(jS,jdm):
    if HH[jj,ii] >= 0:
      continue
    x0 = hlon[jj,ii]
    y0 = hlat[jj,ii]
    if y0 < 60.:
      continue
    if y0 < np.min(LAT) or y0 > np.max(LAT):
      continue

    icc += 1
    ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=1.)
    ixx = np.expand_dims(ixx, axis=0)
    jxx = np.expand_dims(jxx, axis=0)

    if icc == 0:
      INDX = ixx.copy()
      JNDX = jxx.copy()
    else:
      INDX = np.append(INDX, ixx, axis=0)
      JNDX = np.append(JNDX, jxx, axis=0)

    IMOM.append(ii)
    JMOM.append(jj)

IMOM = np.array(IMOM)
JMOM = np.array(JMOM)

fgmapi  = f'NSIDC_sis2_NEP_gmapi_{jdm}x{idm}.pkl'
dfgmapi = os.path.join(pthsis, fgmapi)

print(f'Saving gmapi --> {dfgmapi}')
with open(dfgmapi, 'wb') as fid:
  pickle.dump([IMOM, JMOM, INDX, JNDX], fid)





