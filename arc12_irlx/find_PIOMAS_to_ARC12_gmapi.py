"""
  Find gmapi indices: 4 vertices of PIOMAS grid for each ARC12 grid
  for bilinear interpolation of PIOMAS --> ARC12
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

# ARC12 grid:
ptharc  = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'
dflarc  = os.path.join(ptharc,'ocean_hgrid.nc')
dfltopo = os.path.join(ptharc,'ocean_topog.nc')

ds_topo = xarray.open_dataset(dfltopo)
HH = ds_topo['depth'].data
jdm, idm = HH.shape

hlon, hlat = mmom6.read_mom6grid(dflarc, grdpnt='hgrid') 

# PIOMAS LON/LAT:
import mod_regmom as mrmom
fgmapi  = f'PIOMAS_mom6_ARC12_gmapi_{jdm}x{idm}.npz'
pthgmapi = '/work/Dmitry.Dukhovskoy/ARC12/irlx'
dfgmapi = os.path.join(pthgmapi, fgmapi)

pthdata = '/work/Dmitry.Dukhovskoy/data/PIOMAS_ice'
varthck = 'heff'
varconc = 'area'

print('Searching gmapi for PIOMAS interpolation onto MOM6 ARC12')
jS = 0
icc = -1
IMOM = []
JMOM = []
yr0 = 2010
flthck  = f'piomas_heff{yr0}_v21.nc'
flconc  = f'piomas_area{yr0}_v21.nc'
dflthkn = os.path.join(pthdata, flthck)
dflconc = os.path.join(pthdata, flconc)

ds_thkn = xarray.open_dataset(dflthkn)
LAT  = ds_thkn['lat_scaler'].data
LON  = ds_thkn['lon_scaler'].data

# COnvert to -180, -180 grid:
LON = np.where(LON >= 180., LON-360., LON)

tindx = 0
H2d   = ds_thkn[varthck].data[tindx,:].squeeze()  # thikness, m

for ii in range(idm):
  if ii%50 == 0:
    print(f' icc={icc} {ii/idm*100:.2f}% done ...')
  for jj in range(jS,jdm):
    if HH[jj,ii] >= 0:
      continue
    x0 = hlon[jj,ii]
    y0 = hlat[jj,ii]
    if y0 < 50.:
      continue
    if y0 < np.min(LAT) or y0 > np.max(LAT):
      continue

    icc += 1
    ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=1.)
    # Ignore bndry points:
    if len(ixx) == 0 or len(jxx) == 0:
      continue
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

print(f'Saving gmapi --> {dfgmapi}')
np.savez(dfgmapi, IMOM=IMOM, JMOM=JMOM, INDX=INDX, JNDX=JNDX)

#with open(dfgmapi, 'wb') as fid:
#  pickle.dump([IMOM, JMOM, INDX, JNDX], fid)


