# Calculate monthly mean potential densities
# at z= z0 depth
# NEP MOM6 simulation
#
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#import torch
import sys
import pdb
import netCDF4
import importlib
from netCDF4 import Dataset as ncFile
import timeit
import pickle
import yaml


PPTHN = []
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
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mom6util
import mod_misc1 as mmsc1

YR       = 1993
pfld     = 'ssh'
zz0      = -200.
f_derive = False
expt     = 'NEP_BGCphys'
pthpkl   = '/work/Dmitry.Dukhovskoy/data/mom6_nep_tmp/'
fldump   = f"{pthpkl}{expt}_meanSSH_{YR}.pkl"

btx = 'hovm_dltssh.py'

with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom)
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom)
LONM, LATM  = mom6util.read_mom6grid(dfgrid_mom)
HHM         = mom6util.read_mom6depth(dftopo_mom)
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]

#DX, DY = mhycom.dx_dy(LONM, LATM)

import mod_swstate as msw
Rho_mn = []
icc    = 0
for jday in range(1,366):
  dnmb   = mtime.jday2dnmb(YR, jday)
  dv     = mtime.datevec(dnmb)
  ds     = mtime.datestr(dnmb)
  MM     = dv[1]
  DM     = dv[2]

  print(f"Reading {YR}/{MM}/{DM}")
  pthout = f'/work/Dmitry.Dukhovskoy/run_output/NEP_BGCphys/{YR}/{MM:02d}/'
  flmom  = f'ocean_{YR}_{jday:03d}_12.nc'
  dflmom = os.path.join(pthout,flmom)
  if not os.path.isfile(dflmom):
    print(f"File not found skipping: {dflmom}")
    continue

  nc   = ncFile(dflmom,'r')
  if icc == 0:
    ZL = nc.variables['zl'][:].data.squeeze()
    if np.min(ZL) > 0.: ZL = -ZL
# Find closest depth
    dmm = np.abs(ZL-zz0)
    iz0 = np.argmin(dmm)

#  dP = nc.variables['h'][:].data.squeeze()
#  ssh = nc.variables['ssh'][:].data.squeeze()
  T2d = nc.variables['potT'][iz0,:,:].data.squeeze()
  S2d = nc.variables['salt'][iz0,:,:].data.squeeze()
  Rho2d = msw.sw_denso0(S2d, T2d)
 


