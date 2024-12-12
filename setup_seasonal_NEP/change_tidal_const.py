"""
  Check tides imposed at the OB
  OBC_TIDE_CONSTITUENTS = "M2,S2,N2,K2,K1,O1,P1,Q1,MM,MF"
  MF - The near-fortnightly tide Mf, of period 13.66 d, is the largest of  
       the zonally symmetric, long-period tides. Like all the long-period lunar tides, 
       it may be thought of as a time-varying modulation of the Earth permanent tide M0.
       The MF tide is generally 2 cm or less.
  For testing, change amplitudes and U of specific (or all) tidal constituents

  MM  - monthly tides, amplitudes ~0.014 m
"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
#import pickle
import matplotlib.pyplot as plt
from yaml import safe_load

import mod_utils_ob as mutob
importlib.reload(mutob)


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
sys.path.append('./seasonal-workflow')
from boundary import Segment
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
from mod_utils_fig import bottom_text

varnm  = 'ssh'

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)

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
dftopo_mom  = os.path.join(pthtopo, ftopo_mom)

segments = [ Segment(1, 'north', hgrid, output_dir=outdir),
             Segment(2, 'east',  hgrid, output_dir=outdir),
             Segment(3, 'south', hgrid, output_dir=outdir),
             Segment(4, 'west',  hgrid, output_dir=outdir)]

nOB = len(segments)

pthtides='/work/Dmitry.Dukhovskoy/NEP_input/tides_OBs/'

# Any ampl. < 0 - leave unchanged
TAMPL = np.zeros((10)) - 1.
TU    = np.zeros((10)) - 1.
TV    = np.zeros((10)) - 1.
# Remove the fortnight and monthly tides:
TAMPL[8] = 0.0
TAMPL[9] = 0.0
TU[8]    = 0.0
TU[9]    = 0.0
TV[8]    = 0.0
TV[9]    = 0.0

# Amplitudes
for isgm in range(4):
  nsgm = isgm + 1
  const=10
  flnm = f'tz_{nsgm:03d}.nc'
  dfinput = os.path.join(pthtides, flnm)
  dset = xarray.open_dataset(dfinput)
  varnm = f'zamp_segment_{nsgm:03d}'
  AA = dset[varnm].data
  id1, id2, id3, id4 = AA.shape[:4]

  for icnst in range(len(TAMPL)):
    A0 = TAMPL[icnst]
    ampl = AA[0,icnst,:].squeeze()
    ampl = np.where(ampl>0, A0, ampl)
    if A0 >= 0.0:
      if id3 > 1:
        AA[0,icnst,:,0] = ampl
      else:
        AA[0,icnst,0,:] = ampl

  dset = dset.drop_vars(varnm)
  dimy = f'ny_segment_{nsgm:03d}'
  dimx = f'nx_segment_{nsgm:03d}'
  dset = dset.assign(**{varnm: (['time','constituent',dimy, dimx], AA)})

  flnm_out = f'tz_{nsgm:03d}_test.nc'
  flout = os.path.join(pthtides, flnm_out)

  print(f"Writing modified tides --> {flout}")
  dset.to_netcdf(
      flout,
      format='NETCDF3_64BIT',
      engine='netcdf4',
      encoding={'time': {'dtype': 'float64', 'calendar': 'gregorian'}},
      unlimited_dims='time'
  )

# Velocities
for isgm in range(4):
  nsgm = isgm + 1
  const=10
  flnm = f'tu_{nsgm:03d}.nc'
  dfinput = os.path.join(pthtides, flnm)
  dset = xarray.open_dataset(dfinput)
  unm = f'uamp_segment_{nsgm:03d}'
  vnm = f'vamp_segment_{nsgm:03d}'
  UU  = dset[unm].data
  VV  = dset[vnm].data
  id1, id2, id3, id4 = UU.shape[:4]

  for icnst in range(len(TAMPL)):
    U0 = TU[icnst]
    V0 = TV[icnst]
    ampl = UU[0,icnst,:].squeeze()
    ampl = np.where(ampl>0, U0, ampl)
    if U0 >= 0.0:
      if id3 > 1:
        UU[0,icnst,:,0] = ampl
      else:
        UU[0,icnst,0,:] = ampl

    ampl = VV[0,icnst,:].squeeze()
    ampl = np.where(ampl>0, V0, ampl)
    if V0 >= 0.0:
      if id3 > 1:
        VV[0,icnst,:,0] = ampl
      else:
        VV[0,icnst,0,:] = ampl


  dset = dset.drop_vars(unm)
  dimy = f'ny_segment_{nsgm:03d}'
  dimx = f'nx_segment_{nsgm:03d}'
  dset = dset.assign(**{unm: (['time','constituent',dimy, dimx], UU)})

  dset = dset.drop_vars(vnm)
  dset = dset.assign(**{vnm: (['time','constituent',dimy, dimx], VV)})

  flnm_out = f'tz_{nsgm:03d}_test.nc'
  flout = os.path.join(pthtides, flnm_out)

  print(f"Writing modified tides --> {flout}")
  dset.to_netcdf(
      flout,
      format='NETCDF3_64BIT',
      engine='netcdf4',
      encoding={'time': {'dtype': 'float64', 'calendar': 'gregorian'}},
      unlimited_dims='time'
  )


print("All Done")


