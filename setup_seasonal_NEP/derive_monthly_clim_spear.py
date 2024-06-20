"""
  after extracting domain from SPEAR monthly fields
  calculate monthly climatology 

  For the climatology, 1 ensemble member started at 1 month is used for now

"""
import numpy as np
import os
from pathlib import Path
import importlib
from glob import glob
import spear_path
from spear_path import get_spear_paths
import subprocess
import xarray
import argparse
from yaml import safe_load
from pathlib import Path, PurePath
from functools import partial

# Average over all existing files
#YR1     = 1993
#YR2     = 2006
mstart  = 1
ens     = 1
regn    = 'NEP'
fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)

# Find all years:
pthtmp = config['filesystem']['spear_month_ens']
fyrs   = glob(os.path.join(f'{pthtmp}/*/ens{ens:02d}'))
YRS    = []
for ff in fyrs:
  yr = int(ff.split("/")[-2])
  YRS.append(yr)
YR1 = min(YRS)
YR2 = max(YRS)


print(f"Computing monthly climatology from SPEAR for {YR1}-{YR2}")
def find_tmpdir(yr, ens, mstart):
  pthtmp = config['filesystem']['spear_month_ens']
  dirtmp = PurePath(pthtmp) / f'{yr}' / f'ens{ens:02d}' 
  return dirtmp

#YRS = [x for x in range(YR1, YR2+1)]
fun = partial(find_tmpdir, ens=ens, mstart=mstart)
PATHS = [fun(yy) for yy in YRS]

dir_outp = Path(config['filesystem']['spear_month_ens']) / 'monthly_clim'
dir_outp.mkdir(exist_ok=True)    

VARS = ['SSH','so','thetao','vo','uo']

for varnm in VARS:
  if varnm == 'SSH':
    fspear_in = glob(os.path.join(f'{pthtmp}/*/ens{ens:02d}/NEP_spear_*{mstart:02d}.ssh.nc'))
  else:
    fspear_in = glob(os.path.join(f'{pthtmp}/*/ens{ens:02d}/NEP_spear_*{mstart:02d}.{varnm}.nc'))
#  dset = xarray.open_mfdataset(fspear_in, combine='nested', concat_dim='time') 

  assert len(fspear_in) > 0, "no files found"

  print(f"Processing {varnm}") 
  icc = 0 
  for flnm in fspear_in:
    dset = xarray.open_dataset(flnm)
    dset = dset.drop_vars(['average_DT','average_T1','average_T2','time_bnds'])

    if icc == 0: 
      dset_mean = dset.copy()
      Asum = dset[f'{varnm}'].data
    else:
      Asum = Asum + dset[f'{varnm}'].data

    icc += 1

  Asum = Asum / icc
  dset[f'{varnm}'].values = Asum

  foutp = os.path.join(dir_outp,f'spear_moclim.mstart{mstart:02d}.{YR1}-{YR2}.{varnm}.nc')
  if varnm == 'SSH':
    dset = dset.rename({'SSH': 'ssh'})
    foutp = os.path.join(dir_outp,f'spear_moclim.mstart{mstart:02d}.{YR1}-{YR2}.ssh.nc')

  print(f'Saving --> {foutp}')
  dset.to_netcdf(foutp)


