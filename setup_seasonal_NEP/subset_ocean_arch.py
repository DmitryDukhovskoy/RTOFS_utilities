"""
  Subset specific fields from seasonal forecast standard output archive files
  tos, zos, tauu, tauv, uv & uv,

  All arguments are optional 
  usage:  run subset_ocean_arch.py --yr=1993 --mo=4 --ens=10

"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import argparse

MNTHS = [1,4,7,10]
YR = 2012
ENSMB = [x for x in range(1,11)]
VARDROP = ['sshmin', 
           'sshmax',
           'tosmin',
           'tosmax',
           'tossq',
           'sossq',
           'sos',
           'sob',
           'omldamax',
           'average_T1',
           'average_T2',
           'average_DT',
           'time_bnds']

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="Forecast initialization year: 1993, ..., 2020", type=int)
parser.add_argument("--mo", help="Forecast initialization month: 1,4,7,10", type=int)
parser.add_argument("--ens", help="Ensemble number: 1, ..., 10", type=int)
args = parser.parse_args()

if args.yr:
  if args.yr > 1990:
    YR = args.yr
  else:
    raise Exception(f"year should be > 1990")

if args.mo: 
  if args.mo > 0 and args.mo <= 12:
    MNTHS = [args.mo]

if args.ens:
  if args.ens > 0 and args.ens <= 10:
    ENSMB = [args.ens]

print(f'Extracting variables for YR={YR} MNTHS={MNTHS} ENSMB={ENSMB}')

expt_name = 'NEPphys_frcst_dailyOB-expt02'
patharch = f'/archive/Dmitry.Dukhovskoy/fre/NEP/seasonal_daily/{expt_name}/'
pathout  = f'/collab1/data_untrusted/Dmitry.Dukhovskoy/{expt_name}/{YR}'

if not os.path.exists(pathout):
    print(f'Creating {pathout}')
    os.makedirs(pathout)

for MM in MNTHS:
  if YR == 1993 and MM == 1:
    print(f"Requested MM={MM}, F/casts start in MM=4 in 1993, no f/casts for MM=1, skipping ...")
    continue
  if YR ==2020 and MM > 1:
    print(f"Requested MM={MM}, F/casts start only in MM=1 in 2020, no f/casts for MM>1, skipping ...")
    continue

  for ens in ENSMB:
    subdir = f'{YR}-{MM:02d}-e{ens:02d}'
    pthfcst = os.path.join(patharch, subdir, 'history')
    dflin = os.path.join(pthfcst,'ocean_daily.nc')
    flout = f'ocean_dailysub_{subdir}.nc'
    dflout = os.path.join(pathout,flout)

    print(f'subsetting {dflin} --> {dflout}')

    ds = xarray.open_dataset(dflin)
    for vars in VARDROP:
      ds = ds.drop_vars(vars)

      ds.attrs["history"] = f"Created from NEP MOM6-SIS2 seasonal f/cast {subdir}"
      ds.attrs["code"] = f"/home/Dmitry.Dukhovskoy/python/setup_seasonal_NEP/subset_ocean_arch.py"

    print(f'Saving {dflout} ...')
    ds.to_netcdf(dflout, 
                 format='NETCDF3_64BIT',
                 engine='netcdf4',
                 unlimited_dims='time')
 

