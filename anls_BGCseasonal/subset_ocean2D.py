"""
  Subset specific fields from seasonal forecast standard output archive files
  2D ocean fields

  to run for mulitple years/ ens/ months 
  use script: scripts/seasonal_fcst/
  subset2D_NEPbgc_fcst.sh

  Do one field at a time:
  sbatch subset2D_NEPbgc_fcst.sh --yrs 2001 --yre 2010 --ocnm t
  sbatch subset2D_NEPbgc_fcst.sh --yrs 2011 --yre 2024 --icem t
  sbatch subset2D_NEPbgc_fcst.sh --yrs 2011 --yre 2024 --cobbtm t
  sbatch subset2D_NEPbgc_fcst.sh --yrs 2011 --yre 2024 --cobtrc t

Here are the priorities:

First priority: 2D files
ocean_month: zos, tos, tob, sos, sob, tauuo, tauvo, ustar, MLD_002, MLD_restrat (is mlots available, or only mlotsmin/mlotsmax?)
ice_month: siconc
ocean_cobalt_btm: btm_o2, btm_htotal, btm_co3_ion, btm_co3_sol_arag
ocean_cobalt_tracers_int: nsmp_100, nmdp_100, nlgp_100, nsmz_100, nmdz_100, nlgz_100

Second priority: 3D physics
oceanm_yyyy_mm: potT, salt, u, v

Third priority: 3D biogeochemistry
ocean_cobalt_tracers_month_z: htotal, no3, o2, chl, nlg, nmd, nsm, nsmz, nmdz, nlgz

"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import argparse

# List of ocn variables to keep
VAROCN = ['zos',
           'tos',
           'tob',
           'sos',
           'sob',
           'tauuo',
           'tauvo',
           'ustar',
           'MLD_003',
           'MLD_restrat',
           'mlotsmin',
           'mlotsmax']

# List of variables to drop !!!
VARDROP = None

parser = argparse.ArgumentParser()
parser.add_argument("--yrs", help="Start: Forecast initialization year: 1993, ..., 2020", type=int, required=True)
parser.add_argument("--yre", help="End: Forecast initialization year: 1993, ..., 2020", type=int)
parser.add_argument("--mo", help="Forecast initialization month: 1,4,7,10, default: all", type=int)
parser.add_argument("--ensS", help="Start: Ensemble number: 1, ..., 10, default: all", type=int)
parser.add_argument("--ensE", help="End: Ensemble number: 1, ..., 10, default: all", type=int)
args = parser.parse_args()

YRS = args.yrs if args.yrs else None
YRE = args.yre if args.yre else YRS
MM  = args.mo if args.mo else None
ensS = args.ensS if args.ensS else 1
ensE = args.ensE if args.ensE else 10

if MM is None:
  MNTHS = [1,4,7,10]
else:
  MNTHS = [MM]

ENSMB = [x for x in range(ensS,ensE+1)]


print(f'Extracting 2D ocn variables for YR={YRS}-{YRE} MNTHS={MNTHS} ENSMB={ENSMB}')

expt_name = 'NEPbgc_fcst_dailyOB01'
patharch = f'/archive/Dmitry.Dukhovskoy/fre/NEP/forecast_bgc/{expt_name}/'

for YR in range(YRS,YRE+1):
  #pathout  = f'/collab1/data_untrusted/Dmitry.Dukhovskoy/{expt_name}/{YR}'
  pathout = f'/work/Dmitry.Dukhovskoy/tmp/{expt_name}/{YR}'
  if not os.path.exists(pathout):
      print(f'Creating {pathout}')
      os.makedirs(pathout)

  for MM in MNTHS:
    if YR == 1993 and MM == 1:
      print(f"Requested MM={MM}, F/casts start in MM=4 in 1993, no f/casts for MM=1, skipping ...")
      continue
    if YR == 2025 and MM >= 1:
      print(f"Requested MM={MM}, F/casts end in 10/2024, skipping ...")
      continue

    for ens in ENSMB:
      subdir = f'{YR}-{MM:02d}-e{ens:02d}'
      pthfcst = os.path.join(patharch, subdir, 'history')
      dflin = os.path.join(pthfcst,'ocean_month.nc')
      flout = f'ocean2D_month_{subdir}.nc'
      dflout = os.path.join(pathout,flout)

      print(f'subsetting {dflin} --> {dflout}')

      ds = xarray.open_dataset(dflin)
      if VARDROP is None:
        VARDROP = [var for var in ds.data_vars if var not in VAROCN]

      for vars in VARDROP:
        ds = ds.drop_vars(vars)

        ds.attrs["history"] = f"Created from NEP BGC MOM6-SIS2 seasonal f/cast {subdir}"
        ds.attrs["code"] = f"/home/Dmitry.Dukhovskoy/python/setup_seasonal_NEP/subset_ocean3Dflds.py"

      print(f'Saving {dflout} ...')
      ds.to_netcdf(dflout, 
                   format='NETCDF3_64BIT',
                   engine='netcdf4',
                   unlimited_dims='time')
   
      ds.close()

