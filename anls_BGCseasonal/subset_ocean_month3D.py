"""
  Subset specific fields from seasonal forecast standard output archive files
  3D ocean fields

  to run for mulitple years/ ens/ months 
  use script: scripts/seasonal_fcst/
  subset_ocean3D_NEPbgc_fcst.sh

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
VAROCN = ['potT',
          'salt',
          'u',
          'v']

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
ensE = args.ensE if args.ensE else ensS

if MM is None:
  MNTHS = [1,4,7,10]
else:
  MNTHS = [MM]

ENSMB = [x for x in range(ensS,ensE+1)]


print(f'Extracting 3D ocn variables for YR={YRS}-{YRE} MNTHS={MNTHS} ENSMB={ENSMB}')

expt_name = 'NEPbgc_fcst_dailyOB01'
patharch = f'/archive/Dmitry.Dukhovskoy/fre/NEP/forecast_bgc/{expt_name}/'
import glob
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

      fl_pattern = os.path.join(pthfcst, 'oceanm_*_*.nc')
      fl_list = sorted(glob.glob(fl_pattern))

      fl_stamp = f'init{YR}{MM:02d}e{ens:02d}'
      for dflin in fl_list:
        print(f'Processing {dflin}')

        # Extract the file name
        fname = os.path.basename(dflin)  # 'oceanm_1993_07.nc'

        # Parse calendar year and month from file name
        parts = fname.replace('.nc', '').split('_')  # ['oceanm', '1993', '07']
        YR0 = int(parts[1])
        MM0 = int(parts[2])

        flout = f'ocean3Dmonth_{fl_stamp}_{YR0}{MM0:02d}.nc'
        dflout = os.path.join(pathout,flout)

        print(f'subsetting {dflin} --> {dflout}')

        ds = xarray.open_dataset(dflin)
        if VARDROP is None:
          VARDROP = [var for var in ds.data_vars if var not in VAROCN]

        for vars in VARDROP:
          ds = ds.drop_vars(vars)

          ds.attrs["history"] = f"Created from NEP BGC MOM6-SIS2 seasonal f/cast {subdir}"
          ds.attrs["code"] = f"subset_ocean_month3D.py"

        print(f'Saving {dflout} ...')
        ds.to_netcdf(dflout, 
                     format='NETCDF3_64BIT',
                     engine='netcdf4',
                     unlimited_dims='time')
     
        ds.close()

