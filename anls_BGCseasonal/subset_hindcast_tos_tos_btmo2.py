"""
  Subset specific fields from seasonal forecast standard output archive files
  2D ocean fields

Here are the priorities:

First priority: 2D files
ocean_month: zos, tos, tob, sos, sob, tauuo, tauvo, ustar, MLD_002, MLD_restrat (is mlots available, or only mlotsmin/mlotsmax?)
ice_month: siconc
ocean_cobalt_btm: btm_o2, btm_htotal, btm_co3_ion, btm_co3_sol_arag
ocean_cobalt_tracers_int: nsmp_100, nmdp_100, nlgp_100, nsmz_100, nmdz_100, nlgz_100

Note: zos is ssh above geoid, ssh is ssh above the MSL (=0)


Second priority: 3D physics
oceanm_yyyy_mm: potT, salt, u, v

Third priority: 3D biogeochemistry
ocean_cobalt_tracers_month_z: htotal, no3, o2, chl, nlg, nmd, nsm, nsmz, nmdz, nlgz

"""
import datetime as dt
import numpy as np
import xarray
import os
import importlib
import sys
import argparse

# Append custom module paths
PPTHN = None
if 'PPTHN' not in locals() or PPTHN is None:
  cwd = os.getcwd()
  parts = cwd.split(os.sep)
  if 'python' in parts:
    idx = parts.index('python')
    PPTHN = os.sep + os.path.join(*parts[:idx + 1])
  else:
    raise RuntimeError("Directory 'python' not found in current working directory path.")

sys.path.extend([
    os.path.join(PPTHN, 'MyPython', 'hycom_utils'),
    os.path.join(PPTHN, 'MyPython', 'draw_map'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])


from mod_utils_fig import bottom_text
import mod_time as mtime


# List of ocn variables to keep
VAROCN = [
         'tos',
         'tob',
         ]

# List of variables to drop
VARDROP = None

# History files: saved 3 months forward, e.g. 20100101, saved time: Jan 15, Feb 15, Mar 15 2010
mnths_outp  = 3  # how many months saved in 1 history file
INIT_MO = np.array([1,4,7,10])  # init months

parser = argparse.ArgumentParser()
parser.add_argument(
    "--years",
    help="Years of h/cast run to extract: 1993, ..., 2024",
    type=int,
    nargs="+",
    required=True
)
parser.add_argument(
    "--mm",
    help="Months of h/cast, list; default=1,...,12",
    type=int,
    nargs="+",
)
args = parser.parse_args()

YEARS = args.years
MNTHS = args.mm if args.mm else [x for x in range(1,13)]
print(f'Extracting 2D ocn variables for YR={YEARS[0]}-{YEARS[-1]} MNTHS={MNTHS}')

expt_name = 'NEPbgc_hindcast'
patharch = '/archive/Dmitry.Dukhovskoy/fre/NEP/hindcast_bgc/NEPbgc_nudged_hindcast02/history'

for YR in YEARS:
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

    # Find corresponding 3-mo interval where MM belongs:
    indx_init = (MM - 1) // mnths_outp
    MMI = INIT_MO[indx_init]
    indx_mo = (MM - 1) % mnths_outp
    print(f"/nYR={YR} MM={MM} INIT MONTH={MMI} indx_mo={indx_mo}")

    subdir = f'{YR}{MMI:02d}01'
    pthfcst = os.path.join(patharch, subdir)
    dflbtmo2 = os.path.join(pthfcst,'ocean_cobalt_btm.nc')
    dfltocn = os.path.join(pthfcst,'ocean_month.nc')
    flout = f'BGChcast_ocean2D_month_{YR}{MM:02d}.nc'
    dflout = os.path.join(pathout,flout)

    print(f'Reading btm_o2 {dflbtmo2}')
    # Keep as data array and keep the metadata (do not use *.data):
    with xarray.open_dataset(dflbtmo2) as ds_btm:
      BTMO2 = ds_btm["btm_o2"]  # keep as DataArray

    print(f'subsetting {dfltocn} --> {dflout}')


    ds = xarray.open_dataset(dfltocn)
    VARDROP = [var for var in ds.data_vars if var not in VAROCN]
    ds = ds.drop_vars(VARDROP)

    # Select correct month but keep time dimension
    ds = ds.isel(time=slice(indx_mo, indx_mo + 1))
    BTMO2 = BTMO2.isel(time=slice(indx_mo, indx_mo + 1))

    # Add BTMO2:
    # Align just in case
    BTMO2 = BTMO2.sel(time=ds.time)
    ds["btm_o2"] = BTMO2

   
    ds.attrs["history"] = f"Ocean surf/btm T and btm O2 from NEP BGC MOM6-SIS2 hindcast {subdir}"
    ds.attrs["code"] = "subset_hindcast_tos_tos_btmo2.py"

    print(f'Saving {dflout} ...')
    ds.to_netcdf(dflout, 
                 format='NETCDF3_64BIT',
                 engine='netcdf4',
                 unlimited_dims='time')
 
    ds.close()

