"""

  Code is based on the code by Andrew C. Ross
  Modified and edited

  I ran in python for several years at a time, didn't use sbatch due to time limit (1 hr)

  From Andrew:
  the forecasts are forced from the atmosphere by different ensemble members 
  of GFDL's SPEAR seasonal prediction system. The SPEAR data is stored on 
  GFDL's archive at /archive/l1j/spear_med/rf_hist/fcst/s_j11_OTA_IceAtmRes_L33. 
  Within my seasonal workflow repository, the script write_spear_atmos.py will 
  pull the data for a forecast from archive and format it for use in MOM6. 
  There is also a .sh file that you can sbatch instead to run multiple setup 
  jobs at once (though the filesystem sometimes fails if you have more than 
  a few running at one time). 

  ens: ensemble member; either an integer, to get a single member, or 
       "pp_ensemble" to get the post-processed ensemble mean

"""
import numpy as np
import os
from pathlib import Path
from spear_path import get_spear_paths
import subprocess
import xarray
import argparse
from yaml import safe_load

from utils import pad_ds, write_ds

# To run this code cdo is required
#subprocess.run(['module load cdo'], shell=True)
try:
  subprocess.run(["which cdo"], shell=True, check=True, capture_output=True) 
except subprocess.CalledProcessError as err:
  print("Need to load module cdo prior to python session")
  print(err)
  raise Exception("cdo module is missing")


# Specify year, month, ensembles to 
# generate atm fields from SPEAR 
yr1     = 2006
yr2     = 2009
MM      = [1,4,7,10]
ensmb   = [1,2,3,4,5,6,7,8,9,10]
fconfig = 'config_nep.yaml'
#fconfig = 'config_nwa12.yaml'

print(f"Creating SPEAR atmos fields for NEP, {yr1}-{yr2}")
print(f"N ensembles: {len(ensmb)}, Months: {MM}")
with open(fconfig) as ff:
  config = safe_load(ff)
# Note conversion from [-180, 180] to [0, 360]
lon_slice = slice(config['domain']['west_lon'] % 360, config['domain']['east_lon'] % 360)
lat_slice = slice(config['domain']['south_lat'], config['domain']['north_lat'])
work_dir = Path(config['filesystem']['model_input_data']) / 'atmos'
work_dir.mkdir(exist_ok=True)


def extract_spear_atmos(ystart, mstart, ens, outdir=None):
  if outdir is None:
    tmp = Path(os.environ['TMPDIR'])
    if ens == 'pp_ensemble':
        outdir = tmp / f'{ystart}-{mstart:02d}-{ens}_raw'
    else:
        outdir = tmp / f'{ystart}-{mstart:02d}-e{ens:02d}_raw'
    outdir.mkdir(exist_ok=True)

  print(f'outdir = {tmp}')
  files = get_spear_paths(
      ['slp', 't_ref', 'u_ref', 'v_ref', 'q_ref', 'lwdn_sfc', 'swdn_sfc', 'precip'],
      ystart, mstart, 'atmos_daily', 'daily', ens=ens
  )
# Having troubles with subprocess using the whole list of files etc
#  file_strings = ' '.join(map(lambda x: x.as_posix(), files))
#  print(f'Unstaging {len(files)} SPEAR files from tape, may take a while ...')
#  subprocess.run(['dmget ' + file_strings], shell=True, check=True)
#  subprocess.run(['gcp --sync ' + file_strings + ' ' + outdir.as_posix()], shell=True, check=True)
  nn = len(files)
  for ff in files:
    file_spear = ff.as_posix()
    print(f"Unstaging {file_spear} ...")
    subprocess.run(["dmget", file_spear])
    print(f"Copying {file_spear} --> {outdir.as_posix()}")
    subprocess.run(["gcp","--sync",file_spear,outdir.as_posix()], check=True)
    subprocess.run([f'gcp --sync {file_spear} {outdir.as_posix()}'], shell=True) 

  new_files = [outdir / f.name for f in files]

  return new_files

    
for ystart in range(yr1,yr2+1):
  for mstart in MM:
    for ens in ensmb:
      print(f'Processing {ystart} {mstart} ens={ens}')
#      write_atmos(yr0, mm, ens0, work_dir, yslice, xslice)

      out_dir = work_dir / f'{ystart}-{mstart:02d}-e{ens:02d}'
      if not out_dir.is_dir():
        # Read mask for flooding
        static = xarray.open_dataset('/work/acr/spear/atmos.static.nc')
        is_ocean = np.invert(static.land_mask.astype('bool'))

        print(f'extracting SPEAR atmos {ystart}, {mstart}, {ens}')
        extracted_files = extract_spear_atmos(ystart, mstart, ens)

        print('pad and save')
        tmpdir = Path(os.environ['TMPDIR']) / 'atmos_raw' / f'{ystart}-{mstart:02d}-e{ens:02d}'
        tmpdir.mkdir(exist_ok=True, parents=True)
        for f in extracted_files:
          print(f'{str(f)}')
          #open
          ds = xarray.open_dataset(f).sel(lat=lat_slice, lon=lon_slice)

          # Need to mask just the variable of interest and not the
          # coordinate/metadata variables 
          main_var = [x for x in ds.data_vars if len(ds[x].dims)==3][0]
          ds[main_var] = ds[main_var].where(is_ocean)
          padded = pad_ds(ds)
          fout = tmpdir / f.name
          write_ds(padded, fout)
          # cdo doesn't like if the input is also the output here
# setmisstodis: Set missing value to the distance-weighted average i
#               of the nearest neighbors

          subprocess.run([f'cdo -O replace {fout} -setmisstodis,3 -selvar,{main_var} {fout} {fout}.new'], shell=True, check=True)
          # Rename the file generated by cdo to the output file 
          fout.with_suffix(fout.suffix + '.new').rename(fout)

        print(f'gcp atmos*nc --> {out_dir}')
        out_dir.mkdir(exist_ok=True)
        subprocess.run([f'gcp {tmpdir}/atmos*.nc {out_dir}'], shell=True, check=True)
      else:
        print(f'Already found data for {out_dir.as_posix()}')







