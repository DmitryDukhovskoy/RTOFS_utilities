"""
Combine north / south interpolated NSIDC ice concenatrtion fields
created by interp_NSIDC_iconc_mesh025.py
"""
import os
import numpy as np
import sys
import importlib
import xarray as xr
from yaml import safe_load
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


import mod_time as mtime
import mod_misc1 as mmisc
import mod_mom6 as mmom6
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)


parser = argparse.ArgumentParser()
parser.add_argument("--yrs", help=f"Start year of NSIDC data", required=True, type=int)
parser.add_argument("--yre", help=f"End year of NSIDC data, default=start year", type=int)
parser.add_argument("--ms", help=f"Start month of NSIDCS data to combine", required=True, type=int)
parser.add_argument("--me", help=f"End month of NSIDCS data to combine, default=start",  type=int)
parser.add_argument("--clean", 
  help="=1: After combining regional files, remove original files, default=0", 
  choices=[0,1],
  default=0,
  type=int)
args = parser.parse_args()

YRS = args.yrs
YRE = args.yre if args.yre is not None else args.yrs
MS = args.ms
ME = args.me if args.me is not None else args.ms
rm_orig = args.clean == 1

assert ME >= MS, f"End month {ME} must be >= {MS}"

if rm_orig:
  print("\n    Original regional NSIDC files will be removed after being combined\n")

pthr = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/NRT_NOAA_NSIDC_seaconc/'


def combine_month(pth, fln, fls, flg, LON, LAT, remove_original=rm_orig):
  fl_north = os.path.join(pth, fln)
  fl_south = os.path.join(pth, fls)
  fl_global = os.path.join(pth, flg)

  # Do not recreate existing global file
  if os.path.exists(fl_global):
    print(f"Global file already exists: \n{fl_global}, skipping ...")
    return fl_global

  # Check input files
  for fl in (fl_north, fl_south):
    if not os.path.isfile(fl):
      #raise FileNotFoundError(fl)
      print(f"Does not exist: {fl}")

  # Combine fields
  varnm = 'ice_conc'
  with xr.open_dataset(fl_north) as dsn, \
       xr.open_dataset(fl_south) as dss:

    An = dsn[varnm].values
    As = dss[varnm].values

    assert An.shape == As.shape, "North/south shapes differ"
    assert dsn[varnm].dims == dss[varnm].dims, \
        "North/south dimensions differ"

    # Grid dimensions
    jdname, idname = dsn[varnm].dims[-2:]
    jdim, idim = An.shape[-2:]

    assert LON.shape == (jdim, idim)
    assert LAT.shape == (jdim, idim)

    # Combine fields
    # LAT > 0: north; LAT <= 0: south
    north_mask = LAT > 0
    Ag = np.where(north_mask[None, :, :], An, As)

    # Check number of ice grid points
    Nice_north = np.count_nonzero(
        (An > 0) & np.isfinite(An) & north_mask[None, :, :]
    )

    Nice_south = np.count_nonzero(
        (As > 0) & np.isfinite(As) & ~north_mask[None, :, :]
    )

    Nice_global = np.count_nonzero((Ag > 0) & np.isfinite(Ag))

    assert Nice_global == Nice_north + Nice_south, \
        "Number of ice grid points not preserved"

    print(f"Ice points: N={Nice_north}, S={Nice_south}, "
          f"Global={Nice_global}")

    # Create output dataset
    dsout = dsn.copy(deep=True)
    dsout[varnm].values[:] = Ag

    # Add lon/lat single precision
    dsout["longit"] = xr.DataArray(
                LON.astype(np.float32),
                dims=(jdname, idname)
            )

    dsout["latit"] = xr.DataArray(
                LAT.astype(np.float32),
                dims=(jdname, idname)
            )

    dsout["longit"].attrs["long_name"] = "Longitudes of mesh025"
    dsout["latit"].attrs["long_name"] = "Latitudes of mesh025"
    dsout.attrs.update({
        "title": "NRT NSIDC v6 sea ice conc daily interpolated "
                 "onto mesh025 grid, N and S polar regions",
        "source": "combine_interpNSIDC_regn2global_mesh025.py",
        "region": "global",
    })
    # Delete unneeded attributes:
    dsout.attrs.pop("institution", None)
    dsout.attrs.pop("contact", None)
    dsout.attrs.pop("Grid_idm_jdm", None)


  # Write netcdf:
  print(f"Writing combined fields --> {fl_global}")
  dsout.to_netcdf(fl_global, format="NETCDF4", engine='netcdf4')
  dsout.close()

  # Verify output after closing input datasets
  with xr.open_dataset(fl_global) as dscheck:
    assert dscheck[varnm].shape == An.shape

    # Verify the actual combined field
    np.testing.assert_array_equal(dscheck[varnm].values, Ag)

    
  # Remove original files only after successful output
  if remove_original:
      os.remove(fl_north)
      os.remove(fl_south)
      print("Original files removed ...")

  return fl_global

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

# Get MOM6 grid
pthgrid = pths_ufs["gaea"]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')


fltmp = "NSIDC_iconc_interp_mesh025_1080x1440_"
for YR in range(YRS, YRE+1):
  pthiconc =os.path.join(pthr,f"{YR}")
  for MM in range(MS, ME+1):
    fln = f"{fltmp}{YR}{MM:02d}_north.nc"
    fls = f"{fltmp}{YR}{MM:02d}_south.nc"
    flg = f"{fltmp}{YR}{MM:02d}_global.nc"
    
    print(f"Processing {YR}/{MM}")
    fl_global = combine_month(pthiconc, fln, fls, flg, hlon, hlat, remove_original=rm_orig)


print("    All Done\n")

