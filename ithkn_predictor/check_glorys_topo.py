"""
Version 2

  Check static predictor: GLORYS bottom depth at training grid points
  GLORYS

  Use time stamps and J,I grid points from the response
  variable (ithkn), that need to be run first

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray as xr
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, cm
from yaml import safe_load
import argparse
from pathlib import Path

#ROOT = Path(__file__).resolve().parent

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
import mod_glorys as mglr 
import mod_icepredict as micpr

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="Region to process", choices=['north','south'], 
                    required=True, type=str)
args = parser.parse_args()

regn       = args.regn
fld_name   = 'bottom_depth'
varnm      = "deptho"

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

# Load parameters:
regn_name    = config_predictor["regn"][regn]["name"]
lat0         = config_predictor["regn"][regn]["lat_bnd"]
tstep_era    = config_predictor["params"]["tstep"]     # ERA5 daily fields, skipping days
dxy          = config_predictor["params"]["dxy"]
YS           = config_predictor["params"]["ys"]
YS           = config_predictor["params"]["ys"]
YE           = config_predictor["params"]["ye"]
davrg        = config_predictor["params"]["davrg"]
sqrt_frzdays = config_predictor["params"]["sqrt_frz"]   # use sqrt(IFDD) to better fit Zubov relation
intgr_time   = config_predictor["params"]["intgr_frz"]  # Time for freeze degree days accumulation
Tfrz         = config_predictor["params"]["Tfrz"]       # ocean freezing T


DIRS = {
  "pthout"   : config_predictor["linregr"]["pthout"],
  "pthtopo"  : config_predictor["linregr"]["pthgltopo"],
  "gltopo"   : config_predictor["linregr"]["gltopo"],
  }

# Get time steps, grid points:
pthinfo = config_predictor["params"]["pthinfo"]
fltime = config_predictor["params"]["fltime"].format(regn=regn, tstep=tstep_era)
flij   = config_predictor["params"]["flij"].format(regn=regn)
flgrid = config_predictor["params"]["flgrid"]
dfltime = os.path.join(pthinfo, fltime)
dflij   = os.path.join(pthinfo, flij)
dflgrid = os.path.join(pthinfo, flgrid)

assert os.path.isfile(dfltime), f"Time steps file is missing: {dfltime}"
assert os.path.isfile(dflij), f"Subsample grid points file is missing: {dflij}"

# Read saved date numbers and 
# Read time array ad J,I sample grid points
DNMB = np.load(dfltime)
A = np.load(dflij)
JG = A["JG"]
IG = A["IG"]
npnts = IG.shape
nrecs = len(DNMB)

print(f"N of grid points: {npnts}, N time steps: {nrecs}")

# GLORYS grid:
A = np.load(dflgrid)
hlon = A["LON"]
hlat = A["LAT"]
LMsk = A["LMsk"]

# 0 <= lon < 360
hlon = (hlon + 360) % 360

pthgltopo = config_predictor["linregr"]["pthgltopo"]
flgltopo  = config_predictor["linregr"]["gltopo"]
dgltopo = os.path.join(pthgltopo, flgltopo)

print(f"Reading GLORYS topo from {dgltopo}")
assert os.path.isfile(dgltopo), f"GLORYS topo not found {dgltopo}"
with xr.open_dataset(dgltopo) as dstopo:
  A2d = dstopo[varnm].values.squeeze()

# Convert to negative depths:
if np.nanmin(A2d) > 0:
  HH = -A2d

topo_pnts = HH[JG, IG]

# Normalize: deepest ocean = -1, land = +1.
LMocn = np.isfinite(topo_pnts)

topo_norm = np.ones_like(topo_pnts, dtype=float)
topo_min = np.nanmin(topo_pnts)  # Most negative depth

topo_norm[LMocn] = topo_pnts[LMocn] / abs(topo_min)


f_check = False
if f_check:
  plt.ion()

  rmin = -1
  rmax = 1.

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
  ax1.contour(LMsk, [0.99], linestyles='solid', colors=[(0.4,0.4,0.4)])

  #ax1.pcolormesh(A2d, cmap='jet', vmin=rmin, vmax=rmax)

  sc = ax1.scatter(
    IG, JG,
    c=topo_norm,
    cmap='jet',
    s=20,          # marker size
    vmin=rmin,        # optional color scale limits
    vmax=rmax
 )

  jdim, idim = LMsk.shape
  ax1.axis('scaled')
  ax1.set_ylim(jdim // 2, jdim)

  ax1.set_title(f"Grid points for ML training\n region={regn} N pnts={len(JG)}")


  plt.colorbar(sc, ax=ax1, label='Depth, m', orientation='horizontal')

  btx = f' @{machine}: check_glorys_topo.py'
  bottom_text(btx, pos=[0.1,0.1], fsz=10)


