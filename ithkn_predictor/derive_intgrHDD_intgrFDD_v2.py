"""
Version 2.:
  Calculate integrated Heat Degree Days or Freeze Degree Days
  Use hourly ERA5 --> daily

Code for both regions Arctic / Antarctic
  Grid points and time steps are prepared in define_time_IJpnts.py
  Grid points with no ice cover during multiple years are eliminated

  Heat degree days
  in analogy to Zubov's IFDD, but opposite
  which may be important for melt season

  For faster processing, run unstaging script before this:
  /home/Dmitry.Dukhovskoy/scripts/GLORYS_anls/unstage_glorys.sh
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
from scipy.interpolate import interp1d

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
from mod_misc1 import dist_sphcrd
from mod_mom6 import dx_dy
import mod_icepredict as micepr

#from MyPython.mod_cice6_utils import change_base_template, flname_replace_date

parser = argparse.ArgumentParser()
parser.add_argument(
  "--regn", 
  help="Region to process", 
  choices=['north','south'], 
  required=True
)
parser.add_argument(
  "--field",
  help="Output field: ifdd (Freeze Degr. Days) or ihdd (Heat Degr. Days)",
  choices=['ifdd','ihdd'],
  required=True
)
parser.add_argument(
  "--load", 
  help="Load saved ifdd/ihdd tmp file, continue from last record (1), start from time 0 (0 default)", 
  choices=[0,1], 
  default=0, 
  type=int
)
parser.add_argument(
  "--debug", 
  help="1 - run in debug mode, nothing saved, default = 0",
  default=0,
  type=int
)

args = parser.parse_args()

regn       = args.regn
fld_name   = args.field
load_saved = args.load == 1
run_debug  = args.debug == 1

fyaml = 'config_ithkn_predictor.yaml'
with open(fyaml) as ff:
  config_predictor = safe_load(ff)

# Load parameters:
regn_name = config_predictor["regn"][regn]["name"]
lat0      = config_predictor["regn"][regn]["lat_bnd"]
tstep     = config_predictor["params"]["tstep"]
dxy       = config_predictor["params"]["dxy"]
YS        = config_predictor["params"]["ys"]
YS        = config_predictor["params"]["ys"]
YE        = config_predictor["params"]["ye"]
tstep_era = config_predictor["params"]["tstep"]

intgr_time = 90  # Time for freeze degree days accumulation, back from current time
dump_tstp = 50
Tfrz = -1.85    # ocea freezing T

DIRS = {
  "pthithkn" : config_predictor["linregr"]["pthithkn"],
  "pthiconc" : config_predictor["linregr"]["pthiconc"],
  "pthsst"   : config_predictor["linregr"]["pthsst"],
  "pthssh"   : config_predictor["linregr"]["pthssh"],
  "pthui"    : config_predictor["linregr"]["pthui"],
  "pthvi"    : config_predictor["linregr"]["pthvi"],
  "ptht2m"   : config_predictor["linregr"]["ptht2m_1hr"],
  "pthgmapi" : config_predictor["linregr"]["pthgmapi"],
  "pthout"   : config_predictor["linregr"]["pthout"],
  "ithkntmp" : config_predictor["linregr"]["ithkntmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "iconctmp" : config_predictor["linregr"]["iconctmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "ssttmp"   : config_predictor["linregr"]["ssttmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "divutmp"  : config_predictor["linregr"]["divutmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "sattmp"   : config_predictor["linregr"]["sattmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "dfrztmp"  : config_predictor["linregr"]["dfrztmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "heattmp"  : config_predictor["linregr"]["heattmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  "dfrztmp"  : config_predictor["linregr"]["dfrztmp"].format(YS=YS, YE=YE, dxy=dxy, regn=regn),
  }

# Get time steps, grid points:
pthinfo = config_predictor["params"]["pthinfo"]
fltime = config_predictor["params"]["fltime"].format(regn=regn, tstep=tstep)
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


# Add previous days for integrating SAT
# integrating heat deegre days
# Previous (to start) year should exist !
dnmbS = DNMB[0]   # actual start day
#YRS, MMS, DDS = mtime.datevec(dnmbS)[:3]
dnmbP = dnmbS - intgr_time - 1  # previous intgr time preiod, start day
Ypr, Mpr, Dpr = mtime.datevec(dnmbP)[:3]
DNMBprv = micepr.derive_time(Ypr, Ypr, tstep_era)

# Find closest time:
idx0 = max(np.argmin(abs(DNMBprv - dnmbP)) - 2, 0) # add extra index
nrec_prev = len(DNMBprv) - idx0     # how many records to keep for intgr Tfrz
# prepand first days for integrating Tfrz before the start:
DNMB_run = DNMB.copy()
DNMB = np.concatenate((DNMBprv[idx0:], DNMB))



# Load gmapi:
pthindx = DIRS["pthgmapi"]
flout = f"gmapi_ERA5_to_GLORYS_{regn}.nc"
dflout = os.path.join(pthindx, flout)
with xr.open_dataset(dflout) as ds:
  LONE = ds["era_longit"].values
  LATE = ds["era_latit"].values
  IGLR = ds["glorys_indx"].values
  JGLR = ds["glorys_jndx"].values
  IERA = ds["era_indx"].values
  JERA = ds["era_jndx"].values

LONE = (LONE + 360) % 360

def debug_ifdd(ipp, Tsurf, hlat, hlon, SATprv, Tfrz, fld_pnts):
  """
  Debugging Intgr Freeze Dgr Days
  """
  T2m_ij = np.asarray(Tsurf)
  LATI = hlat[JG,IG]
  LONI = hlon[JG,IG]
  xpp = LONI[ipp]
  ypp = LATI[ipp]
  t2m_prv = SATprv[ipp,:]  
  ndays_frz = np.count_nonzero(t2m_prv < Tfrz)
  ihdd = np.nansum(Tfrz - t2m_prv[t2m_prv < Tfrz])
  print(f"  Check pnt: x={xpp:.2f}W, y={ypp:.2f}N, N integr. days: {intgr_time}")
  print(f"  N days T < {Tfrz:.2f}: {ndays_frz}, min/max T: {np.min(t2m_prv):.2f} / {np.max(t2m_prv):.2f}"
        f" IntgrFreeze: {ihdd:.2f}") 
  print(f"  ALL: Min / max T2m: {np.min(T2m_ij):.2f} / {np.max(T2m_ij):.2f}")
  print(f"  ALL: min/max IFDD: {np.min(fld_pnts):.1f} / {np.max(fld_pnts):.1f}")    


def debug_ihdd(ipp, Tsurf, hlat, hlon, SATprv, Tfrz, fld_pnts):
  """
  Debugging Intgr Heat Dgr Days
  """
  T2m_ij = np.asarray(Tsurf)
  LATI = hlat[JG,IG]
  LONI = hlon[JG,IG]
  xpp = LONI[ipp]
  ypp = LATI[ipp]
  t2m_prv = SATprv[ipp,:]  
  ndays_pos = np.count_nonzero(t2m_prv > Tfrz)
  ihdd = np.nansum(t2m_prv[t2m_prv > Tfrz] - Tfrz)
  print(f"  Check pnt: x={xpp:.2f}W, y={ypp:.2f}N, N integr. days: {intgr_time}")
  print(f"  N days T > {Tfrz:.2f}: {ndays_pos}, min/max T: {np.min(t2m_prv):.2f} / {np.max(t2m_prv):.2f}"
        f" IntgrHeat: {ihdd:.2f}") 
  print(f"  ALL: Min / max T2m: {np.min(T2m_ij):.2f} / {np.max(T2m_ij):.2f}")
  print(f"  ALL: min/max IHDD: {np.min(fld_pnts):.1f} / {np.max(fld_pnts):.1f}")    


# Construct predictor sst time series for all locations, 
# Or load previously saved
pthout = DIRS["pthout"]
fltmp = DIRS["dfrztmp"] if fld_name == 'ifdd' else DIRS["heattmp"]
dfltmp = os.path.join(pthout, fltmp)

print(f"\n   Deriving {fld_name} \n")

irec_start = 0
YY = None
if load_saved:
  print(f"Loading saved {dfltmp}, will start from last saved record")
  if not os.path.isfile(dfltmp):
    print(f"Missing tmp file {dfltmp}\n  start from time = 0")
  else:
    # Do not load saved JG, IG from this file: 
    # Keep using those saved in ithkn - should be identical
    data = np.load(dfltmp)
    YY = data["YY"]
    DNMBprv = data["DNMBprv"]
    DNMB_check = data["DNMB"]  

    # Check that this is the right time series:
    dtmp = np.floor(np.abs(DNMB_run - DNMB_check))
    assert np.max(dtmp) == 0, "Check DNMB - dates do not match with saved time series"

    DNMB = np.concatenate((DNMBprv, DNMB_check))
    #Find last saved record, no nans in the column:
    processed = np.all(np.isfinite(YY), axis=0)
    irec_start = np.count_nonzero(processed)
    print(f"Next record to start {irec_start}")


# Find GLORYS - ERA5 pairs:
flg2e = f"glorys2era_pairs_{regn}.npz"
dflg2e = os.path.join(pthinfo, flg2e)
if os.path.isfile(dflg2e):
  print(f"Reading ERA5 indices corresponding GLORYS grid points")
  A = np.load(dflg2e)
  JE = A["JE"]
  IE = A["IE"]

else:
  print("Finding JERA, IERA to match JGLR, IGLR")
  assert len(JG)==len(IG)
  JE = np.empty(len(JG), dtype=int)
  IE = np.empty(len(IG), dtype=int)
  # Alternative to distance-approach, build a lookup table:
  gmapi = {(jg,ig):(je,ie)
           for jg,ig,je,ie in zip(JGLR,IGLR,JERA,IERA)}

  for k, (jj, ii) in enumerate(zip(JG, IG)):
    if k > 0 and k % 500 == 0:
      prc = k/len(JG)*100.
      print(f"  {prc:.2f}% processed")
    key = (int(jj), int(ii))
    if key not in gmapi:
      raise ValueError(f"Missing gmapi entry for GLORYS index {key}")
    JE[k], IE[k] = gmapi[key]


def check_era_glorys_coord(LATE, LONE, IE, JE, hlon, hlat, IG, JG, dmax=27e3):
  """
   Check that gmapi indices are correct
   dmax ~ ERA5 res. (0.25 degree)
  """
  print("Checking gmapi indices ...")
  DD_max = 0
  for k, (jj, ii) in enumerate(zip(JG, IG)):
    jje = JE[k]
    iie = IE[k]
    lon_era = LONE[iie]
    lat_era = LATE[jje]
    lon_glr = hlon[jj,ii]
    lat_glr = hlat[jj,ii]

    DD = dist_sphcrd(lat_era, lon_era, lat_glr, lon_glr)
    if k > 0 and k % 500 == 0:
      prc = k/len(JG)*100.
      print(f"  {prc:.2f}% processed, k={k} dist = {DD:.4f} m")

    DD_max = np.max([DD_max, DD])

    if DD > dmax:
      print(f"ERA point is {DD}m apart from GLORYS, k={k}")
      print(f"ERA lon={lon_era} lat={lat_era}")
      print(f"GLORYS lon={lon_glr} lat={lat_glr}")
      raise Exception("ERR: check not passed")

  print(f"Checked gmapi: OK, overall max dist = {DD_max} m")


check_gmapi = True
if check_gmapi:
  check_era_glorys_coord(LATE, LONE, IE, JE, hlon, hlat, IG, JG)

SATprv = np.zeros((len(IE), nrec_prev)) # should include intgr. period
days_frz = np.zeros(nrec_prev)    # time stamps of saved SAT
"""
  Construct time series of response variable SAT
  2D: locations x time

  First N records will be skipped before actual start date
  To populate SAT array with temp for integrating Heat Dgr Days
"""
iStart = np.where(DNMB == dnmbS)[0][0]
npnts = len(JG)
nrecs = len(DNMB_run)
irec = 0
YRold = 1900
if YY is None:
  YY = np.zeros((npnts, nrecs), dtype=float)*np.nan
for irec0, dnmb0 in enumerate(DNMB):
  YR, MM, DD = mtime.datevec(dnmb0)[:3]

  print(f"Reading ERA5 SAT {YR}/{MM:02d}/{DD:02d}")

  rdate = int(YR*1e4 + MM*100 + DD)

  #flt2m = f"era5_2mTemp_daily7day_Arctic_{YR}.nc"
  #flt2m = f"era5_2mTemp_daily{tstep_era}day_{regn_name}_{YR}.nc"
  flt2m = f"ERA5_reanalysis_sLevels_1hr_0.25x0.25_2m-temperature_{YR}.nc"
  ptht2m = DIRS['ptht2m']
  dflt2m = os.path.join(ptht2m, flt2m)

  T2d = read_era5_T2m(dflt2m, dnmb0)

  Tsurf = []
  for jje, iie in zip(JE, IE):
    Tsurf.append(T2d[jje,iie])

  # Update SAT with previous SATs:
  # add current day at the end and delete 1st column 
  SATprv[:, :-1] = SATprv[:, 1:]
  SATprv[:,-1] = Tsurf
  days_frz[:-1] = days_frz[1:]
  days_frz[-1] = dnmb0

  # Cycle over N previous days 
  # Until the start date
  if dnmb0 < dnmbS:
    continue

  irec = irec0 - iStart
  assert irec >= 0, f"ERROR irec0={irec0} should be > {iStart}"
  assert irec < YY.shape[1]
  if irec < irec_start:
    continue
  
  fld_pnts = []
  assert np.all(days_frz > 0), "days_frz not populated, there are 0s"
  assert np.all(np.diff(days_frz)>0), "days_frz not increasing" 
  
  if fld_name == 'ifdd':
    fld_pnts = micepr.intgr_Tfrz(SATprv, intgr_time, Tfrz, days_frz)
  else:
    fld_pnts = micepr.intgr_HeatDgr(SATprv, intgr_time, Tfrz, days_frz)

  YY[:,irec] = np.asarray(fld_pnts)

  if run_debug:
    ipp = 2085 if regn == 'north' else 6726
    if fld_name == 'ifdd':
      debug_ifdd(ipp, Tsurf, hlat, hlon, SATprv, Tfrz, fld_pnts)
    else:
      debug_ihdd(ipp, Tsurf, hlat, hlon, SATprv, Tfrz, fld_pnts)

  if (irec + 1) % dump_tstp == 0 and not run_debug: 
    print(f"TMP step: Saving {fld_name} time series  --> {dfltmp}")
    np.savez(dfltmp,
           YY=YY,
           JG=JG,
           IG=IG,
           JE=JE,
           IE=IE,
           DNMBprv=DNMB[:iStart],
           DNMB=DNMB[iStart:])

if irec_start < len(DNMB) or run_debug:
  # No need to save if already everything processed
  print(f"END TMP step: Saving {fld_name} time series and IG, JG --> {dfltmp}")
  np.savez(dfltmp,
         YY=YY,
         JG=JG,
         IG=IG,
         JE=JE,
         IE=IE,
         DNMBprv=DNMB[:iStart],
         DNMB=DNMB[iStart:])



f_check = False
if f_check:
  DNMB_saved = DNMB[iStart:]

  plt.ion()

  # Plot ERA5 SAT at the sample locations:
  ir0 = 1920    # Jan 1 2025
  dnmbP = DNMB_saved[ir0]
  YRp, MMp, DDp = mtime.datevec(dnmbP)[:3]

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
  ax1.contour(LMsk, [0.99], linestyles='solid', colors=[(0.5,0.8,1)])
  #ax1.scatter(IG, JG, s=5, color=(0.8,0.4,0), marker='.')

  AA = YY[:,ir0]
  sc = ax1.scatter(
    IG, JG,
    c=AA,
    cmap='jet',
    s=20,          # marker size
    vmin=0,        # optional color scale limits
    vmax=1500
  )
  
  ax1.set_title(f"HeatDgr Days, {YRp}/{MMp:02d}/{DDp:02d}")

  plt.colorbar(sc, ax=ax1, label='HeatDgrDays')


  # Check time series for 1 point:
  xp = 120
  yp = -58
  LATI = LAT[JG,IG]
  LONI = LON[JG,IG]
  #ipp = np.argmin((LATI-yp)**2 + (LONI-xp)**2)
  ipp = np.argmin((LATI-yp)**2)

  tser = YY[ipp,:]
  ax1.cla()
  ax1.plot(tser)
  ax1.set_xlim(1860, len(tser))



