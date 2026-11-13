import numpy as np
import mod_time as mtime
import xarray as xr
import os
import glob
from yaml import safe_load
from scipy.interpolate import interp1d
#import sys
#import matplotlib.pyplot as plt
#import importlib
from mod_misc1 import dist_sphcrd


def full_model_name(ml_model):
  """
    Linear regression models:
      OLS 1: trained over 1993-2002, no predictor for interannual trend
      OLS 2: trained 1993-2025, added predictor for interannual trend
    Random Forest (all 1993-2025):
      RF1:  Ntrees 100
      RF2:  Ntrees 200
      RF3:  Ntrees 5
      RF4:  Ntrees 2
   Added predictor: integrated Heat Degree Days (IHDD)
      RF11: Ntrees 50
      RF12: Ntrees 100
      RF13: Ntrees 10
    Hist. Grad Boost Regressor (decision tree)
      GBR1: max_leaf = 8,  max_iter=500,  max_depth=8
      GBR2: max_leaf = 31, max_iter=1000, max_depth=10
      GBR3: max_leaf = 63, max_iter=1500, max_depth=15
   Added predictor: integrated Heat Degree Days (IHDD)
      GBR11: max_leaf = 63, max_iter=1500, max_depth=15 
  """
  if ml_model[:4] == 'clim':
    mname = "clim"
  elif ml_model[:3] == 'gbr':
    nmb = int(ml_model[3:])
    mname = f"GradBoostRegr_model{nmb:02d}"
  elif ml_model[:2] == 'rf':
    nmb = int(ml_model[2:])
    mname = f"RF_model{nmb:02d}"
  else:
    print(f"Unrecognized ML model: {ml_model}")
    return None

  return mname

def find_oras_ice_file(dpath, rdate):
  """
  Find ORAS ocean / ice file  file containing rdate.
  Ocean field: sst only

  Filename format:
    0001_1d_YYYYMMDD_YYYYMMDD_icemod.nc
  """
  for fname in glob.glob(os.path.join(dpath, "*_1d_*_icemod.nc")):
    parts = os.path.basename(fname).split("_")
    date_start, date_end = int(parts[2]), int(parts[3])

    if date_start <= rdate <= date_end:
      return fname

  return None

def check_oras_mesh025_coord(LATE, LONE, IE, JE, hlon, hlat, IG, JG, dmax=28e3, vocal=False):
  """
  For 2 sets of (J,I) indices on 2 mesh grids
  Check that they correspond to the grid points
  That are within dmax distance 

  Check that gmapi indices are correct
  dmax ~ ORAS res. (0.25 degree ~= 27 km)
  """
  print("Checking gmapi indices ...")
  DD_max = 0
  lat_min_oras = np.max(LATE[0,:])
  for k, (jj, ii) in enumerate(zip(JG, IG)):
    jje = JE[k]
    iie = IE[k]
    lon_oras = LONE[jje, iie]
    lat_oras = LATE[jje, iie]
    lon_025 = hlon[jj,ii]
    lat_025 = hlat[jj,ii]

    # Skip grid points south of the ORAS5 southern boundary:
    if lat_025 < lat_min_oras:
      continue

    DD = dist_sphcrd(lat_oras, lon_oras, lat_025, lon_025)
    if vocal and (k > 0 and k % 20000) == 0:
      prc = k/len(JG)*100.
      print(f"  {prc:.2f}% processed, k={k} dist = {DD:.4f} m")

    DD_max = np.max([DD_max, DD])

    if DD > dmax:
      print(f"ORAS5 point is {DD}m apart from mesh025, k={k}")
      print(f"ORAS5 lon={lon_oras} lat={lat_oras}, jje={jje}, iie={iie}")
      print(f"mesh025 lon={lon_025} lat={lat_025}, jj={jj}, ii={ii}")
      raise Exception(f"ERR: check not passed, DD={DD} > {dmax}")

  print(f"Checked gmapi: OK, overall max dist = {DD_max} m")

def read_daily_era5(dfl, dnmb0):
  """
  Derive daily fields from ERA5 hourly
  """
  yr, mm, dd = mtime.datevec(dnmb0)[:3]
  rdate0 = yr*10000 + mm*100 + dd

  with xr.open_dataset(dfl) as ds:
    # Check that all 24 hourly records are present
    if ds.sizes["Time"] != 24:
      raise ValueError(f"Expected 24 hourly records, found {ds.sizes['Time']}")

    # Check that all records belong to the requested day rday0=YYYYMMDD
    dates = ds.time.dt.strftime("%Y%m%d").astype(int).values

    if not np.all(dates == rdate0):
      raise ValueError(
          f"Time mismatch: requested {rdate0}, "
          f"file contains {np.unique(dates)}"
      )

    # Daily averages
    T2m_mean = ds["temperature"].mean(dim="Time", skipna=True).values

    # Convert K --> C if needed:
    if np.max(T2m_mean) > 270 and np.min(T2m_mean) > 150:
      T2m_mean -= 273.15

  return T2m_mean

def check_era_mesh025_coord(LATE, LONE, IE, JE, hlon, hlat, IG, JG, dmax=27e3, vocal=False):
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
    if vocal and (k > 0 and k % 20000) == 0:
      prc = k/len(JG)*100.
      print(f"  {prc:.2f}% processed, k={k} dist = {DD:.4f} m")

    DD_max = np.max([DD_max, DD])

    if DD > dmax:
      print(f"ERA point is {DD}m apart from mesh025, k={k}")
      print(f"ERA lon={lon_era} lat={lat_era}")
      print(f"mesh025 lon={lon_glr} lat={lat_glr}")
      raise Exception("ERR: check not passed")

  print(f"Checked gmapi: OK, overall max dist = {DD_max} m")

def construct_ydays(DNMB, npnts, order_fast="time"):
  """
    Predictor year day represent as
    cos(2*pi*365/jday) and sin(2*pi*365/jday) 
    to avoid end year  / start year discontinuity
    when jday = 365 --> jday 1

    Note that alfa1*cos(theta) + alfa2*sin(theta) = B*sin(theta + phase)
    where B = sqrt(alfa1**2 + alfa2**2), i.e.
    correctly represents the amplitude and phase

    Construct 1D array (time x npnts)
    order_fast - dimension that changes fast
    order_fast: time = loc1: 1, 2, 3, ..., nrecs, loc2: 1,2,3, ...nrecs, ... 
           coord = time1: 1, 2, 3, ..., npnts, time2: 1, 2, 3, ..., npnts

  """
  nrecs = len(DNMB)
  cosD = np.empty(nrecs)
  sinD = np.empty(nrecs)

  DV = mtime.datevec2D(DNMB)
  years = DV[:,0]
  years_range = np.arange(years[0], years[-1]+1)

  for YR in years_range:
    mask  = years == YR
    ndays_year = mtime.year_days(YR)
    dnmb_subset = DNMB[mask]

    assert dnmb_subset.size > 0, "year {YR} has no records"

    # Compute day of the year wrt to Jan 1:
    dnmbJ1 = mtime.datenum([YR,1,1])
    jdays = dnmb_subset - dnmbJ1
    c2rad = 2 * np.pi / ndays_year
    cos_days = np.cos(c2rad * jdays)
    sin_days = np.sin(c2rad * jdays)

    cosD[mask] = cos_days
    sinD[mask] = sin_days

  assert len(cosD) == len(sinD) == nrecs, \
    "cos days or sin days numerb of records is incorrect"

  if order_fast == "time":
    cosD = np.tile(cosD, npnts)
    sinD = np.tile(sinD, npnts)

  elif order_fast == "coord":
    cosD = np.repeat(cosD, npnts)
    sinD = np.repeat(sinD, npnts)

  else:
    raise ValueError(f"Unknown order: {order}")

  return cosD, sinD

def construct_coord_sphere(hlon, hlat, IG, JG, nrecs, order_fast="time"):
  """
  Represent geographic coordinates as Cartesian coordinates
  on the unit sphere to avoid longitude discontinuities.
  """
  lon = np.deg2rad(hlon[JG, IG])
  lat = np.deg2rad(hlat[JG, IG])

  Xcrd = np.cos(lat) * np.cos(lon)
  Ycrd = np.cos(lat) * np.sin(lon)
  Zcrd = np.sin(lat)

  if order_fast == "coord":
    Xcrd = np.tile(Xcrd, nrecs)
    Ycrd = np.tile(Ycrd, nrecs)
    Zcrd = np.tile(Zcrd, nrecs)

  elif order_fast == "time":
    Xcrd = np.repeat(Xcrd, nrecs)
    Ycrd = np.repeat(Ycrd, nrecs)
    Zcrd = np.repeat(Zcrd, nrecs)

  else:
    raise ValueError(f"Unknown order: {order}")

  return Xcrd, Ycrd, Zcrd

def construct_mean_ithkn(npnts, DNMB, dflmni, order_fast="time", Mavrg=3, fld='ithkm'):
  """
    Construct a predictor from monthly mean ice thickness (or ice volume)
    for the specified dates.

    For each date, the predictor is the average over the previous Mavrg
    months, excluding the current month.

    Example:
        June, Mavrg=3  --> average(Mar, Apr, May)

    For the beginning of the record (1993), where fewer than Mavrg previous
    months are available, average over all available previous months.

    npnts = number of grid points used for creating stat model
    DNMB  = 1D array of date numbers for stat. model
    dflmni = path/file.npz with numpy arrays of monthly mean spatially averaged  ice statistics
  """
  # Read saved monthly ice volume and ice mean thkn:
  print(f"Loading ice vol and mean ice thickness --> {dflmni}")
  data = np.load(dflmni)
  if fld == 'ithkm':
    ITHK = data['ITHKM']
  elif fld == 'ivol':
    ITHK = data['IVOL']
  else:
    raise ValueError(
      f"Unsupported field '{fld}'. Expected 'ithkm' or 'ivol'."
    )

  DNMB_mo = data['DNMB']  # monthly date stamps
  DVM = mtime.datevec2D(DNMB_mo)[:,:3]
  ITHKN_comb = np.zeros((len(DNMB)))*np.nan   # array of combined mean ice thickness, averaged over Mavrg
  for irec, dnmb0 in enumerate(DNMB):
    YR, MM, DD = mtime.datevec(dnmb0)[:3]
    indx = np.where( (DVM[:,0] == YR) & (DVM[:,1] == MM))[0]
    assert len(indx) == 1, (
      f"Expected one monthly record for {YR}/{MM:02d}, "
      f"found {len(indx)}."
    )

    # Previous-month averaging window
    iend = indx[0] - 1
    if iend < 0:
      # No previous months available (first month in record)
      mith = ITHK[0]
    else:
      istart = max(0, iend - Mavrg + 1)
      mith = np.mean(ITHK[istart:iend + 1])

    print(f"Mean ice thickness {YR}/{MM}/{DD}: {mith} m")
    ITHKN_comb[irec] = mith

  assert np.all(np.isfinite(ITHKN_comb)), (
      "ITHKN_comb contains NaN or infinite values."
  )

  # Repeat for all locations
  if order_fast == "time":
    ITHKN_crd = np.tile(ITHKN_comb, npnts)
  elif order_fast == "coord":
    ITHKN_crd = np.repeat(ITHKN_comb, npnts)
  else:
    raise ValueError(f"Unknown order: {order}")

  return ITHKN_crd

def calc_divU(dltI, dltJ, ii, jj, U2d, V2d, Acell, DX, DY):
  """
    Compute average div u ice over specified region
    Assuming output fieds are at the cell centers
  """
  divU = 0
  jdm, idm = U2d.shape

  Intgr = 0.0
  # Define the box around the grid point:
  iS = int(ii - dltI)
  iE = int(ii + dltI)
  iS = np.max([iS, 0])
  iE = np.min([iE, idm-1])

  jS = int(jj - dltJ)
  jE = int(jj + dltJ)
  # Boundaries, better - for global grid
  # use grid points at the opposite side
  jS = np.max([0, jS])
  jE = np.min([jE, jdm-1])

  # Integrate along the boundary:
  # Note different sign of outward norm vectors 
  # along the box sides
  Intgr = (
    np.sum(U2d[jS:jE+1, iE] * DY[jS:jE+1, iE])
    - np.sum(U2d[jS:jE+1, iS] * DY[jS:jE+1, iS])
    + np.sum(V2d[jE, iS:iE+1] * DX[jE, iS:iE+1])
    - np.sum(V2d[jS, iS:iE+1] * DX[jS, iS:iE+1])
  )

  # subtract half of the four corner contributions
  Intgr -= 0.5 * (
    U2d[jS, iE] * DY[jS, iE]
    + U2d[jE, iE] * DY[jE, iE]
    - U2d[jS, iS] * DY[jS, iS]
    - U2d[jE, iS] * DY[jE, iS]
    + V2d[jE, iS] * DX[jE, iS]
    + V2d[jE, iE] * DX[jE, iE]
    - V2d[jS, iS] * DX[jS, iS]
    - V2d[jS, iE] * DX[jS, iE]
  )

  # Space-Average divergence:
  Area = np.sum(Acell[jS:jE+1, iS:iE+1])
  assert Area > 0, f"ii={ii}, jj={jj}, dltI={dltI}, dltJ={dltJ}, Area = {Area}"
  div_uice = Intgr / Area

  return div_uice

def subset_cice6_divu(dflui, IG, JG, hlon, hlat, dxy):
  """
    Subset and compute spatial-average divU ice for specified grid points 
    CICE6 analysis
    dflui:  CICE6 dir/file name with ice velocity fields
    IG, JG:  grid points where divU is computed
    dxy:     distance between the subsampled grid points used 
             for deriving linear regr. parameters (km)

    Note no time averaging here, 
    see model training scripts if averaging is needed
  """
  from mod_mom6 import dx_dy

  DX, DY = dx_dy(hlon, hlat)
  Acell = DX*DY

  with xr.open_dataset(dflui) as dsice:
    A2d = dsice['uvel'].values

  assert A2d.shape == hlon.shape, f"uvel CICE6 shape does not conform {A2d.shape}"

  # Treat nans as no ice grid cells
  U2d = np.nan_to_num(A2d, nan=0.0)

  with xr.open_dataset(dflui) as dsice:
    A2d = dsice['vvel'].values

  # Treat nans as no ice grid cells
  V2d = np.nan_to_num(A2d, nan=0.0)

  fld_pnts = []
  ncheck = 90000
  print(f"Calculating ice divergence averaged over +/- dxy={dxy}km, {len(JG)} pnts ...")
  for icnt, (jj, ii) in enumerate(zip(JG, IG)):
    if icnt > 0 and icnt % ncheck == 0:
      print(f"    {icnt/len(JG)*100.:.2f}% done ...")
    # Estimate box size based on min distance criterion
    dltX = DX[jj,ii]*1e-3  # km
    dltY = DY[jj,ii]*1e-3  # km
    dltI = int(np.ceil(dxy / dltX))
    dltJ = int(np.ceil(dxy / dltY))
    div_uice = calc_divU(dltI, dltJ, ii, jj, U2d, V2d, Acell, DX, DY)
    fld_pnts.append(div_uice)

  return fld_pnts


def intgr_Tfrz(SATprv, intgr_time, Tfrz, days_frz):
  """
    Use Zubov definition of accumulated freeze degree days
    sum(Tfrz-T), when SAT < Tfrz
    Linear interpolation is safer
  """
  Tintrp = np.arange(days_frz[-1] - intgr_time, days_frz[-1]+1)
  #cs = CubicSpline(days_frz, SATprv, axis=1)
  # SATs during the requested previous Ndays
  #SATi = cs(Tintrp)

  interp = interp1d(days_frz, SATprv,
                  axis=1,
                  kind='linear')
  SATi = interp(Tintrp)

  T2frz = np.where(SATi < Tfrz, SATi, np.nan)
  intgrFDD = np.nansum(Tfrz - T2frz, axis=1)  # integrated Freeze degree days

  f_chck = False
  if f_chck:
    iC=15
    Ti = SATi[iC,:]
    T0 = SATprv[iC,:]
    ax1.cla()
    ax1.plot(Tintrp, Ti)
    ax1.plot(days_frz, T0,'.-')

  return intgrFDD

def subset_oras5_sst(dflsst, IORsub, JORsub, rdate):
  """
  Subset SST fields for specified grid points
  from ORAS5 fields
  It may keep > 1 days

  dflsst:     SST dir + file name
  IORsub, JORsub : mapping indices from ORAS5 --> mesh025 at analysis grid points
  """
  with xr.open_dataset(dflsst) as dsst:
    # Read time in YYYYMMDD format
    Time = dsst["time_counter"].dt.strftime("%Y%m%d").astype(int).values
    if not np.any(Time == rdate):
      raise ValueError(
          f"Time mismatch: requested {rdate}, "
          f"file contains {Time}"
      )

    itime = np.where(Time == rdate)[0][0]
    SST_oras = dsst["isstempe"].isel(time_counter=itime).values

  # Subample ORAS grid points to the analysis grid points on mesh025
  SST = SST_oras[JORsub, IORsub]

  return SST

def subset_oras5_iconc(dflice, IORsub, JORsub, rdate, varnm='ileadfra'):
  """
  Subset ice conc fields for specified grid points
  from ORAS5 fields
  It may keep > 1 days

  dflsst:     SST dir + file name
  IORsub, JORsub : mapping indices from ORAS5 --> mesh025 at analysis grid points
  """
  with xr.open_dataset(dflice) as dsice:
    # Read time in YYYYMMDD format
    Time = dsice["time_counter"].dt.strftime("%Y%m%d").astype(int).values
    if not np.any(Time == rdate):
      raise ValueError(
          f"Time mismatch: requested {rdate}, "
          f"file contains {Time}"
      )

    itime = np.where(Time == rdate)[0][0]
    A2d = dsice[varnm].isel(time_counter=itime).values

  # Treat nans as no ice grid cells
  A2d = np.nan_to_num(A2d, nan=0.0)

  Iconc = A2d[JORsub, IORsub]
  Iconc[Iconc > 1] = 1.
  Iconc[Iconc < 0] = 0.

  return Iconc

def subset_gdas_ifdd_ihdd(dfifdd, IG, JG, idim, jdim):
  """
  Extract IFDD or IHDD at selected mesh025 grid points.

  dfifdd - dir/file_name.npz of IFDD or IHDD fields
  IFDD and IHDD fields are derived from daily-averaged GDAS T2m fields
  and interpolated onto the mesh025 grid.

  See: derive_intgrFDD_mesh025
  """
  # Find out if requested IFDD or IHDD
  ifdd_case = 'ifdd' in dfifdd.lower()
  varname = "IFDD" if ifdd_case else "IHDD"

  DATA = np.load(dfifdd)
  IDD = DATA[varname]
  JP  = DATA["JG"]
  IP  = DATA["IG"]

  jdim_saved = DATA["jdim"]
  idim_saved = DATA["idim"]

  assert idim_saved == idim and jdim_saved == jdim, (
      f"{varnm} grid {jdim_saved} x {idim_saved} "
      f"does not match jdim={jdim}, idim={idim}"
  )

  # Check dims:
  assert np.all((JP >= 0) & (JP < jdim)), "Saved JP indices out of bounds"
  assert np.all((IP >= 0) & (IP < idim)), "Saved IP indices out of bounds"

  # Insert saved IFDD / IHDD values back onto the full mesh025 grid
  A2d = np.full((jdim, idim), np.nan)
  A2d[JP, IP] = IDD

  # Extract the requested ML grid points
  IDD = A2d[JG, IG]

  assert np.all(np.isfinite(IDD)), (
    "Some requested points are outside the saved {varnm} domain or contain NaN"
  )

  return IDD

def subset_gdas_t2m(dflt2m, IG, JG, idim, jdim):
  """
    Subset daily T2m GDAS onto ML grid points JG, IG
    daily T2m GDAS interpolated onto mesh025 in interp_GDAS_T2m_daily_mesh025_Ndays.py
  """
  with xr.open_dataset(dflt2m) as dt2m:
    A2d = dt2m["temp_2m"].isel(time=0).values

  assert A2d.shape == (jdim, idim), f"Expected shape does not conform T2m saved: {A2d.shape}"

  T2m = A2d[JG,IG]

  return T2m 

def construct_bottom_topo(
    HH, IG, JG, nrecs, order_fast="time", Hmin=-5000
  ):
  """
  Derive bottom depth at training points (IG,JG)

  Assumed:
    - Ocean depth is negative downward in the source file.
    - Land is >= 0

  order_fast="coord": coordinates vary fastest in the record sequence.
  order_fast="time":  time varies fastest in the record sequence.

  Depth are normalized to vary -1 <= h <= 1
  using min depth Hmin, which should match Hmin
  in the training setup 
  Depth < Hmin, truncated to Hmin (to have -1 as the lower bound)
  """
  if order_fast not in ("coord", "time"):
    raise ValueError("order_fast must be 'coord' or 'time'.")

  if Hmin is not None and abs(Hmin) < 1.e-10:
    Hmin = None

  normlz = Hmin is not None

  assert np.nanmin(HH) < 0, f"Bottom depths: Ocean points mush have negative depths"
  LMocn = np.isfinite(HH) & (HH < 0)

  topo_pnts = HH[JG, IG]

  if normlz:
    topo_norm = np.maximum(topo_pnts / abs(Hmin), -1.0)
  else:
    topo_norm = topo_pnts

  # Expand over records according to record ordering.
  if order_fast == "coord":
    Btopo = np.tile(topo_norm, nrecs)
  else:
    Btopo = np.repeat(topo_norm, nrecs)

  return Btopo

def construct_predictors_day(
    dnmb0, fyaml, MLPARAM, PRED_NAMES, regn,
    order_fast="time", Mavrg=3
  ):
  """
  Construct predictors for 1 day forecast
  PRED_NAMES - predictors used in this model
  predictor array must have exactly the same columns, in exactly the same order,
  as when ML was trained.
  
  hlon, hlat - grid where prediction is done (mesh025)
  IG, JG - grid point indices where prediction is done
  MLYAML - YAML dictionary with input / output directories, etc.
  PREAD_NAMES - list of predictors in the right order
  sqrt_frzdays - True: use sqrt of integrated freezing degree days
  sst_max - max sst threshold for possible sea ice
  regn - region of prediction
  intgr_frz - for freeze degree days, integration period, days
  dxy - ice length scale, used for calc. ice predictors and gird point subset
  """

  hlon      = MLPARAM["LON"]
  hlat      = MLPARAM["LAT"]
  IG        = MLPARAM["IG"]    # analysis grid points on mesh025
  JG        = MLPARAM["JG"]
  IORsub    = MLPARAM["IORsub"]  # mapping indx ORAS5 --> mesh025
  JORsub    = MLPARAM["JORsub"]
  IERsub    = MLPARAM["IERsub"]  # mapping ERA5 --> mesh025
  JERsub    = MLPARAM["JERsub"]  # mapping ERA5 --> mesh025
  HH        = MLPARAM["topo"]
  sst_max   = MLPARAM["sst_max"]
  intgr_frz = MLPARAM["intgr_frz"]
  Tfrz      = MLPARAM["Tfrz"]  
  sqrt_frzdays = MLPARAM["sqrt_frzdays"]
  dxy       = MLPARAM["dxy"]
  Hbtm_min  = MLPARAM["Hbtm_min"]

  with open(fyaml) as ff:
    MLYAML = safe_load(ff)


  DNMB = np.asarray([dnmb0])
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  rdate = YR0*10000 + MM0*100 + DD0
  Iconc = None
  SST = None 
  jdim, idim = hlon.shape 
    
  raw = {} # predcitors with not stand. values, can be in any order
  # Coord --> polar coord:
  if 'Xcrd' in PRED_NAMES or 'Ycrd' in PRED_NAMES or 'Zcrd' in PRED_NAMES:
    Xcrd, Ycrd, Zcrd = construct_coord_sphere(
         hlon, hlat, IG, JG, len(DNMB), order_fast=order_fast
         )
      
    raw['Xcrd'] = Xcrd
    raw['Ycrd'] = Ycrd
    raw['Zcrd'] = Zcrd
  
  # Time:
  if 'cosD' in PRED_NAMES or 'sinD' in PRED_NAMES:
    cosD, sinD = construct_ydays(DNMB, len(IG), order_fast=order_fast)
  
    raw['cosD'] = cosD
    raw['sinD'] = sinD

  if 'btmtopo' in PRED_NAMES:
    Btopo = construct_bottom_topo(
      HH, IG, JG, len(DNMB), Hmin=Hbtm_min, order_fast=order_fast
    )
    print(f"btmtopo min/max: {np.min(Btopo):.3f} / {np.max(Btopo):.3f}")
    raw['btmtopo'] = Btopo

  # Mean ice thickness over ice area during previous N months:
  # The predictor accounts for the long-term trend in ice thicknesses
  # Check Mavrg - should match train_linregr_ithkn.py
  # Interannual trend: mean ice thickness previous N months
  # Use GLORYS mean ice thicknesses (prepared on PPAN), SOCA ithkn is unrealistically thin
  # which results in very biased ice thickness prediction
  if 'mnithkn' in PRED_NAMES:
    pthglr = MLYAML["PRED"]["pthglr"]
    if YR0 < 2026:
      fltmp = f"GLORYS_monthly_icevol_ithknmn_{regn}_1993_2025.npz"
    elif YR0 == 2026: 
      fltmp = f"GLORYS_monthly_icevol_ithknmn_{regn}_{YR0}_{YR0}.npz"

    dflmni = os.path.join(pthglr, fltmp)
    assert os.path.isfile(dflmni), f"File is missing: {dflmni}"
    mnithkn = construct_mean_ithkn(len(IG), DNMB, dflmni, order_fast=order_fast, Mavrg=Mavrg)
      
    print(f"mnithkn min/max: {np.min(mnithkn):.3f} / {np.max(mnithkn):.3f}")
    raw['mnithkn'] = mnithkn
      
  # SST
  SST = None
  if 'sst' in PRED_NAMES:
    # ORAS5
    #YRm, MMm, DDm, HRm = mtime.datevec(dnmb0)[:4]
    #rdate = int(YRm)*10000 + int(MMm)*100 + int(DDm)

    print("\nDeriving ORAS5 sst")
    pthsst = MLYAML["PRED"]["pthsst"]
    dflsst = find_oras_ice_file(pthsst, rdate)

    if not os.path.isfile(dflsst):
        raise FileNotFoundError(f"MOM6 analysis SST Not found {dflsst}")

    SST = subset_oras5_sst(dflsst, IORsub, JORsub, rdate)

    assert len(SST) == len(JG), f"Subset of SST does not match expected N of grid points {len(JG)}"

    print(f"SST min/max: {np.min(SST):.3f} / {np.max(SST):.3f}")
    raw['sst'] = SST

  # Ice conc
  if 'iconc' in PRED_NAMES:
    print("\nDeriving ORAS5 iconc")
    varnm   = 'aicen'
    pthice  = MLYAML["PRED"]["pthiconc"]
    dflice  = find_oras_ice_file(pthsst, rdate)

    if dflice is None or not os.path.isfile(dflice):
      raise FileNotFoundError(f"Not found {dflice}")

    Iconc = subset_oras5_iconc(dflice, IORsub, JORsub, rdate)

    # Eliminate ice in the warm ocean:
    if SST is not None and Iconc.shape == SST.shape:
      Iconc[SST > sst_max] = 0.

    raw['iconc'] = Iconc

  # divU ice
  if 'divu' in PRED_NAMES:
    # Computed in derive_divUice_oras5_mesh025.py
    print("\nLoading ORAS5 divu ice")
    pthprd   = MLYAML["PRED"]["pthprd"]
    flnm0    = MLYAML["PRED"]["fldivui"].format(dxy=dxy, rdate=rdate, regn=regn)

    fldivui = f"{flnm0}.npz"
    dfliceout = os.path.join(pthprd, fldivui)
   
    if dfliceout is None or not os.path.isfile(dfliceout):
      raise FileNotFoundError(f"Not found divU {dfliceout}")

    DATA = np.load(dfliceout)
    divU = DATA["DIVUI"]
    I025 = DATA["I025"]
    J025 = DATA["J025"]

    assert divU.shape == IG.shape, f"len(divU)={len(divU)}, expected len(IG)={len(IG)}"
    assert np.max(abs(IG-I025) + abs(JG-J025)) == 0, f"Check analysis points for divU"    
 
    print(f"divU min/max: {np.min(divU):.2e} / {np.max(divU):.2e}")
    raw['divu'] = divU

  # Freeze days and heat days
  need_frzdays  = 'frzdays'  in PRED_NAMES and 'frzdays'  not in raw
  need_heatdays = 'heatdays' in PRED_NAMES and 'heatdays' not in raw

  if need_frzdays:
    print("\nReading ERA5 IFDD ")
    pthifdd = MLYAML["PRED"]["pthifdd"]
    flifdd  = MLYAML["PRED"]["flifdd"].format(Ndays=intgr_frz, dxy=dxy, rdate=rdate, regn=regn)
    flnpz = f"{flifdd}.npz"
    dfifdd = os.path.join(pthifdd, flnpz)
  
    if not os.path.isfile(dfifdd):
        raise FileNotFoundError(f"Not found {dfifdd}")
   
    DATA = np.load(dfifdd)
    frzdays = DATA["IFDD"]
    I025 = DATA["IG"]        # analysis points on mesh025
    J025 = DATA["JG"]        # analysis points on mesh025

    assert frzdays.shape == IG.shape, f"IFDD shape does not agree with {IG.shape}"
    assert np.max(abs(IG-I025) + abs(JG-J025)) == 0, f"Check analysis points for IFDD" 
 
    if sqrt_frzdays:
      frzdays = np.sqrt(frzdays)

    print(f"IFDD min/max: {np.min(frzdays):.3f} / {np.max(frzdays):.3f}")
    raw['frzdays'] = frzdays

  if need_heatdays:
    print("\nReading ERA5 IHDD ")

    pthihdd = MLYAML["PRED"]["pthihdd"]
    flihdd  = MLYAML["PRED"]["flihdd"].format(Ndays=intgr_frz, dxy=dxy, rdate=rdate, regn=regn)
    flnpz = f"{flihdd}.npz"
    dfihdd = os.path.join(pthihdd, flnpz)

    if not os.path.isfile(dfihdd):
        raise FileNotFoundError(f"Not found {dfihdd}")

    DATA = np.load(dfihdd)
    heatdays = DATA["IHDD"]
    I025 = DATA["IG"]
    J025 = DATA["JG"]

    assert heatdays.shape == IG.shape, f"IHDD shape does not agree with {IG.shape}" 
    assert np.max(abs(IG-I025) + abs(JG-J025)) == 0, f"Check analysis points for IHDD" 

    print(f"IHDD min/max: {np.min(heatdays):.3f} / {np.max(heatdays):.3f}")
    raw['heatdays'] = heatdays

  # SAT
  if 'sat' in PRED_NAMES:
    print("\nDeriving ERA5 T2m")

    pthera = MLYAML["ERA5"]["pthera"].format(YR=YR0)
    flera = f"ERA5_forcing_{YR0}-{MM0:02d}-{DD0:02d}.nc"
    dflera = os.path.join(pthera, flera)
    if not os.path.isfile(dflera):
      raise FileNotFoundError(f"Not found: {dflera}")

    T2m = read_daily_era5(dflera, dnmb0)
    # Subet to mesh025 analysis grid points:
    SAT = T2m[JERsub, IERsub]

    print(f"SAT min/max: {np.min(SAT):.3f} / {np.max(SAT):.3f}")
    raw['sat'] = SAT

  # Transform and Standardize predictors
  # Use only active predictors
  # Dict. with standardized arrays in the PRED_NAMES order:
  std_arr = {}
  # Construct predictors dictionary where each
  # predictor is linked to the raw array
  # The order of the predictor should match the order in PRED_NAMES
  # Will create dict like this:
  #   'cosD'    : cosD,  # not standardized
  #   'sinD'    : sinD,  # not standardized
  #    ...
  #    'Zcrd'    : raw['Zcrd']  # not standardized
  #     ...
  #
  pred_dict = {}
  for name in PRED_NAMES:
    pred_dict[name] = raw[name]

  # Combine all standardized or not standardized predictors into a list:
  PRED_FINAL = []
  PRED_FINAL = [pred_dict[name] for name in PRED_NAMES]

  return PRED_FINAL, Iconc, SST


