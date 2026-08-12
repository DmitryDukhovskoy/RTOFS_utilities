import numpy as np
import mod_time as mtime
import xarray as xr
import os
from scipy.interpolate import interp1d
#import sys
#import matplotlib.pyplot as plt
#import importlib

import mod_glorys as mglr

def models_info():
  MODEL_NAMES = {
    0  : "clim",
    1  : "OLS_model1",
    2  : "OLS_model1_1993_2025",
    3  : "RF_model01_north",
    4  : "RF_model02_north",
    5  : "RF_model03_north"
    }

  return MODEL_NAMES

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

def check_dnmb_array(DNMB, print_months=True):
  """
    Check Number of records per year
  """
  DVM = mtime.datevec2D(DNMB)[:,:3]
  YR = DVM[:,0]
  MM = DVM[:,1]

  for year in range(YR[0],YR[-1]+1):
    #nyr = len(np.where(YR == year)[0])
    nyr = np.count_nonzero(YR == year)
    print(f"YEAR={year}, N records={nyr}")
    if print_months:
      for month in range(1,13):
        #nmo = len(np.where( (YR == year) & (MM == month) )[0])
        nmo = np.count_nonzero((YR == year) & (MM == month))
        print(f"      MM={month:02d}, N records={nmo}")

  print(f"Total N records = {len(DNMB)}")

  return


def subset_glorys_iconc(dflice, IG, JG):
  """
    Subset ice conc fields for specified grid points 
    GLORYS daily 
    dflice - GLORYS dir/iconc_file_name.nc
  """
  with xr.open_dataset(dflice) as dsice:
    A2d = dsice['siconc'].isel(time=0).values.squeeze()

  # Treat nans as no ice grid cells
  A2d = np.nan_to_num(A2d, nan=0.0)

  Iconc = A2d[JG,IG]
  Iconc[Iconc > 1] = 1.
  Iconc[Iconc < 0] = 0.

  return Iconc

def subset_glorys_sst(dflice, IG, JG):
  """
    Subset SST fields for specified grid points 
    GLORYS daily 
    dflice - GLORYS dir/sst_file_name.nc
  """
  with xr.open_dataset(dflice) as dsice:
    A2d = dsice['thetao'].isel(time=0, depth=0).values.squeeze()

  SST = A2d[JG,IG]

  return SST

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


def subset_glorys_divu(dfui, dfvi, IG, JG, hlon, hlat, dxy):
  """
    Subset and compute spatial-average divU ice for specified grid points 
    GLORYS daily 
    dflice - GLORYS dir/sst_file_name.nc
    dxy - distance between the subsampled grid points used 
          for deriving linear regr. parameters (km)
  """
  from mod_mom6 import dx_dy  

  DX, DY = dx_dy(hlon, hlat)
  Acell = DX*DY

  with xr.open_dataset(dfui) as dsice:
    A2d = dsice['usi'].isel(time=0).values.squeeze()

  # Treat nans as no ice grid cells
  U2d = np.nan_to_num(A2d, nan=0.0)

  with xr.open_dataset(dfvi) as dsice:
    A2d = dsice['vsi'].isel(time=0).values.squeeze()

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

def derive_time(YS, YE, ptht2m, regn_name, ndays_era):
  """
    Derive time array of available ERA5 fields
    for the full year within YS - YE
  """
  DNMB = None
  time_stmp = []

  print(f"Deriving time array from ERA5 fields")
  # Derive Time from saved atm. fields:
  #ptht2m = DIRS['ptht2m']
  for YR in range(YS,YE+1):
    flnm = f"era5_2mTemp_daily{ndays_era}day_{regn_name}_{YR}.nc"
    dflnm = os.path.join(ptht2m, flnm)
    assert os.path.isfile(dflnm), f"Missing ERA5: {dflnm}, check ndays flag"

    dnmb0 = mtime.datenum([YR,1,1])
    with xr.open_dataset(dflnm, decode_times=False) as ds:
      Time = ds["valid_time"].values
      DYR = dnmb0 + Time
    time_stmp.append(DYR)

  DNMB = np.concatenate(time_stmp)
  return DNMB

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

def subset_era_frzdays(IG, JG, dnmbS, intgr_time, ndays_era, ptht2m, dfgmapi, regn, Tfrz=-1.85):
  """
    Derive dynamic predictor: sqrt of the number of freeze degree days
    Following Zubov's relation: h2 + 50h = 8 IFDD, 
    IFDD = sum of (Tfrz - Tair), when Tair < Tfrz
    For 1 day = dnmbS: integrate (Tfrz-SAT) back to dnmbS-intgr_time
  """
  regions = {
      "north": ("Arctic", 65.0),
      "south": ("Antarctic", -60.0),
  }
  regn_name, lat0 = regions[regn]

  # Load gmapi:
  with xr.open_dataset(dfgmapi) as ds:
    LONE = ds["era_longit"].values
    LATE = ds["era_latit"].values
    IGLR = ds["glorys_indx"].values
    JGLR = ds["glorys_jndx"].values
    IERA = ds["era_indx"].values
    JERA = ds["era_jndx"].values

  LONE = (LONE + 360) % 360

  # Find GLORYS - ERA5 pairs:
  # Alternative to distance-approach, build a lookup table:
  # every glorys JG,IG ---> ERA5 indicex JE,IE
  print("Finding JERA, IERA to match JGLR, IGLR")
  assert len(JG)==len(IG)
  JE = np.empty(len(JG), dtype=int)
  IE = np.empty(len(IG), dtype=int)
  gmapi = {(jg,ig):(je,ie)
           for jg,ig,je,ie in zip(JGLR,IGLR,JERA,IERA)}

  for k, (jj, ii) in enumerate(zip(JG, IG)):
    if k > 0 and k % 50000 == 0:
      prc = k/len(JG)*100.
      print(f"  {prc:.2f}% processed")
    #JE[k], IE[k] = find_era_indx(jj, ii, JERA, IERA, JGLR, IGLR)
    key = (int(jj), int(ii))
    if key not in gmapi:
      raise ValueError(f"Missing gmapi entry for GLORYS index {key}")
    JE[k], IE[k] = gmapi[key]

  # Construct days for integrating Freeze deegre days
  # Note that not every daily fields may be saved
  YR0, MM0, DD0 = mtime.datevec(dnmbS)[:3]
  dnmbP = dnmbS - intgr_time - (ndays_era // 2)  # previous intgr time preiod, start day, add extra
  Ypr, Mpr, Dpr = mtime.datevec(dnmbP)[:3]
  DNMBall = derive_time(Ypr, YR0, ptht2m, regn_name, ndays_era)

  # Find closest time to start - end the integration:
  idxS = max(np.argmin(abs(DNMBall - dnmbP)) - 1, 0)
  idxE = np.argmin(abs(DNMBall - dnmbS)) 
  assert int(DNMBall[idxE]) == int(dnmbS), f"Check idxE - should return {dnmbS}"
  assert dnmbS-DNMBall[idxS] >= intgr_time, \
     f"Integration period is not covered by selected idxS={idxS}"
  DNMB = DNMBall[idxS:idxE+1]

  # Construct 2D array with previous SAT for integration
  npnts = len(JG)
  YRold = 1900
  ds_t2m = None

  SAT = np.zeros((len(IE), len(DNMB)))
  days_frz = np.zeros(len(DNMB))
  print("Deriving SAT for integration ...")
  for irec0, dnmb0 in enumerate(DNMB):
    YR, MM, DD = mtime.datevec(dnmb0)[:3]
    rdate = int(YR*1e4 + MM*100 + DD)
    flt2m = f"era5_2mTemp_daily{ndays_era}day_{regn_name}_{YR}.nc"
    dflt2m = os.path.join(ptht2m, flt2m)
    if YR != YRold:
      if ds_t2m is not None:
          ds_t2m.close()
      YRold = YR
      ds_t2m = xr.open_dataset(dflt2m, decode_times=False)
      Time = ds_t2m['valid_time'].values
      dnmb_day1 = mtime.datenum([YR,1,1])
      TM = (Time + dnmb_day1).astype(int)

    idx = np.where(TM == int(dnmb0))[0]
    if len(idx) == 0:
      raise ValueError(f"No matching day for {dnmb0} {YR}/{MM}/{DD}")
    iday = idx[0]
    A2d = ds_t2m['t2m'].isel(valid_time=iday).values.squeeze()
    T2d = A2d - 273.15  # K --> C

    SAT[:, irec0] = T2d[JE, IE]
    days_frz[irec0] = dnmb0

  assert np.all(days_frz > 0), "days_frz not populated, there are 0s"
  assert np.all(np.diff(days_frz)>0), "days_frz not increasing, required for time interpolation"
  fld_pnts = intgr_Tfrz(SAT, intgr_time, Tfrz, days_frz)

  if ds_t2m is not None:
    ds_t2m.close()

  return fld_pnts

def subset_era_sat(dflt2m, IG, JG, dfgmapi, dnmb0):
  """
    ERA5 2m air temp
  """
  YR, MM, DD = mtime.datevec(dnmb0)[:3]
  # Load gmapi:
  with xr.open_dataset(dfgmapi) as ds:
    LONE = ds["era_longit"].values
    LATE = ds["era_latit"].values
    IGLR = ds["glorys_indx"].values
    JGLR = ds["glorys_jndx"].values
    IERA = ds["era_indx"].values
    JERA = ds["era_jndx"].values

  LONE = (LONE + 360) % 360

  # Find GLORYS - ERA5 pairs:
  # Alternative to distance-approach, build a lookup table:
  # every glorys JG,IG ---> ERA5 indicex JE,IE
  print("Finding JERA, IERA to match JGLR, IGLR")
  assert len(JG)==len(IG)
  JE = np.empty(len(JG), dtype=int)
  IE = np.empty(len(IG), dtype=int)
  gmapi = {(jg,ig):(je,ie)
           for jg,ig,je,ie in zip(JGLR,IGLR,JERA,IERA)}

  for k, (jj, ii) in enumerate(zip(JG, IG)):
    if k > 0 and k % 50000 == 0:
      prc = k/len(JG)*100.
      print(f"  {prc:.2f}% processed")
    #JE[k], IE[k] = find_era_indx(jj, ii, JERA, IERA, JGLR, IGLR)
    key = (int(jj), int(ii))
    if key not in gmapi:
      raise ValueError(f"Missing gmapi entry for GLORYS index {key}")
    JE[k], IE[k] = gmapi[key]

  with xr.open_dataset(dflt2m, decode_times=False) as dsice:
    Time = dsice['valid_time'].values
    dnmbJ1 = mtime.datenum([YR,1,1])
    TM = (Time + dnmbJ1).astype(int)
    idx = np.where(TM == int(dnmb0))[0]
    if len(idx) == 0:
      raise ValueError(f"No matching day for {dnmb0} {YR}/{MM}/{DD}")
    iday = idx[0]
    A2d = dsice['t2m'].isel(valid_time=iday).values.squeeze()

  T2d = A2d - 273.15  # K --> C

  fld_pnts = T2d[JE,IE]

  return fld_pnts


def sens_tests_colors():
  # Line colors:
  CLRS      = np.array([
      [0.00, 0.45, 0.70],  # blue
      [0.90, 0.17, 0.31],  # red
      [0.00, 0.62, 0.38],  # green
      [0.90, 0.60, 0.00],  # orange
      [0.95, 0.90, 0.25],  # yellow
      [0.80, 0.47, 0.65],  # pink/violet
      [0.35, 0.70, 0.90],  # light blue
      [0.50, 0.39, 0.64],  # purple
      [0.55, 0.63, 0.79],  # steel blue
      [0.40, 0.60, 0.00],  # olive green
      [0.70, 0.30, 0.00],  # brown
      [0.70, 0.00, 0.30],  # wine red
      [0.00, 0.55, 0.75],  # teal
      [0.75, 0.75, 0.75],  # light gray
      [0.30, 0.30, 0.30],  # dark gray
      [0.10, 0.80, 0.60],  # aqua green
      [0.55, 0.20, 0.60],  # plum purple
      [0.20, 0.70, 0.30],  # jade green
      [0.80, 0.55, 0.35],  # tan
      [0.25, 0.25, 0.55]   # deep indigo
  ])
  return CLRS

def construct_predictors_day(
        hlon, hlat, IG, JG, dnmb0, DIRS, PRED_NAMES, 
        sqrt_frzdays, sst_max, standz, regn, YS, YE,
        ndays_era, intgr_time, Tfrz,
        order_fast="time", Mavrg=3, dxy=50,
        PRED_MEAN = None, PRED_STDEV = None
   ):
  """
    Construct predictors for 1 day forecast
    PRED_NAMES - predictors used in this model
    predictor array must have exactly the same columns, in exactly the same order,
    as when RF was trained !!!

    hlon, hlat - grid where prediction is done (GLORYS)
    IG, JG - grid point indices where prediction is done
    DIRS - dictionary with input / output directories
    PREAD_NAMES - list of predictors in the right order
    sqrt_frzdays - True: use sqrt of integrated freezing degree days
    sst_max - max sst threshold for possible sea ice
    standz - True: standardize predictors (False for RF)
    regn - region of prediction
    YS, YE - training period, start/end years
    ndays_era - time step in ERA5 atm. fields subsets
    intgr_time - for freeze degree days, integration period, days
    dxy - ice length scale, used for calc. ice predictors and gird point subset

  """
  DNMB = np.asarray([dnmb0])
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  rdate = YR0*10000 + MM0*100 + DD0
  Iconc = None
  SST = None 

  regions = {
      "north": ("Arctic", 65.0),
      "south": ("Antarctic", -60.0),
  }
  regn_name, _ = regions[regn]


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

  # Mean ice thickness over ice area during previous N months:
  # Check Mavrg - should match train_linregr_ithkn.py
  # Interannual trend: mean ice thickness previous N months:
  if 'mnithkn' in PRED_NAMES:
    pthout = DIRS["pthout"]
    fltmp = f"GLORYS_monthly_icevol_ithknmn_{regn}_{YS}_{YE}.npz"
    if YR0 > 2025:
      fltmp = f"GLORYS_monthly_icevol_ithknmn_north_{YR0}_{YR0}.npz"
    dflmni = os.path.join(pthout, fltmp)
    assert os.path.isfile(dflmni), f"File is missing: {dflmni}"
    mnithkn = construct_mean_ithkn(len(IG), DNMB, dflmni, order_fast=order_fast, Mavrg=Mavrg)

    raw['mnithkn'] = mnithkn

  # SST
  SST = None
  if 'sst' in PRED_NAMES:
    print("\nDeriving GLORYS sst")
    pthice = os.path.join(DIRS["pthsst"],f"{YR0}")
    dflice = mglr.find_file(rdate, pthice)
    if dflice is None:
        raise FileNotFoundError(f"Not found {dflice}")
    SST = subset_glorys_sst(dflice, IG, JG)

    raw['sst'] = SST

  # Ice conc
  if 'iconc' in PRED_NAMES:
    print("\nDeriving GLORYS iconc")
    pthice = os.path.join(DIRS["pthiconc"],f"{YR0}")
    dflice = mglr.find_file(rdate, pthice)
    if dflice is None:
        raise FileNotFoundError(f"Not found {dflice}")
    Iconc = subset_glorys_iconc(dflice, IG, JG)

    # Eliminate ice in the warm ocean:
    if SST is not None:
      Iconc[SST > sst_max] = 0.

    raw['iconc'] = Iconc
  # divU ice
  if 'divu' in PRED_NAMES:
    print("\nDeriving GLORYS divu ice")
    pthu = os.path.join(DIRS["pthui"],f"{YR0}")
    dfui = mglr.find_file(rdate, pthu)
    pthv = os.path.join(DIRS["pthvi"],f"{YR0}")
    dfvi = mglr.find_file(rdate, pthv)
    divU = subset_glorys_divu(dfui, dfvi, IG, JG, hlon, hlat, dxy)

    raw['divu'] = divU

  # Freeze days
  if 'frzdays' in PRED_NAMES:
    print("\nDeriving GLORYS Freeze degree days")
    ptht2m = DIRS['ptht2m']
    pthgmapi = DIRS["pthgmapi"]
    flout = f"gmapi_ERA5_to_GLORYS_{regn}.nc"
    dfgmapi = os.path.join(pthgmapi, flout)

    frzdays = subset_era_frzdays(IG, JG, dnmb0, intgr_time, ndays_era,
                                        ptht2m, dfgmapi, regn, Tfrz=Tfrz)
    if sqrt_frzdays:
      frzdays = np.sqrt(frzdays)

    raw['frzdays'] = frzdays

  # SAT
  if 'sat' in PRED_NAMES:
    print("\nDeriving GLORYS SAT")
    flt2m = f"era5_2mTemp_daily{ndays_era}day_{regn_name}_{YR0}.nc"
    ptht2m = DIRS['ptht2m']
    dflt2m = os.path.join(ptht2m, flt2m)

    SAT = subset_era_sat(dflt2m, IG, JG, dfgmapi, dnmb0)

    raw['sat'] = SAT

  # Transform and Standardize predictors
  # Use only active predictors
  # Dict. with standardized arrays in the PRED_NAMES order:
  std_arr = {}
  if standz:
    for name, mu, sigma in zip(PRED_NAMES, PRED_MEAN, PRED_STDEV):
      if np.isnan(mu):
        continue
      std_arr[name] = (raw[name] - mu) / sigma

  # Construct predictors dictionary where each
  # predictor is linked to the standardized or raw array
  # The order of the predictor should match the order in PRED_NAMES
  # Will create dict like this:
  #   'cosD'    : cosD,  # not standardized
  #   'sinD'    : sinD,  # not standardized
  #    ...
  #    'Zcrd'    : std_arr['Zcrd'],  # standardized
  #    'mnithkn' : std_arr['mnithkn'],  # standardized 
  # OR:
  #    'Zcrd'    : raw['Zcrd']  # not standardized
  #     ...
  #
  pred_dict = {}
  for name in PRED_NAMES:
    if standz and name in std_arr:
      pred_dict[name] = std_arr[name]
    else:
      pred_dict[name] = raw[name]

  # Combine all standardized or not standardized predictors into a list:
  PRED_FINAL = []
  PRED_FINAL = [pred_dict[name] for name in PRED_NAMES]  

  return PRED_FINAL, Iconc, SST


