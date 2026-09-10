import numpy as np
import mod_time as mtime
import xarray as xr
import os
from yaml import safe_load
from scipy.interpolate import interp1d
#import sys
#import matplotlib.pyplot as plt
#import importlib

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

def subset_mom6anls_sst(dflsst, IG, JG, idim, jdim):
  """
    Subset SST fields for specified grid points
    from MOM6 analysis fields
    dflsst:     SST dir + file name
    IG, JG :     array-like I- and J-indices of the desired grid points 
    idim, jdim : int Expected X and Y dimensions of the SST field
  """
  with xr.open_dataset(dflsst) as dsst:
    A2d = dsst["Temp"].isel(time=0, z_l=0).values

  assert A2d.shape == (jdim, idim), f"SST shape does not conform {A2d.shape}"

  SST = A2d[JG,IG]

  return SST

def subset_cice6_iconc(dflice, IG, JG, idim, jdim, varnm='aicen'):
  """
    Subset ice conc fields for specified grid points
    from CICE6 analysis fields
    dflice:      ice conc dir + file name
    IG, JG :     array-like I- and J-indices of the desired grid points 
    idim, jdim : int Expected X and Y dimensions of the SST field
  """
  with xr.open_dataset(dflice) as dsice:
    aicen = dsice[varnm].values

  # Aggregate over categories:
  A2d = np.nansum(aicen, axis=0)

  # Treat nans as no ice grid cells
  A2d = np.nan_to_num(A2d, nan=0.0)

  assert A2d.shape == (jdim, idim), f"SST shape does not conform {A2d.shape}"

  Iconc = A2d[JG,IG]
  Iconc[Iconc > 1] = 1.
  Iconc[Iconc < 0] = 0.

  return Iconc


def subset_gdas_frzdays(dfifdd, IG, JG, idim, jdim):
  """
  Extract IFDD at selected mesh025 grid points.

  IFDD fields are derived from daily-averaged GDAS T2m fields
  and interpolated onto the mesh025 grid.

  See: derive_intgrFDD_mesh025
  """

  DATA = np.load(dfifdd)

  IFDDs = DATA["IFDD"]   # IFDD at selected grid points in polar region
  JP    = DATA["JG"]
  IP    = DATA["IG"]

  jdim_saved = DATA["jdim"]
  idim_saved = DATA["idim"]

  assert idim_saved == idim and jdim_saved == jdim, (
      f"IFDD grid {jdim_saved} x {idim_saved} "
      f"does not match jdim={jdim}, idim={idim}"
  )

  # Insert saved IFDD values back onto the full mesh025 grid
  A2d = np.full((jdim, idim), np.nan)
  A2d[JP, IP] = IFDDs

  # Extract the requested ML grid points
  IFDD = A2d[JG, IG]

  assert np.all(np.isfinite(IFDD)), (
      "Some requested points are outside the saved IFDD domain or contain NaN"
  )

  return IFDD

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
 

def construct_predictors_day(
        hlon, hlat, IG, JG, dnmb0, fyaml, PRED_NAMES,
        sqrt_frzdays, sst_max, standz, regn,
        intgr_time, Tfrz,
        order_fast="time", Mavrg=3, dxy=50, 
        PRED_MEAN = None, PRED_STDEV = None
   ):
  """
    Construct predictors for 1 day forecast
    PRED_NAMES - predictors used in this model
    predictor array must have exactly the same columns, in exactly the same order,
    as when RF was trained !!!
    
    hlon, hlat - grid where prediction is done (mesh025)
    IG, JG - grid point indices where prediction is done
    MLYAML - YAML dictionary with input / output directories, etc.
    PREAD_NAMES - list of predictors in the right order
    sqrt_frzdays - True: use sqrt of integrated freezing degree days
    sst_max - max sst threshold for possible sea ice
    standz - True: standardize predictors (False for RF)
    regn - region of prediction
    intgr_time - for freeze degree days, integration period, days
    dxy - ice length scale, used for calc. ice predictors and gird point subset
  
  """

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

  # Mean ice thickness over ice area during previous N months:
  # Check Mavrg - should match train_linregr_ithkn.py
  # Interannual trend: mean ice thickness previous N months
  # Use GLORYS mean ice thicknesses, SOCA ithkn is too bad
  if 'mnithkn' in PRED_NAMES:
    pthglr = MLYAML["PRED"]["pthglr"]
    if YR0 < 2026:
      fltmp = f"GLORYS_monthly_icevol_ithknmn_{regn}_1993_2025.npz"
    elif YR0 == 2026: 
      fltmp = f"GLORYS_monthly_icevol_ithknmn_north_{YR0}_{YR0}.npz"

    dflmni = os.path.join(pthglr, fltmp)
    assert os.path.isfile(dflmni), f"File is missing: {dflmni}"
    mnithkn = construct_mean_ithkn(len(IG), DNMB, dflmni, order_fast=order_fast, Mavrg=Mavrg)
      
    raw['mnithkn'] = mnithkn
      
  # SST
  SST = None
  if 'sst' in PRED_NAMES:
    # SOCA MOM6 rdate, 6hr before the f/cast time for IAU:
    fhr = 6  # IAU time window
    dnmb_mom = dnmb0 - 6/24
    YRm, MMm, DDm, HRm = mtime.datevec(dnmb_mom)[:4]
    mom_rdate = int(YRm)*10000 + int(MMm)*100 + int(DDm)
    print("\nDeriving MOM6 sst")
    pthsst = MLYAML["PRED"]["pthsst"].format(mom_rdate=mom_rdate)
    flsst  = MLYAML["PRED"]["flsst"].format(hr=int(HRm), fhr=fhr)
    dflsst = os.path.join(pthsst, flsst)

    if dflsst is None:
        raise FileNotFoundError(f"MOM6 analysis SST Not found {dflsst}")
    SST = subset_mom6anls_sst(dflsst, IG, JG, idim, jdim)

    raw['sst'] = SST

  # Ice conc
  if 'iconc' in PRED_NAMES:
    print("\nDeriving SOCA CICE6 iconc")
    varnm = 'aicen'
    pthice = MLYAML["PRED"]["pthiconc"].format(YR=YR0, MM=MM0, DD=DD0)
    ficonc = MLYAML["PRED"]["fliconc"].format(YR=YR0, MM=MM0, DD=DD0)
    dflice = os.path.join(pthice, ficonc)

    if dflice is None or not os.path.isfile(dflice):
        raise FileNotFoundError(f"Not found {dflice}")

    Iconc = subset_cice6_iconc(dflice, IG, JG, idim, jdim)

    # Eliminate ice in the warm ocean:
    if SST is not None:
      Iconc[SST > sst_max] = 0.

    raw['iconc'] = Iconc

  # divU ice
  if 'divu' in PRED_NAMES:
    print("\nDeriving SOCA CICE6 divu ice")
    pthui = MLYAML["PRED"]["pthui"].format(YR=YR0, MM=MM0, DD=DD0)
    flui  = MLYAML["PRED"]["flui"].format(YR=YR0, MM=MM0, DD=DD0)
    dflui = os.path.join(pthui, flui)

    if not os.path.isfile(dflui):
        raise FileNotFoundError(f"Not found {dflui}")

    divU = subset_cice6_divu(dflui, IG, JG, hlon, hlat, dxy)

    raw['divu'] = divU

  # Freeze days
  if 'frzdays' in PRED_NAMES:
    print("\nDeriving GDAS GFSv17 Freeze degree days")
    pthifdd = MLYAML["PRED"]["pthifdd"]
    flifdd  = MLYAML["PRED"]["flifdd"].format(Ndays=intgr_time, dxy=dxy, rdate=rdate, regn=regn)
    flnpz = f"{flifdd}.npz"
    dfifdd = os.path.join(pthifdd, flnpz)

    if not os.path.isfile(dfifdd):
        raise FileNotFoundError(f"Not found {dfifdd}")

    frzdays = subset_gdas_frzdays(dfifdd, IG, JG, idim, jdim)

    if sqrt_frzdays:
      frzdays = np.sqrt(frzdays)

    raw['frzdays'] = frzdays

  # SAT
  if 'sat' in PRED_NAMES:
    print("\nDeriving GDAS T2m")
    ptht2m = MLYAML["PRED"]["pthit2m"].format(regn=regn, YR=YR0)
    flt2m = MLYAML["PRED"]["flt2m"].format(YR=YR0, MM=MM0, DD=DD0, regn=regn)
    dflt2m = os.path.join(ptht2m, flt2m)

    if not os.path.isfile(dflt2m):
      raise FileNotFoundError(f"Not found {dflt2m}")

    SAT = subset_gdas_t2m(dflt2m, IG, JG, idim, jdim)

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





