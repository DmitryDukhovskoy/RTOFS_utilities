import numpy as np
import mod_time as mtime
import xarray as xr
import os
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

