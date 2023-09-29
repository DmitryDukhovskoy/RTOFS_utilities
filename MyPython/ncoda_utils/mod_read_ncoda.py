"""
  Read NCODA binary
  
  Dmitry Dukhovskoy, NOAA NCEP
"""
import os
from netCDF4 import Dataset as ncFile
import numpy as np
import sys

def read_ncoda_inc(fina,IDM,JDM,KDM,rLayer=None):
  """
    Read NCODA increment files used in ncoda_archv_inc.f
    files are binary stream access
    rLayer - read 1 layer, otherwise read all KDM
  """
  IJDM = IDM*JDM
#  fgb.seek(0)


  print('Reading HYCOM :{0} '.format(fina))
  print('Grid: IDM={0}, JDM={1}'.format(IDM,JDM))
  try: 
    fga = open(fina,'rb')
  except:
    raise Exception('Could not open '+fina)

  if rLayer:
    lr1 = rLayer
    lr2 = lr1
  else:
    lr1 = 1
    lr2 = KDM

  F = []
  ccL = -1
  for ii in range(lr1,lr2+1):
    print('Reading layer {0} '.format(ii))
    fga.seek(0)
    k0 = ii-1
    fga.seek(k0*IJDM*4,0)
    dmm = np.fromfile(fga, dtype='>f',count=IJDM) # read 1 layer
    dmm = dmm.reshape((JDM,IDM))
    ccL += 1
#   print('lr={0}, ccL={1}'.format(ii,ccL))
#   breakpoint()
    if ccL == 0:
      F = np.copy(dmm)
    else:
      F = np.dstack((F,dmm))


  fga.close()
#
# Reshape 3D array:
# Make layers in 0 axis
  if len(np.shape(F)) == 3 and F.shape[2] == KDM:
    F = np.transpose(F, (2,0,1))
    
  return F


def ncoda_depths(zneg=True):
  """
    NCODA interface and mid-point depths
    zi = from hycom/ncoda_archv_inc/zi.txt
  """
  print(' Reading NCODA depths')
  ZI = np.array([
      0.00,
      1.00,
      3.00,
      5.00,
     11.00,
     21.00,
     35.00,
     53.00,
     67.00,
     85.00,
     99.00,
    117.00,
    131.00,
    149.00,
    171.00,
    189.00,
    211.00,
    229.00,
    251.00,
    269.00,
    291.00,
    309.00,
    371.00,
    389.00,
    451.00,
    469.00,
    531.00,
    589.00,
    651.00,
    709.00,
    771.00,
    829.00,
    971.00,
   1029.00,
   1271.00,
   1329.00,
   1571.00,
   1629.00,
   1871.00,
   2129.00,
   2371.00,
   2629.00])

  kzi = ZI.shape[0]
#
# Follow ncoda_archv_inc.f to define the mid-depths points in NCODA
# Note that in the ncoda code zz is NCODA z-level mid-depths points
# and zi is array with z-level interface depths
  ZM = np.zeros((kzi-1))
  for kk in range(1,kzi-1):
    ZM[kk] = 0.5*(ZI[kk-1]+ZI[kk])

  kzm = ZM.shape[0]

  if zneg:
    ZI = -ZI
    ZM = - ZM

  return ZI, kzi, ZM, kzm


def anls2bckground_date(rdate):
  """
  Define background date given analysis date as a string
  2022061600 or 20220616 then 00 hrs assumed
  """
  import datetime

  ll = len(rdate)
  yr = int(rdate[0:4])
  mo = int(rdate[4:6])
  md = int(rdate[6:8])
  if ll<10:
    hr = 0
  else:
    hr = int(rdate[8:10])

  time_anls = datetime.datetime(yr,mo,md,hr,0,0)
  time_bgrd = time_anls - datetime.timedelta(days=1)
  rbgrd     = time_bgrd.strftime('%Y%m%d%H')

  return rbgrd

def adddays_date(rdate,ndays):
  """
  Add/subtract n days from rdate
  """
  import datetime

  ll = len(rdate)
  yr = int(rdate[0:4])
  mo = int(rdate[4:6])
  md = int(rdate[6:8])
  if ll<10:
    hr = 0
  else:
    hr = int(rdate[8:10])

  time_anls = datetime.datetime(yr,mo,md,hr,0,0)
  time_bgrd = time_anls + datetime.timedelta(days=ndays)
  rbgrd     = time_bgrd.strftime('%Y%m%d%H')

  return rbgrd


def rdate2julian(rdate):
  """
  Given string with rdate derive year day
  """
  import datetime

  ll = len(rdate)
  yr = int(rdate[0:4])
  mo = int(rdate[4:6])
  md = int(rdate[6:8])

  time_anls = datetime.datetime(yr,mo,md,0,0,0)
  time_jan1 = datetime.datetime(yr,1,1,0,0,0)

  jday = (time_anls-time_jan1).days+1

  return jday

def datenum(ldate0,ldate_ref=[1,1,1]):
  """
  Given list [YY,MM,DD] - current date 
  compute days wrt to reference date - optional
  """
  import datetime

  ll = len(ldate0)
  YR = ldate0[0]
  MM = ldate0[1]
  DD = ldate0[2]
  if ll > 3:
    HH = ldate0[3]
    MN = ldate0[4]
  else:
    HH = 0
    MN = 0

  lr = len(ldate_ref)
  YRr = ldate_ref[0]
  MMr = ldate_ref[1]
  DDr = ldate_ref[2]
  if lr > 3:
    HHr = ldate_ref[3]
    MNr = ldate_ref[4]
  else:
    HHr = 0
    MNr = 0


  time0 = datetime.datetime(YR,MM,DD,HH,MN,0)
  timeR = datetime.datetime(YRr,MMr,DDr,HHr,MNr,0)

  dnmb = (time0-timeR).days+1+(HH-HHr)/24.+(MN-MNr)/1440.

  return dnmb


def parse_rdate(rdate):
  """
  Derive date from rdate string
  """

  ll = len(rdate)
  yr = int(rdate[0:4])
  mo = int(rdate[4:6])
  md = int(rdate[6:8])
  if ll<10:
    hr = 0
  else:
    hr = int(rdate[8:10])

  return yr, mo, md, hr

def date_type_output(rdate,toutp):
  """
  Define date of the saved output given forecast date D
  and type of saved output: 
  incup - incremental update run = D-1-6 hrs from f/cast time
  anls  - hindcast analysis run (note: HYCOM is free run, atm. forc - anlysis)
          D-1:00hr - D-1:24hr
  fcast - true forecast = D
  """
  import datetime

  ll = len(rdate)
  yr = int(rdate[0:4])
  mo = int(rdate[4:6])
  md = int(rdate[6:8])
  if ll<10:
    hr = 0
  else:
    hr = int(rdate[8:10])

  time_fcst = datetime.datetime(yr,mo,md,hr,0,0)
  if toutp == "fcast":
    time_out = time_fcst
  elif toutp == "incup":
    time_out = time_fcst - datetime.timedelta(hours=30)
  elif toutp == "fcast":
    time_out = time_fcst + datetime.timedelta(days=1)

  date_outp = time_out.strftime('%Y%m%d%H')

  return date_outp
  



  



