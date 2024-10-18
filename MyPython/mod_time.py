"""
  Functions computing calendar days, time conversions 
"""
import datetime
import time
import numpy as np

def datenum(ldate0,ldate_ref=[1,1,1,0,0]):
  """
  Given list [YY,MM,DD] - current date 
  compute days wrt to reference date - optional
  Hours and Minutes  - optional
  [YY,MM,DD,HR]
  [YY,MM,DD,HR,MN]
  """

  ll = len(ldate0)
  YR = ldate0[0]
  MM = ldate0[1]
  DD = ldate0[2]
  HR = 0
  MN = 0
  if ll == 4:
    HR = ldate0[3]
    MN = 0
  elif ll == 5:
    HR = ldate0[3]
    MN = ldate0[4]

  lr = len(ldate_ref)
  YRr = ldate_ref[0]
  MMr = ldate_ref[1]
  DDr = ldate_ref[2]
  HRr = 0
  MNr = 0
  if lr == 4:
    HRr = ldate_ref[3]
    MNr = 0
  elif lr == 5:
    HRr = ldate_ref[3]
    MNr = ldate_ref[4]

  YR = int(YR)
  MM = int(MM)
  DD = int(DD)
  HR = int(HR)
  MN = int(MN)
  time0 = datetime.datetime(YR,MM,DD,HR,MN,0)
  timeR = datetime.datetime(YRr,MMr,DDr,HRr,MNr,0)

  dnmb = float((time0-timeR).days)+1.+(HR-HRr)/24.+(MN-MNr)/1440.

  return dnmb

def adddays_date(rdate,ndays):
  """
  Add/subtract n days from rdate
  rdate is in the format YYYYMMDD[HR]
  """

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

def jday2dnmb(YR, jday, ldate_ref=[1,1,1]):
  """
    Convert Year and year day (jday) to date number dnmb
    wrt to ref. date ldate_ref

    type of YR and jday should be the same (int, list or array)
  """
#  if isinstance(YR, float): YR=int(YR)
#  if isinstance(jday, float): jday=int(jday)
#  print(f'YR={YR} jday={jday}')
#  print(type(YR))
  if isinstance(YR, list): YR = np.array(YR)
  if isinstance(jday, list): jday = np.array(jday)
  if isinstance(YR, np.generic): 
    YR = YR.item()
  if isinstance(jday, np.generic): 
    jday =  jday.item()
#  print(type(YR))

  if isinstance(YR, np.ndarray):
    nyr = len(YR)
    dnmb = np.zeros((nyr))
    for ii in range(nyr):
#      print(f'ii={ii} {YR[ii]} {jday[ii]}')
      dnmb0 = datenum([YR[ii],1,1])
      dnmb[ii]  = dnmb0 + jday[ii]-1
  else:
    dnmb0 = datenum([YR,1,1])
    dnmb  = dnmb0 + jday-1
      
  return dnmb 

def dnmb2jday(dnmb):
  """
    Obtain year day from a date number
  """
  DV = datevec(dnmb)
  YR = DV[0]
  dnmb0 = datenum([YR,1,1])
  jday = dnmb - dnmb0 + 1

  return YR, jday


def rdate2jday(rdate):
  """
  Given string with rdate YYYYMMDD or YYYYMMDDHH derive year day
  """

  ll = len(rdate)
  yr = int(rdate[0:4])
  mo = int(rdate[4:6])
  md = int(rdate[6:8])

  time_anls = datetime.datetime(yr,mo,md,0,0,0)
  time_jan1 = datetime.datetime(yr,1,1,0,0,0)

  jday = (time_anls-time_jan1).days+1

  return jday

def date2jday(ldate0):
  """
  Given list [YY,MM,DD] - current date 
  compute year day
  Hours and Minutes  - optional
  [YY,MM,DD,HR]
  [YY,MM,DD,HR,MN]
  """
  YR     = ldate0[0]
  dnmbJ1 = datenum([YR,1,1]) 
  dnmb0  = datenum(ldate0)
  jday   = dnmb0 - dnmbJ1 + 1.

  return jday
 

def parse_rdate(rdate):
  """
  Derive date fields from rtofs string
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


def rdate2date(rdate):
  """
  Convert rtofs date YYYYMMDD or YYYYMMDDHR to YY, MM, DD, HR
  return list
  """
  YR     = int(rdate[0:4])
  MM     = int(rdate[4:6])
  DD     = int(rdate[6:8])

  if len(rdate) > 9:
    HR = int(rdate[8:10])
  else:
    HR = 0

  Ldate = [YR,MM,DD,HR]
  return Ldate


def rdate2datenum(rdate):
  """
  Convert rtofs date YYYYMMDD or YYYYMMDDHR to
  matlab-type datenum
  """
  YR     = int(rdate[0:4])
  MM     = int(rdate[4:6])
  DD     = int(rdate[6:8])

  if len(rdate) > 9:
    HR = int(rdate[8:10])
  else:
    HR = 0

  Ldate = [YR,MM,DD,HR]
  dnmb = datenum(Ldate)

  return dnmb

def datevec(dnmb,ldate_ref=[1,1,1]):
  """
  For datenum computed wrt to reference date - see datenum
  convert datenum back to [YR,MM,DD,HR,MN]
  """
  if isinstance(dnmb, np.generic): 
    dnmb = dnmb.item()

  if not (isinstance(dnmb, int) or isinstance(dnmb, float)):
    raise Exception('dnmb should be int or float, for array use datevec2D')
  lr = len(ldate_ref)
  YRr = ldate_ref[0]
  MMr = ldate_ref[1]
  DDr = ldate_ref[2]
  if lr > 3:
    HRr = ldate_ref[3]
    MNr = ldate_ref[4]
  else:
    HRr = 0
    MNr = 0

  timeR = datetime.datetime(YRr,MMr,DDr,HRr,MNr,0)
  dfrct = dnmb-np.floor(dnmb)
  if abs(dfrct) < 1.e-6:
    HR = 0
    MN = 0
  else:
    HR = int(np.floor(dfrct*24.))
    MN = int(np.floor(dfrct*1440.-HR*60.))

  ndays = int(np.floor(dnmb))-1
  time0 = timeR+datetime.timedelta(days=ndays, seconds=(HR*3600 + MN*60))
  YR = time0.year
  MM = time0.month
  MD = time0.day
  HR = time0.hour
  MN = time0.minute

  dvec = [YR,MM,MD,HR,MN]

  return dvec


def datevec2D(DNMB0,ldate_ref=[1,1,1]):
  """
  For datenum computed wrt to reference date - see datenum
  convert datenum back to [YR,MM,DD,HR,MN]

  dnmb is an array for N dates: [dnmb(1), dnmb(2), ...]

  output: DV - 2D array
  [YR1, MM1, DD1],
  [YR2, MM2, DD2], ...
                                 
  """

  lr = len(ldate_ref)
  YRr = ldate_ref[0]
  MMr = ldate_ref[1]
  DDr = ldate_ref[2]
  if lr > 3:
    HRr = ldate_ref[3]
    MNr = ldate_ref[4]
  else:
    HRr = 0
    MNr = 0

  nrec = DNMB0.shape[0]
  DV = np.zeros((nrec,5), dtype=int)
  for irec in range(nrec):
    dnmb  = DNMB0[irec]
    timeR = datetime.datetime(YRr,MMr,DDr,HRr,MNr,0)
    dfrct = dnmb-np.floor(dnmb)

    if abs(dfrct) < 1.e-6:
      HR = 0
      MN = 0
    else:
      HR = int(np.floor(dfrct*24.))
      MN = int(np.floor(dfrct*1440.-HR*60.))

    ndays = int(np.floor(dnmb))-1
    time0 = timeR+datetime.timedelta(days=ndays, seconds=(HR*3600 + MN*60))
    YR = time0.year
    MM = time0.month
    MD = time0.day
    HR = time0.hour
    MN = time0.minute

    dvec = np.array([YR,MM,MD,HR,MN], dtype=int)
    DV[irec,:] = dvec

  return DV


def datevec1D(dnmb,ldate_ref=[1,1,1]):
  """
  For datenum computed wrt to reference date - see datenum
  convert datenum back to [YR,MM,DD,HR,MN]
  Input is 1D numpy array

  """

  lr = len(ldate_ref)
  YRr = ldate_ref[0]
  MMr = ldate_ref[1]
  DDr = ldate_ref[2]
  if lr > 3:
    HRr = ldate_ref[3]
    MNr = ldate_ref[4]
  else:
    HRr = 0
    MNr = 0

  timeR = datetime.datetime(YRr,MMr,DDr,HRr,MNr,0)
  dfrct = dnmb-np.floor(dnmb)
  HRi   = np.floor(dfrct*24.).astype(int)
  MNi   = np.floor(dfrct*1440.-HRi*60.).astype(int)

  ndays = (np.floor(dnmb)-1).astype(int)

  YR = []
  MM = []
  MD = []
  HR = []
  MN = []
  for it in range(np.shape(ndays)[0]):
    time0 = timeR+datetime.timedelta(days=ndays.item(it), \
                   seconds=(HRi.item(it)*3600 + MNi.item(it)*60))
    YR.append(time0.year)
    MM.append(time0.month)
    MD.append(time0.day)
    HR.append(time0.hour)
    MN.append(time0.minute)

  YR = np.array(YR)
  MM = np.array(MM)
  MD = np.array(MD)
  HR = np.array(HR)
  MN = np.array(MN)
  dvec = [YR,MM,MD,HR,MN]

  return dvec

def datestr(dnmb,ldate_ref=[1,1,1], show_hr=True):
  """
  For datenum computed wrt to reference date - see datenum
  convert datenum back to [YR,MM,DD,HR,MN]
  print the date
  """
# Check if this is an array or a scalar:
  if hasattr(dnmb, '__len__'):
    farray = True
  else:
    farray = False

  lr = len(ldate_ref)
  YRr = ldate_ref[0]
  MMr = ldate_ref[1]
  DDr = ldate_ref[2]
  if lr > 3:
    HRr = ldate_ref[3]
    MNr = ldate_ref[4]
  else:
    HRr = 0
    MNr = 0

  timeR = datetime.datetime(YRr,MMr,DDr,HRr,MNr,0)

  def get_dstr(dnmb):
    dfrct = dnmb-np.floor(dnmb)
    if abs(dfrct) < 1.e-6:
      HR = 0
      MN = 0
    else:
      HR = int(np.floor(dfrct*24.))
      MN = int(np.floor(dfrct*1440.-HR*60.))

    ndays = int(np.floor(dnmb))-1
    time0 = timeR+datetime.timedelta(days=ndays, seconds=(HR*3600 + MN*60))

    if show_hr:
      dstr = time0.strftime('%Y/%m/%d %H:%M')
    else:
      dstr = time0.strftime('%Y/%m/%d')
    
    return dstr

  if farray:
    nrec = len(dnmb)
    DSTR = []
    for irec in range(nrec):
      dstr = get_dstr(dnmb[irec])
      DSTR.append(dstr)
  else:
    DSTR = get_dstr(dnmb)

  return DSTR


def dnumb2rdate(dnmb, ihours=True):
  """
  Convert mat-type date number to rtofs date string
  YYYYMMDD or YYYYMMDDHR if ihours
  """

  dvec = datevec(dnmb)
#  dvec = [YR,MM,MD,HR,MN
  YR = dvec[0]
  MM = dvec[1]
  MD = dvec[2]
  HR = dvec[3]
  if ihours:
    rdate_out = '{0:4d}{1:02d}{2:02d}{3:02d}'.format(YR,MM,MD,HR)
  else:
    rdate_out = '{0:4d}{1:02d}{2:02d}'.format(YR,MM,MD)

  return rdate_out 
  
def dnumbTrue(dnmb0,sfx):
  """
    From RTOFS date - forecast date at n00
    determine the actual date/time of the output
    using sfx: n-30, ..., n-24, ..., n00, f01, ....
    n-24 - incrementally updated fields
    n-24 - n00 - "hindcast" RTOFS forced with atm analysis
    n00 - initial fields of the actual forecast
    f?? - forecasts
  """
  if sfx == 'n-24':
    dnmbP = dnmb0-1
  elif sfx[0] == 'f':
    hr = int(sfx[1:])
    dnmbP = dnmb0+float(hr)/24.

  return dnmbP

def month_days(imo, YR):
  """
    Define the number of days in a month
  """
  dnmb1 = datenum([YR,imo,1])
  dv2   = datevec(dnmb1 + 32)
  imoN  = dv2[1]
  YRN   = dv2[0]
  dnmb2 = datenum([YRN, imoN,1])
  ndays = int(dnmb2-dnmb1)

  return ndays


