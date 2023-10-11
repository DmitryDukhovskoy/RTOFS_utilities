"""
  Functions for Gulfstream analysis

  Dmitry Dukhovskoy NOAA/NCEP/EMC April 2023

"""
import os
import numpy as np
import sys
import importlib
import matplotlib.pyplot as plt
import datetime

def derive_contour(AA, tz0=12., xFS=2565, yFS=1809):
  """
    AA is T field with regions that are outside the 
   study region being masked out
   index space - for 0.08 global RTOFS tripolar grid 4500 x 3298 
  """

#  plt.ion()
  figA = plt.figure(10,figsize=(8,8))
  plt.clf()

  axA1 = plt.axes([0.1, 0.2, 0.7, 0.7])
  CS   = axA1.contour(AA,[tz0])

  axA1.axis('equal')
  axA1.set_xlim([2520,3050])
  axA1.set_ylim([1800,2250])

  SGS  = CS.allsegs[0]  # should be only 1 contoured value
  nsgs = len(SGS)

# Delete all closed contours
  CNTR = []  
  for isg in range(nsgs):
    XY = SGS[isg]
    X  = XY[:,0]
    Y  = XY[:,1]

    dEnd = np.sqrt((X[0]-X[-1])**2+(Y[0]-Y[-1])**2)
    if dEnd < 1.:
      continue

    CNTR.append(XY)

# Aranage all segments in order
# First find segment that starts in Fl. Str
  nC  = len(CNTR)
  TCNT = []

  for ii in range(nC):
    if ii == 0:
      x0 = xFS
      y0 = yFS
    else:
      x0 = TCNT[-1,0]
      y0 = TCNT[-1,1]

    xsgm, ysgm, imin, jmin = arange_1segm(CNTR,x0,y0)  

# Remove selected segment:
    CNTR.pop(imin)

    if ii == 0:
      TCNT = np.transpose(np.array((xsgm,ysgm)))
    else:
      aa   = np.transpose(np.array((xsgm,ysgm)))
      TCNT = np.append(TCNT, aa, axis=0)


#X = TCNT[:,0]
#Y = TCNT[:,1]
#axA1.plot(X,Y,'.-') 
  plt.close(figA)

  return TCNT

def adjust_rtofs_gsnw(Ir, Jr, In, Jn, nskip=0):
  """
  Adjust RTOFS Gulf Stream north wall contour
  to NAVO contour by choping off segments 
  that are outside the NAVO contour

  Also, if the contour has too many points
  skip nskip points along the contour

  if nskip < 0 - choose nskip based on NAVO contour length
  to make rtofs contour approximately same # of points
  """
  in0   = In[0]
  jn0   = Jn[0]
  dd    = np.sqrt((Ir-in0)**2 + (Jr-jn0)**2)
  icutS = np.where(dd == min(dd))[0][0]
  in0   = In[-1]
  jn0   = Jn[-1]
  dd    = np.sqrt((Ir-in0)**2 + (Jr-jn0)**2)
  icutE = np.where(dd == min(dd))[0][0]
  lnavo = np.shape(In)[0]

  if nskip == 0:
    Irr = Ir[icutS:icutE+1]
    Jrr = Jr[icutS:icutE+1]
  elif nskip > 0:
    nstp = nskip+1
    Irr = Ir[icutS:icutE+nskip:nstp]
    Jrr = Jr[icutS:icutE+nskip:nstp]
  else:
    npnts = icutE-icutS+1
    nskip = int(npnts/lnavo)
    nstp  = nskip+1
    Irr = Ir[icutS:icutE+nskip:nstp]
    Jrr = Jr[icutS:icutE+nskip:nstp]

  return Irr, Jrr

def arange_1segm(CNTR,x0,y0):
  """
    Find segment closest to x0, y0 
    arange the orientation of the segment
    so that it starts from the pnt closest to x0, y0
    CNTR - list with segments X,Y as np arrays
  """
  nC  = len(CNTR)
  DFS = np.zeros((nC,2))*np.nan
  for isg in range(nC):
    XY = CNTR[isg]
    X  = XY[:,0]
    Y  = XY[:,1]

    d1 = np.sqrt((X[0]-x0)**2+(Y[0]-y0)**2)
    d2 = np.sqrt((X[-1]-x0)**2+(Y[-1]-y0)**2)

    DFS[isg,0] = d1
    DFS[isg,1] = d2   
 
  imin = np.argmin(np.min(DFS, axis=1))
  jmin = np.argmin(np.min(DFS, axis=0))

  xsgm = CNTR[imin][:,0]
  ysgm = CNTR[imin][:,1]

  if jmin == 1:
    xsgm = np.flip(xsgm)
    ysgm = np.flip(ysgm)

  return xsgm, ysgm, imin, jmin


def str2month(sMO):
  """
    Convert string months to numbers
  """
  mo2num = {
    "JAN" : 1,
    "FEB" : 2,
    "MAR" : 3, 
    "APR" : 4, 
    "MAY" : 5, 
    "JUN" : 6, 
    "JUL" : 7,
    "AUG" : 8,
    "SEP" : 9,
    "OCT" : 10,
    "NOV" : 11, 
    "DEC" : 12
  }

  try:
    MM = mo2num[sMO]
  except:
    print("Converting NAVO month {0} to number failed".format(sMO))

  return MM
 

def read_navogs(rdate, pthnavo, missing=1):
  """
    Read Gulfstream northwall position from the NAVO forecast
    From Gwen:
    On Hera, you can find the NAVO Gulf Stream north wall data 
    files at 
    ##/scratch2/NCEPDEV/ovp/Lichuan.Chen/DCOM/$YYYYMMDD/wtxtbul.  
    /scratch2/NCEPDEV/ovp/Samira.Ardani/DCOM
    We use the gs*.sub to generate the plots at 
    https://polar.ncep.noaa.gov/global/fronts/.  
    You can also find the data files on WCOSS2.  
    They are located at /lfs/h1/ops/prod/dcom/$YYYYMMDD/wtxtbul 
    (e.g., /lfs/h1/ops/prod/dcom/20230407/wtxtbul/gs091nw.sub).  

    INPUT:
      rdate - YYYMMDD date of the Gulfstream NW position
              Forecast position is 1 day f/ward

      missing - # of days searching backward/forward
                if the forecast is missing for this date

  """
  import mod_misc1 as mmisc  
  import mod_time as mtime
  
  rnmb    = mtime.rdate2datenum(rdate)
#  rdnavo     = mtime.adddays_date(rdate,-1)
#  navonmb = mtime.rdate2datenum(rdnavo)
  navonmb0 = rnmb - 1
  rdnavo0  = mtime.dnumb2rdate(navonmb0, ihours=False)


  gsFound = False
  lfind = [0]
  for ii in range(1,missing+1):
    lfind.append(ii)
    lfind.append(-ii)

  for ii in range(len(lfind)):
    if gsFound:
      break
    ifnd = lfind[ii]
    navonmb = navonmb0 + ifnd
    rdnavo  = mtime.dnumb2rdate(navonmb, ihours=False)
    pthdat  = pthnavo + rdnavo + '/wtxtbul/'
    if not os.path.isdir(pthdat):
      print('NAVO dir is missing ' + pthdat)
      continue
# Find file:
# Look for gs*.sub file with Gulf Stream front:
    for ffl in os.listdir(pthdat):
      if ffl.startswith('gs') and ffl.endswith('sub'):
        fpnavo = pthdat + ffl
        print(rdnavo + ' ok ')
        gsFound = True
        break

  if not gsFound:
    print('NAVO gs*.sub file not found for {1} +/- {0}'.format(missing,rdanvo0))
    raise Exception(" NAVO not found, try increase missing search range")
      
  fid = open(fpnavo,'r')
  fpos = fid.tell()
  fid.seek(0,2)
  fend = fid.tell()
# Read header then data
  fid.seek(0)

  XNW = []
  YNW = []
  fdata = False
  while fpos < fend:
    dmm = fid.readline().split()
    if len(dmm)==0:
      continue
    elif len(dmm[0].strip()) == 0:
      print('Empty string, skipping ...')
      continue

    smm = str(dmm[0])
    fpos = fid.tell()
    if fpos == fend:
      break

    nss = '**'
    nnm = '**'
    if len(dmm) >= 9:
      nnm = str(dmm[3])+str(dmm[4])
      nss = str(dmm[2])+str(dmm[3])      

    if nss == 'SOUTHWALL':
      fdata = False
      break    

    if nnm == 'NORTHWALL':
      fdata = True
      MD  = int(str(dmm[7]))
      sMO = str(dmm[8])
      smm = str(dmm[9])
      YR  = 2000+int(smm[0:2])
      MM  = str2month(sMO)
      fcstnmb = mtime.datenum([YR,MM,MD]) 
      continue

    if fdata:
      nrc = len(dmm)
      for irc in range(nrc):
        aa  = dmm[irc]
        iN  = aa.index('N')
        iW  = aa.index('W')
        lat = float(aa[0:iN])
        lon = -abs(float(aa[iN+1:iW]))

        XNW.append(lon)
        YNW.append(lat)

    
  fid.close()

  return XNW, YNW, fcstnmb
    

