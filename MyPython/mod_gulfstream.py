"""
  Functions/ subroutines for Gulfstream analysis

  Dmitry Dukhovskoy NOAA/NCEP/EMC 

"""
import os
import numpy as np
import sys
import importlib
import matplotlib.pyplot as plt
import datetime

def derive_contour(AA, tz0=12., xFS=2565, yFS=1809, \
                   xl1=2520, yl1=1800, xl2=3050, yl2=2250):
  """
    AA is T field with regions that are outside the 
   study region being masked out
   index space - for 0.08 global RTOFS tripolar grid 4500 x 3298 
   xFS, yFS - control point in Florida Straits
  """

#  plt.ion()
  figA = plt.figure(10,figsize=(8,8))
  plt.clf()

  axA1 = plt.axes([0.1, 0.2, 0.7, 0.7])
  CS   = axA1.contour(AA,[tz0])

  axA1.axis('equal')
  axA1.set_xlim([xl1,xl2])
  axA1.set_ylim([yl1,yl2])

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
      x0   = xFS
      y0   = yFS
      dltD = 100.
    else:
      x0   = TCNT[-1,0]
      y0   = TCNT[-1,1]
      dltD = 10.

    xsgm, ysgm, imin, jmin = arange_1segm(CNTR,x0,y0, dltD=dltD)  

    if len(xsgm) == 0:
      continue

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

def arange_1segm(CNTR, x0, y0, dltD=50.):
  """
    Find segment closest to x0, y0 
    arange the orientation of the segment
    so that it starts from the pnt closest to x0, y0
    CNTR - list with segments X,Y as np arrays
    Ignore contours that are > dltD points from the previous segment
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

# Disconnected segment - ignore:
  if np.min(DFS > dltD):
    xsgm = []
    ysgm = []
    imin = []
    jmin = []

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
  navonmb0 = rnmb - 1 # NAVO file naming: day X keeps f/cast of GW for X+1
  rdnavo0  = mtime.dnumb2rdate(navonmb0, ihours=False)

  print(f"NAVO GS: requested date {rdate}")

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
    rfcst   = mtime.dnumb2rdate(navonmb+1, ihours=False)
    pthdat  = pthnavo + rdnavo + '/wtxtbul/'
    if not os.path.isdir(pthdat):
      print('NAVO dir is missing ' + pthdat)
      continue
# Find file:
# Look for gs*.sub file with Gulf Stream front:
    for ffl in os.listdir(pthdat):
      if ffl.startswith('gs') and ffl.endswith('sub'):
        fpnavo = pthdat + ffl
        print(f"{rdnavo}:  ok --> f/cast for {rfcst}")
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
    

def derive_GSssh(fina, finb, HH, II, JJ, ISH, JSH, SSH=[],\
                 xlim1=2520, xlim2=3020, ylim1=1760, ylim2=2250, sshg=0.05):
  """
    Derive GS contour using SSH fields
    for 0.08 HYCOM fields: RTOFS, GOFS3.1
    xlim1/xlim2 ylim1/ylim2 - defines GS region 
    sshg - GS contour
  """
  import mod_misc1 as mmisc
  import mod_read_hycom as mhycom

  print(f'Deriving GS contour from {fina}')

  huge = 1.e20
  rg   = 9806.

  if len(SSH) == 0:
    SSH, _, _, _ = mhycom.read_hycom(fina,finb,'srfhgt')
    SSH[SSH>huge] = np.nan
    SSH = SSH/9.806

    ii1 = min(II)
    ii2 = max(II)
    jj1 = min(JJ)
    jj2 = max(JJ)
    SSHmn = np.nanmean(SSH[jj1:jj2,ii1:ii2])

    f_dmean = True
    if f_dmean:
      SSH = SSH-SSHmn

# Find Gulfstream
  IDM = SSH.shape[1]
  JDM = SSH.shape[0]
  X, Y = np.meshgrid(np.arange(IDM), np.arange(JDM))
  MS, IM, JM    = mmisc.inpolygon_v2(X,Y,II,JJ)  # Gulfstream region
  Msh, Ish, Jsh = mmisc.inpolygon_v2(X,Y,ISH,JSH)  # shelves to exclude

  z0 = -100.
  SSHsub = SSH.copy()
#SSHsub[np.where( (HH<0) & (HH>=z0) & (LAT>24.2) )] = -2.
  SSHsub[np.where(MS==0)] = np.nan
  SSHsub[np.where((Msh==1) & (HH>z0))] = np.nan
# Patches near Bahamas
  SSHsub[1808:1815, 2571:2586] = np.nan

  GSCNT = derive_contour(SSHsub, tz0=sshg)
  IGS   = GSCNT[:,0]
  JGS   = GSCNT[:,1]

  return IGS, JGS

def derive_GSssh_anls(pthanls, fanls, IIr, JJr, ISHr, JSHr, \
                      LONr, LATr, HH, sshg=0.05, remap=True):
  """
    Derive GS contour from SSH 2dvar analysis field
    Note: different from RTOFS grid !
    remap - maps onto RTOFS grid
  """
  import mod_misc1 as mmisc
  import mod_read_ncoda as mncoda
  import mod_valid_utils as mvalid
  import mod_read_hycom as mhycom
  importlib.reload(mvalid)
  import mod_rtofs as mrtofs

  print(f'Deriving GS contour from {fanls}')

  fdlon  = mvalid.find_fileindir(pthanls, 'grdlon_sfc_1o3241x2441_*datafld')
  fdlat  = mvalid.find_fileindir(pthanls, 'grdlat_sfc_1o3241x2441_*datafld')
  fdtopo = mvalid.find_fileindir(pthanls, 'depths_sfc_1o3241x2441_*datafld')

  jgrd  = 2441
  igrd  = 3241
  LON   = mncoda.read_2Dbin(fdlon, igrd, jgrd)
  LAT   = mncoda.read_2Dbin(fdlat, igrd, jgrd)
  HH    = mncoda.read_2Dbin(fdtopo, igrd, jgrd)
  if np.min(HH) > -1.e-30:
    HH = -HH

  HH    = np.where(HH>=-1.e-30, 100., HH)
  IDM   = HH.shape[1]
  JDM   = HH.shape[0]

  LON = np.where(LON > 180, LON-360, LON)

# Read SSH 2D analysis:
  fdanls = os.path.join(pthanls,fanls)
  SSH    = mncoda.read_2Dbin(fdanls, igrd, jgrd)
  SSH    = np.where(HH >= 0., np.nan, SSH)

  # Match indicex RTOFS ---> anls grid
  II, JJ = mmisc.match_indices(IIr,JJr,LONr,LATr,LON,LAT)

  ii1 = min(II)
  ii2 = max(II)
  jj1 = min(JJ)
  jj2 = max(JJ)
  SSHmn = np.nanmean(SSH[jj1:jj2,ii1:ii2])

  f_dmean = True
  if f_dmean:
    SSH = SSH-SSHmn

  # Gulfstream region:
  xR1 = 2520
  xR2 = 3020
  yR1 = 1760
  yR2 = 2250
  xy1 = mmisc.match_indices([xR1],[yR1],LONr,LATr,LON,LAT)
  xy2 = mmisc.match_indices([xR2],[yR2],LONr,LATr,LON,LAT)
  xlim1 = xy1[0][0]
  ylim1 = xy1[1][0]
  xlim2 = xy2[0][0]
  ylim2 = xy2[1][0]

  # Find Gulfstream
  ISH, JSH = mmisc.match_indices(ISHr,JSHr,LONr,LATr,LON,LAT)
  X, Y = np.meshgrid(np.arange(IDM), np.arange(JDM))
  MS, IM, JM    = mmisc.inpolygon_v2(X,Y,II,JJ)  # Gulfstream region
  Msh, Ish, Jsh = mmisc.inpolygon_v2(X,Y,ISH,JSH)  # shelves to exclude

  z0 = -100.
  SSHsub = SSH.copy()
  #SSHsub[np.where( (HH<0) & (HH>=z0) & (LAT>24.2) )] = -2.
  SSHsub[np.where(MS==0)] = np.nan
  SSHsub[np.where((Msh==1) & (HH>z0))] = np.nan
  # exclude shallow bahamas
  SSHsub[1437:1458,2527:2555] = np.nan

  # control point, Fl Stratis:
  xFSr = 2565
  yFSr = 1809
  #xyFS = mmisc.match_indices([2585],[1809],LONr,LATr,LON,LAT)
  #xFS = xyFS[0][0]
  #yFS = xyFS[1][0]
  xFS = 2518
  yFS = 1445

  GSCNT = derive_contour(SSHsub, tz0=sshg, xFS=xFS, yFS=yFS, \
                               xl1=xlim1, yl1=ylim1, xl2=xlim2, yl2=ylim2)
  IGS = GSCNT[:,0]
  JGS = GSCNT[:,1]

  # map anls indices onto HYCOM grid:
  if remap:
    IGS0 = IGS.copy()
    JGS0 = JGS.copy()

    print('Remapping analysis GS indices onto 0.08 grid')
    nigs = len(IGS)

    for ii in range(nigs):
      if ii%100 == 0:
        print(f'   {float(ii/nigs)*100:.1f}% done ...')

      ix0 = IGS0[ii]
      jx0 = JGS0[ii]      
      lon0, lat0 = mrtofs.interp_indx2lonlat(ix0, jx0, LON, LAT)

      ix0_hycom, jx0_hycom = mrtofs.interp_lonlat2indx(lon0,lat0,LONr,LATr)

      IGS[ii] = ix0_hycom
      JGS[ii] = jx0_hycom

  return IGS, JGS

