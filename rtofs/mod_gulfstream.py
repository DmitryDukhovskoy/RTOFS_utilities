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
    /scratch2/NCEPDEV/ovp/Lichuan.Chen/DCOM/$YYYYMMDD/wtxtbul.  
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
    print('NAVO gs*.sub file not found for {1} +/- {0}'.format(missing,rdnavo))
    XNW = []
    YNW = []
    fcstnmb = -999. 
    return XNW, YNW, fcstnmb
#    raise Exception(" NAVO not found, try increase missing search range")
      
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
    
def plot_gulfstr(TCNT, A3d, INW_navo, JNW_navo, LON, LAT, \
                 INW_rtofs, JNW_rtofs, sinfo='', btx='', ssinf=''):
  """
    Plot Gulf Stream map
    for checking front
  """
  print('=======   Start Plotting   =========')
  from matplotlib import cm
  import mod_rtofs as mrtofs
  import mod_time as mtime
  import mod_utils_fig as mufig

  #clrs   = cm.get_cmap('viridis',200)
  clrs   = cm.get_cmap('rainbow',200)
  clrs.set_bad(color=[0.3,0.3,0.3])
  #  clrs.set_under(color=[0.8,0.8,0.8])
  #  clrs.set_under(color=[0.8,0.7,1])

  plt.ion()
    
  rmin = 0.
  rmax = 28.
  tz0  = 12.0
    
  Xt = TCNT[:,0]
  Yt = TCNT[:,1]

  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()

  ax1 = plt.axes([0.1, 0.2, 0.7, 0.7])

  Tplt = np.squeeze(A3d[0,:,:])
  im1 = ax1.pcolormesh(Tplt, shading='flat', \
                 cmap=clrs,\
                 vmin=rmin, \
                 vmax=rmax)
  #  im1.set_clim(rmin,rmax)
  #  ax1.contour(HH, [0.0], colors=[(0,0,0)], linewidths=1)
  #  ax1.contour(T12, [tz0], colors=[(0.,0.4,0.6)], linewidths=1)
  clr_rtofs = [0,0.8,1]
  clr_navo  = [0,0,0]

  # Plot Gulf Stream front:
  ln1, = ax1.plot(Xt,Yt,'-',color=clr_rtofs, label="RTOFS")
  # NAVO Gulf Stream North Wall:
  ln2, = ax1.plot(INW_navo, JNW_navo,'-',color=clr_navo, label="NAVO")

  lons = np.linspace(-180,180,73)
  lats = np.linspace(-90,90,37)
  ax1.contour(LON, lons, linestyles='dotted', colors=[(0.8,0.8,0.8)])
  ax1.contour(LAT, lats, linestyles='dotted', colors=[(0.8,0.8,0.8)])

# Plot RTOFS front used for MHD:
  ax1.plot(INW_rtofs,JNW_rtofs,'-', color=[0,0,1])
  #
  xlim1 = 2520
  xlim2 = 3020
  ylim1 = 1760
  ylim2 = 2250 
  #  ax1.axis('equal')
  ax1.axis('scaled')
  ax1.set_xlim([xlim1,xlim2])
  ax1.set_ylim([ylim1,ylim2])
  ax1.set_xticklabels([])
  ax1.set_yticklabels([])
  ax1.set_xticks([])
  ax1.set_yticks([])

  # Actual limits:
  Ylim = ax1.get_ylim()
  Xlim = ax1.get_xlim()
  Yl1  = int(Ylim[0])
  Yl2  = int(Ylim[1])
  Xl1  = int(Xlim[0])
  Xl2  = int(Xlim[1])
  # Put lon/lats on axis:
  lon1 = LON[Yl1,Xl1]
  lon2 = LON[Yl2,Xl2]
  lat1 = LAT[Yl1,Xl1]
  lat2 = LAT[Yl2,Xl2]

  iLN1 = np.min(np.where(lons>=lon1)[0])
  iLN2 = np.max(np.where(lons<=lon2)[0])
  iLT1 = np.min(np.where(lats>=lat1)[0])
  iLT2 = np.max(np.where(lats<=lat2)[0])
  dltx = 1
  dlty = 1
  if iLN2-iLN1 >= 8:
    dltx = 2
  if iLT2-iLT1 >= 8:
    dlty = 2

  # X-axis labels
  for ikk in range(iLN1,iLN2+1,dltx):
    xx0 = lons[ikk]
    yy0 = lat1   # for Mercator part of the grid, lat = const along j=fixed
    ii0, jj0 = mrtofs.find_indx_lonlat(xx0, yy0, LON, LAT)
    jj0 = Yl1-20
    xstl = '{0:3.1f}W'.format(abs(xx0))
    ax1.text(ii0, jj0, xstl,
             fontsize=12,
             horizontalalignment='center')

  # Y-axis labels
  for ikk in range(iLT1,iLT2+1,dlty):
    yy0 = lats[ikk]
    xx0 = lon1
    ii0, jj0 = mrtofs.find_indx_lonlat(xx0, yy0, LON, LAT)
    ii0 = Xl1-10
    ystl = '{0:3.1f}N'.format(abs(yy0))
    if jj0 > Yl2:
      continue
    ax1.text(ii0, jj0, ystl,
             fontsize=12,
             verticalalignment='center',
             horizontalalignment='right')

  ax1.set_title(sinfo)

  ax1.set_xlim(Xlim)
  ax1.set_ylim(Ylim)

  # Select coordinate of the region of interest:
  #f_setrgn = False
  #if f_setrgn:
  ## Bind the button_press_event with the onclick() method
  #  fig1.canvas.mpl_connect('button_press_event', onclick)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='both')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  # Legend:
  if len(ssinf) > 0:
    ax5 = plt.axes([0.7, 0.03, 0.2, 0.13])
    lgd = plt.legend(handles=[ln1,ln2], loc='upper left')
    ax5.text(0,0.01,ssinf)
    ax5.axis('off')

  if len(btx) > 0:
    mufig.bottom_text(btx, pos=[0.1, 0.03])



