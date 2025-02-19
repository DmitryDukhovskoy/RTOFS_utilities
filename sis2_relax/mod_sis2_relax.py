"""
  Subroutine for preparing SIS2 relaxation files
"""
import xarray
import os
import importlib
import numpy as np
import sys
import matplotlib.pyplot as plt
from yaml import safe_load

PPTHN = '/home/Dmitry.Dukhovskoy/python'
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')
sys.path.append('./seasonal-workflow')

import mod_misc1 as mmisc1
import mod_mom6 as mom6util
importlib.reload(mom6util)
from mod_utils_fig import bottom_text


def interp2Dfld(A2d, IMOM, JMOM, INDX, JNDX, LMsk, LONs, LATs, hlon, hlat, \
                eps_err=1.e-2, land_mask=False):
  """
    Interpolate A2d (2D field) from PIOMAS onto MOM6 grid
    IMOM, JMOM - MOM6 indices where fields need to be interpolated
    INDX, JNDX - n x 4 arrays of PIOMAS grid points (gmapi) for bilinear interpolation
    LMsk - land/ocean mask of MOM6 grid
    LONs, LATs - grid PIOMAS
    hlon, hlat - grid MOM6
    land_mask - False: do not mask land with nans, keep filled with ocean values
  """
  import mod_utils_ob as muob
  import mod_bilinear as mblnr
  importlib.reload(mblnr)

# Find basis functions for a reference rectangle:
  phi1,phi2,phi3,phi4 = mblnr.basisFn_RectRef()
  phi_basis           = np.array([phi1, phi2, phi3, phi4]).transpose() # basis funs in columns

  npnts = len(INDX)
  jdm, idm = LMsk.shape

  print(f'Interpolating 2D {npnts} pnts ...')
# Make sure that gmapi is for the right section:
  assert npnts==len(IMOM), "INDX and IMOM mismatch in length"

  Ai = np.zeros((jdm,idm))
  if land_mask:
    Ai = np.where(LMsk==0,np.nan, Ai)

  for ikk in range(npnts):
    if ikk%7500 == 0:
      print(f'   {ikk/npnts*100.:.2f}% done ...')
    imom = IMOM[ikk]
    jmom = JMOM[ikk]
    x0   = hlon[jmom, imom]
    y0   = hlat[jmom, imom]
    II = np.squeeze(INDX[ikk,:])
    JJ = np.squeeze(JNDX[ikk,:])
    ii1, ii2, ii3, ii4 = II
    jj1, jj2, jj3, jj4 = JJ

    xx = LONs[JJ,II]
    yy = LATs[JJ,II]

# Use cartesian coordinates for mapping
    if x0 < 0.: x0 = x0+360.
    xx = np.where(xx<0., xx+360., xx)
    f_repeated= muob.check_repeated_vertices(xx,yy)
    if f_repeated:
      print(f"Bad box with coninciding vertices ikk={ikk}, approximate interpolation")
      xht = 1.e-3
      yht = 1.e-3
    else:
      xref     = x0-0.1
      yref     = y0-0.1
      XV, YV   = mblnr.lonlat2xy_wrtX0(xx, yy, xref, yref)
      x0c, y0c = mblnr.lonlat2xy_pnt(x0, y0, xref, yref)
      xht, yht = mblnr.map_x2xhat(XV, YV, x0c, y0c)   # cartesian coord

    # Fix round off errors for points on the side of the ref. square that are close to +/-1:
    if abs(xht)-1. < eps_err:
      xht = np.round(xht)
    if abs(yht)-1. < eps_err:
      yht = np.round(yht)

  # check if xht and yht <= 1
  # If the grid point inside the box, then mapping is the problem
  # if not - then these are a few cases near singularities of SPEAR I/J axes
  # over land where bounding boxes could not be located 
  # In a few cases, box coordinates may repeat due to SPEAR singular I/J grid
  # that does not change coordinates for I/J near singular regions (boxes are triangulars)
  # THis results in singular matrix A for mappring
  # For very thin  rotated, skewed quadrilaterals mapping does not work well
  # Try to rotate the quadrilateral
    if abs(xht) > 1 or abs(yht) > 1:
      xht0 = xht
      yht0 = yht
      XVr, YVr, x0r, y0r  = muob.rotate_box(XV, YV, x0c, y0c)
      xht, yht = mblnr.map_x2xhat(XVr, YVr, x0r, y0r)
#      print(f"ERR: ikk={ikk} Mapping ref box xht={xht0:6.3f} yht={yht0:6.3f}, fixed " +\
#            f" xht={xht:6.3f}, yht={yht:6.3f}")

    if abs(xht) > 1. or abs(yht) > 1.:
# If nothing works, interpolate into the center
# these should be very rare for locations on land where I / J axes converge 
      print(f"Fixing by rotating ref BOX failed ikk={ikk} " +\
            f"xht={xht:8.5f} yht={yht:8.5f}, approximate xhy, yht as middle pnt")
      xht = 1.e-3
      yht = 1.e-3

    aa1 = np.squeeze(A2d[jj1,ii1])
    aa2 = np.squeeze(A2d[jj2,ii2])
    aa3 = np.squeeze(A2d[jj3,ii3])
    aa4 = np.squeeze(A2d[jj4,ii4])

    HT  = np.array([aa1, aa2, aa3, aa4]).transpose()
    # Typically, land values should be filled
    # in case, they have not:  Get rid off nans
    nnans = len(np.where(np.isnan(HT))[0])
    if nnans == len(HT):
      Ai[jmom,imom] = np.nan
      continue
    else:
      mnv = np.nanmean(HT)
      HT = np.where(np.isnan(HT), mnv, HT)

    hintp  = mblnr.bilin_interp(phi1, phi2, phi3, phi4, xht, yht, HT)

# Check: abs. values of interpolated values <= original data
    mxHT = np.max(abs(HT))
    if mxHT == 0: 
      mxHT = 1.e-20
    dmx  = abs(hintp)/mxHT
    if dmx > 1.1:
      print(f"!!! segm{nsgm} {varnm} Min/Max test violated: ikk={ikk} dlt: {dmx}")
      if dmx > 1.5:
        raise Exception("MinMax test: Interp error Check interpolation")

    Ai[jmom,imom] = hintp

  return Ai


def read_PIOMAS(yr0, mm0, dfpiomas, varnm):
  """
  Derive thikness or conc. fields for yr0, mm0 
  dfpiomas = dir + filename

  monthly fields
  1901 - 2010
  https://psc.apl.uw.edu/research/projects/piomas-20c/

  PIOMAS-20C is a sea ice thickness reconstruction covering the period 1901-2010. 
  It is constructed using a coupled ice-ocean model using atmospheric forcing data from 
  the ECMWF ERA-20C reanalysis to provide atmospheric forcing. Sea ice concentrations 
  from the Hadley Center HadISST v2.0 data set are assimilated to constrain the model 
  at the ice-edge. 

  
  """
  import mod_time as mtime

  ds_piomas = xarray.open_dataset(dfpiomas)
  varconc = 'sic'
  varthck = 'sit'

  print(f'Reading PIOMAS {yr0}/{mm0} {dfpiomas}')

  if not varnm=='sic' and not varnm=='sit':
    raise Exception (f'PIOMAS variables are sic and sit, requested {varnm}')


  # Find record #: days since 1901-01-01 = day=1, index=0
  dnmb0 = mtime.datenum([yr0,mm0,1])
  dnmbR = mtime.datenum([1901,1,1])
  ndays = int(dnmb0-dnmbR) + 1
  #Time  = dset['time'].data  # np datetime array
  Month = ds_piomas['month'].data
  Year  = ds_piomas['year'].data
  D     = np.sqrt((Month-mm0)**2 + (Year-yr0)**2)
  tindx = np.argmin(D)
  A2d   = ds_piomas[varnm].data[tindx,:].squeeze()  # thikness, m

  return A2d

def linear_distr1D(A1d, nav=3):
  """
    Spread out a value over adjacent cells
    using linear distribution function (averaging)
    input: 1D array
  """
  Afltr = A1d.copy()
  npnts = len(A1d)
  nav_hlf = int(np.floor(nav/2))
  assert nav<npnts, f"Number of averaged grid points {nav} should be < {npnts}"
  for ik in range(npnts):
    i1 = ik-nav_hlf
    i2 = ik+nav_hlf
    i1 = max([0,i1])
    i2 = min([i2,npnts-1])+1
    a_mn = np.mean(A1d[i1:i2])
    Afltr[ik] = a_mn

  return Afltr

def gauss_distr1D(A1d, icat0, sgm=1.3, conserve=True):
  """
    Spread out a value over adjacent cells
    using gaussian distribution function (averaging)
    input: 1D array
    sgm - controlls the spread of the gaussian filter
    icat0 - ice category where the initial ice is given
  """
  Afltr = A1d.copy()
  npnts = len(A1d)
  aa0 = A1d[icat0]
  XX = np.arange(npnts)
  Afltr = aa0*np.exp(-(XX-icat0)**2/(sgm**2))

  if conserve:
  # Conserve grid ice concentration:
    Afltr = aa0*Afltr/np.sum(Afltr)

  return Afltr

def adjust_hice(ICAT, hcat, ccat, ncat, hice, eps0=1.e-6):
  """
  # Adjust hice for all cats except for the last one
  ICAT - ice cat min/max thicknesses 
  """
  chcat = hcat*ccat
  for kk in range(ncat-1):
    hmin = ICAT[kk]
    hmax = ICAT[kk+1] - eps0
    hcat_k = chcat[kk]/ccat[kk]
    if hcat_k < hmin:
    # increase hcat and decrease ccat
      hcat_k = hmin + eps0
    elif hcat_k > hmax:
    # decrease hcat and increase ccat
      hcat_k = hmax - eps0
    hcat[kk] = hcat_k

    ctot = np.sum(ccat)
    htot = np.sum(ccat*hcat)
    print(f"ice cat={kk+1} ctot={ctot:.2f} htot={htot:.2f}")

  chcat = hcat*ccat
  htot = np.sum(chcat)
  ctot = np.sum(ccat)
  dlt_chcat = htot - hice  # want dlt_chcat = 0.

  # Adjust last cat to keep hice and cice conserved:
  if abs(dlt_chcat) > eps0:
    dlt_h = dlt_chcat/ccat[ncat-1]
    hcat[ncat-1] = hcat[ncat-1] - dlt_h
  # Check ice cat limits:
  dlt_hice = np.zeros((ncat))
  for kk in range (ncat):
    hmin = ICAT[kk]
    hmax = ICAT[kk+1] - eps0
    if hcat[kk] < hmin:
      dlt_hice[kk] = hcat[kk]-hmin
    elif hcat[kk] > hmax:
      dlt_hice[kk] = hcat[kk]-hmax
    print(f"Up swap: icat={kk+1} h[k] exceeds ice limits by = {dlt_hice[kk]:.6f}")
    print(f"icat={kk+1} hmin={hmin:.3f}/hmax={hmax:.3f} hcat={hcat[kk]:.3f}")

  return hcat, ccat

def redistribute_hice(hice, cice, ICAT0=[], eps0 = 1.e-6, itd_method='equal'):
  """
    Redistribute ice thickness (grid cell mean), 1 pnt, 
    by ice categories such that the mean grid 
    ice thickness is conserved, i.e.
    htot = sum(i_cat)(h_cat(i)*conc_cat(i)) = hice
    ctot = sum(i)(conc_cat(i)) = cice
    eps0 - close to 0 value used to check errors and "0"

    ice thickn. distribution methods for initial thkn and conc distr. by cats:
    - simple = place all ice in the category that corresponds hice (grid cell mean value)
    - gauss = start with ice concentration aplying Gaussian filter on cice in the cat. = hice
    - equal = start with ice concentration equally distributed over the cats.

  """
  if len(ICAT0) == 0:
    ICAT0 = np.array([1.0e-10, 0.1, 0.3, 0.7, 1.1])
  ncat  = len(ICAT0)

  err_lim = 1.e-3  # allow higher error for ITD adjustment to avoid iteration errors

  # Add upper bound for the last cat:
  hi_max = 75.
  ICAT = np.append(ICAT0,[hi_max])

  hcat = np.zeros((ncat))
  ccat = np.zeros((ncat))

  #itmp = ICAT.copy()*0.
  if hice < ICAT[0]:
  # open water
    return hcat, ccat
  if cice < eps0: 
  # open water
    return hcat, ccat

# Find ice cat. where grid cell mean hice falls in:
  indx_cat = 999
  for kk in range(ncat):
    hbnd = ICAT[kk]
    if kk < ncat-1:
      hbnd_up = ICAT[kk+1]
    else:
      hbnd_up = 100.

    if hice >= hbnd and hice < hbnd_up:
      indx_cat = kk
      break    

  assert indx_cat < ncat, f"Could not find ice cat for {hice}"
  print(f"hice={hice:.2f} m, assigned cat={indx_cat+1} hlim={ICAT[indx_cat]:.2f}/{ICAT[indx_cat+1]:.2f}")
  hcat = np.zeros((ncat))
  ccat = np.zeros((ncat))
  ccat[indx_cat] = cice
  #chice = hice*cice      # m3/m2 or [m] - vol/m2 used in CICE chice*rho_ice = kg/m2
  # chice_cat array = hcat(k)*ccat(k) and sum(chice_cat)=sum(hcat(k)*ccat(k)) = hice !
  # so, when only 1 cat. is non-zero: chice_cat[k]=hice 
  chice_cat = np.zeros((ncat))
  chice_cat[indx_cat] = hice  # m3/m2 or [m] - vol/m2 used in CICE chice*rho_ice = kg/m2
 
  if hice/cice > np.max(ICAT):
  # ice will be too thick if redistribute by cats, try simple ITD:
    print(f"hice {hice:.3f} and cice {cice:.3f} will result in ice > {ICAT[-1]:.3f}")
    print(f"Simple ITD")
    hcat[indx_cat] = hice/cice
    return  hcat, ccat

  match itd_method:
    case('simple'):
      chcat = chice_cat.copy()
      #hcat[indx_cat] = hice/cice
      ccat = np.where(ccat<eps0, eps0, ccat)
    case('gauss'):
      chcat = gauss_distr1D(chice_cat, indx_cat, sgm=1.8) 
      ccat = gauss_distr1D(ccat, indx_cat, sgm=2.2)
    case('equal'):
      ccat = np.zeros((ncat)) + cice/ncat
    #case('weighted'):
    #  ccat = weighted_distr1D(chice_cat, indx_cat)

  # no 0 in the ccat:
  ccat = np.where(ccat<eps0, eps0, ccat)

  # First guess of h ice by categories given (hcat(k)*ccat(k)) and sum(ccat(k))=cice
  # Make sure sum(ccat)=cice
  cfctr = np.sum(ccat)/cice  # reduction factor to adjust concentration:
  ccat = ccat/cfctr
  hcat = chcat/ccat  # hice[k] may not be within ice thkn cat. limits, adjust

  # adjust hcat to keep thkn values within the ice categories:
  # readjust ccat to keep chcat[k] unchanged
  ccat_min = 1.e-10  # min concentration to keep in each category

  # Adjust hice for all cats except for the last one
  hcat, ccat = adjust_hice(ICAT, hcat, ccat, ncat, hice)
   
  chcat = hcat*ccat
  htot = np.sum(chcat)
  ctot = np.sum(ccat)
  print(f"Up swap: ctot={ctot:.4f} htot={htot:.4f}")

  # Expected that tot ice and cice are conserved here:
  assert abs(htot-hice) < eps0, f"1. htot={htot} hice={hice} do not match"
  assert abs(np.sum(ccat)-cice) < eps0, f"1. ={htot} hice={hice} do not match"

  # Down swap only if thickiest ice needs adjustment
  # other cats should be good
  # if hcat is too low - move from lower cats
  # if hcat is too high - increase conc, reduce hcat and adjust conc everywhere to keep cice
  hmin = ICAT[ncat-1]
  hmax = ICAT[ncat] - eps0
  hcat_k = hcat[ncat-1]
  adj_hthk = hcat_k < hmin or hcat_k > hmax
  icc = 0
  niter = 50
  while adj_hthk:
    icc += 1
    if icc > niter:
      break
   
    hk_old = hcat[ncat-1]
    print(f"iter={icc} Adjusting thick cat {ncat}: hcat[k]={hcat_k:.3f}")
    if hcat_k < hmin:
      # increase hcat and decrease ccat
      hcat_k = hmin + eps0
    elif hcat_k > hmax:
      # decrease hcat and increase ccat
      hcat_k = hmax - eps0
    hcat[ncat-1] = hcat_k
    #dlt_ccat = chcat[ncat-1]/hcat_k - ccat[ncat-1]
    ccat[ncat-1] = chcat[ncat-1]/hcat_k

    # Adjust ccat to conserve ctot damping or taking ice from the thickest cats
    ctot = np.sum(ccat)
    dlt_cice = ctot-cice
    wght = ccat/np.sum(ccat)
    ccat = ccat - wght*dlt_cice
    ccat = np.where(ccat<=eps0, eps0, ccat)
    hcat = chcat/ccat

    hcat, ccat = adjust_hice(ICAT, hcat, ccat, ncat, hice)

    hcat_k = hcat[ncat-1]
    dlt_hk = hcat_k - hk_old
    adj_hthk = (hcat_k > hmax)

    print(f"new h({ncat})={hcat_k:.3f}, dlt = {dlt_hk:.6f}")
    if abs(dlt_hk)<1.e-6:
      break
    


#   icc=0
#   while abs(dlt_cice) > err_lim:
#     icc += 1
#     if icc > 20:
#       break
#     for kk in range(0,ncat-1):
#       ck_old = ccat[kk]
#       hk_old = hcat[kk]
#       chk_old = ck_old*hk_old
#       cmax = cice - (np.sum(ccat)-ccat[kk]) # max conc for this cat
#       ck_new = ccat[kk] - dlt_cice
#       ck_new = np.max([eps0,ck_new])
#       ck_new = np.min([cmax, ck_new])
#       hk_new = chcat[kk]/ck_new
#       # make sure new hice[k] is withim the limits for the category:
#       hmin = ICAT[kk]
#       hmax = ICAT[kk+1] - eps0
#       hk_new = np.max([hmin,hk_new])
#       hk_new = np.min([hmax,hk_new])
#       hcat[kk] = hk_new
#       ccat[kk] = ck_new
#       htot = np.sum(ccat*hcat) # should be conserved by the algorithm above
#       if abs(htot-hice) > eps0:
#       # readjust ice conc:
#         ck_new = chcat[kk]/hk_new
#         cmax = cice - (np.sum(ccat)-ccat[kk]) # max conc for this cat
#         ck_new = np.max([eps0,ck_new])
#         ck_new = np.min([1.-eps0, ck_new])
#         ccat[kk] = ck_new
# #
#       ctot = np.sum(ccat)
#       dlt_cice = ctot-cice
#       print(f"iter {icc} icat={kk+1} dlt_cice={dlt_cice:.6f} ctot={ctot:.3f}")
#       if abs(dlt_cice) <= err_lim:
#         break

  # Check final distribution: 
  htot = np.sum(ccat*hcat) # should be conserved 
  ctot = np.sum(ccat)
 
  print(f"Input:          hice={hice:.3f}, cice={cice:.3f}")
  print(f"After redistr:  htot={htot:.3f}, ctot={ctot:.3f}")
#  assert abs(htot-hice) < err_lim, f"3. htot={htot} hice={hice} do not match"
#  assert abs(ctot-cice) < err_lim, f"3. ctot={ctot} cice={cice} do not match"

  return hcat, ccat
  



