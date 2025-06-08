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
    if ikk%10000 == 0:
      print(f'   {ikk/npnts*100.:.2f}% done ...')
    imom = IMOM[ikk]
    jmom = JMOM[ikk]
    #print(f'ikk={ikk}') 
    if LMsk[jmom,imom] == 0:
      continue
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
      # For rotated grid boxes, mapping may give singular matrix AA
      # Try to rotate the quadrilateral to orient sides with X and Y axis:
      XVr, YVr, x0r, y0r  = muob.rotate_box(XV, YV, x0c, y0c)
      #xht, yht = mblnr.map_x2xhat(XV, YV, x0c, y0c)   # cartesian coord
      xht, yht = mblnr.map_x2xhat(XVr, YVr, x0r, y0r)

    # Fix round off errors for points on the side of the ref. square that are close to +/-1:
    if abs(xht)-1. < eps_err:
      xht = np.round(xht)
    if abs(yht)-1. < eps_err:
      yht = np.round(yht)

  # Already rotated now, no need in this check
  # check if xht and yht <= 1
  # If the grid point inside the box, then mapping is the problem
  # if not - then these are a few cases near singularities of SPEAR I/J axes
  # over land where bounding boxes could not be located 
  # In a few cases, box coordinates may repeat due to SPEAR singular I/J grid
  # that does not change coordinates for I/J near singular regions (boxes are triangulars)
  # THis results in singular matrix A for mappring
  # For very thin  rotated, skewed quadrilaterals mapping does not work well
  # Try to rotate the quadrilateral
  #  if abs(xht) > 1 or abs(yht) > 1:
  #    xht0 = xht
  #    yht0 = yht
  #    XVr, YVr, x0r, y0r  = muob.rotate_box(XV, YV, x0c, y0c)
  #    xht, yht = mblnr.map_x2xhat(XVr, YVr, x0r, y0r)
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

  print(f'Reading PIOMAS-reconstruct {yr0}/{mm0} {dfpiomas}')

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

def read_PIOMASv21(yr0, mm0, dfpiomas, varnm):
  """
  Derive thikness or conc. fields for yr0, mm0 
  dfpiomas = dir + filename

  monthly fields from PIOMAS v2.1 reanalysis
  1979-present

  monthly fields
  1979-present
  https://pscfiles.apl.washington.edu/zhang/PIOMAS/data/v2.1/

  PIOMASv2.1  is a sea ice reanalysis
  sea ice concentration (edge) is assimilated using sat. ice conc. 
  """
  import mod_time as mtime
  varthck = 'heff'
  varconc = 'area'

  print(f'Reading PIOMASv2.1 {yr0}/{mm0} {dfpiomas}')

  if not varnm==varconc and not varnm==varthck:
    raise Exception (f'PIOMAS variables are {varcon} and {varthck}, requested {varnm}')

  ds_piomas = xarray.open_dataset(dfpiomas)
  Month = ds_piomas['month'].data
  Year  = ds_piomas['year'].data
  D     = np.sqrt((Month-mm0)**2 + (Year-yr0)**2)
  tindx = np.argmin(D)
  assert(D[tindx]==0), f"Requested {yr0}/{mm0} not found in {dflthkn}"
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

def adjust_hcat(ICAT, hcat, ccat, ncat, hice, eps0=1.e-6, verbose=False):
  """
    Try to Adjust ice thkns in all cats
    Dumping all extra ice into the last cat
    This algorithm works mostly when hcat[k] > hlim[k]
    ICAT - ice cat min/max thicknesses 
    Last category - max ice thickness is not bounded
  """
  chcat = hcat*ccat
  for kk in range(ncat):
    hmin = ICAT[kk]
    if kk<ncat-1:
      hmax = ICAT[kk+1] - eps0
    else:
      hmax = 1.e3

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
    if verbose:
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
  if verbose:
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

def adjust_high_thkn(hcat, ccat, ICAT, ncat, hice, eps0, verbose=False):
  """
    Reduce ice in thikest cat
  """
  chcat = hcat*ccat  # preserve ice mass
  hmin = ICAT[ncat-1]
  hmax = ICAT[ncat] - eps0
  hcat_k = hcat[ncat-1]
  adj_high = hcat_k > hmax # adjust excesss ice in thick. cat.

  if not adj_high:
    print(f"no high thkn found, H({ncat})={hcat_k:.3f}")

  icc = 0
  niter = 50
  while adj_high:
    icc += 1
    if icc > niter:
      break
   
    hk_old = hcat[ncat-1]
    print(f"adj_high: iter={icc} Adjusting thick cat {ncat}: hcat[k]={hcat_k:.3f}")
    if hcat_k > hmax:
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

    hcat, ccat = adjust_hcat(ICAT, hcat, ccat, ncat, hice, verbose=False)

    hcat_k = hcat[ncat-1]
    dlt_hk = hcat_k - hk_old
    adj_high = hcat_k > hmax

    print(f"adj_high: iter={icc} new h({ncat})={hcat_k:.3f}, dlt = {dlt_hk:.6f}")
    if abs(dlt_hk)<1.e-6:
      break
    
  return hcat, ccat

def adjust_low_thkn(hcat, ccat, ICAT, ncat, hice, eps0, verbose=False):
  """
    Increase ice in thikest cat
    redistribute ice from lower cats
  """
  chcat = hcat*ccat  # preserve ice mass
  ithk = ncat-1
  hmin_thk = ICAT[ithk] + eps0
  hmax_thk = ICAT[ithk+1] - eps0
  hcat_thk = hcat[ithk]
  adj_low = hcat_thk < hmin_thk # adjust excesss ice in thick. cat.

  if not adj_low:
    print(f"no low thkn found, H({ncat})={hcat_thk:.3f}")

  # Estimate how much ice needs to be added for min thkn and conc:
  cmax_thk  = hice*cice/hmin_thk    # max concentration if all ice keep in thkst cat
  cmin_thk  = 1.e-3*cmax_thk        # some small conc >> eps0
  chmin_thk = hmin_thk*cmin_thk     # min ice mass in thkst cat to have required conc and thkn
  dch_thk   = chmin_thk-chcat[ithk] # change of ice mass in thkst cat

  for kk in range(ithk): 
    # go from 1st to the next to last cat:
    hmin = ICAT[kk] + eps0
    hmax = ICAT[kk+1] - eps0
    cmin = eps0
    hk_old = hcat[kk]
    ck_old = ccat[kk]
    chk_old = chcat[kk] 
    # Check if there is enough ice to adjust the thick. cat
    # Then find how much ice needs to be redistributed from this cat
    if (chk_old-hmin*cmin) >= dlt_ch:
      chk_new = chk_old-dlt_ch
      wgt = chk_new/chk_old
      hk_new = hk_old*wgt
      hk_new = np.max([hmin, hk_new])
      ccat[kk] = chk_new/hk_new
      hcat[kk] = hk_new 
    else:
    # Insufficient ice in this cat, take all leaving min ice mass 
      hcat[kk] = hmin
      ccat[kk] = cmin
      dlt_ch_k = chk_old - hcat[kk]*ccat[kk]
      assert dlt_ch_k > 0, f"dlt_chk should be >0 dlt_chk={dlt_ch_k:.4f}"
      dlt_ch = dlt_ch - dlt_ch_k

    # Adjust thickest cat:
    ccat[ithk]  = cice - np.sum(ccat[:ithk])
    chcat[ithk] = chcat[ithk] + dlt_ch
    hcat[ithk]  = chcat[ithk]/ccat[ithk]
    dlt_ch = ch_min-chcat[ithk]
    if dlt_ch <= eps0:
      break

  return hcat, ccat

def check_hcice(hcat, ccat, hice, cice, eps0=1.e-6, verb=False):
  """
    CHeck cice and hice conserved
  """
  chcat = hcat*ccat
  err_cice = False
  err_hice = False
  if abs(np.sum(chcat)-hice) > eps0:
    err_hice = True
    if verb:
      print(f"hice {hice:.3f} not conserved: {np.sum(chcat):.3f}")

  if abs(np.sum(ccat)-cice) > eps0:
    err_cice = True
    if verb:
      print(f"cice {cice:.3f} not conserved: {np.sum(cice):.3f}")

  return err_hice, err_cice

def check_hcat(hcat, ccat, ICAT, ncat, eps0, icat=-1):
  """
    Check if hcat[k] is within the cat thkn limits
    specify icat to check 1 cat
    returns adj_hcat:
    <0 - below the hmin
    >0 - exceeds hmax
    0 - within the limits
  """
  adj_hcat = 0
  if icat<-1:
    ii1=0
    ii2=ncat
  else:
    ii1=icat
    ii2=ii1+1

  for kk in range(ii0,ii1):
    if ccat[kk] < eps0:
    # No ice
      continue
    hmin = ICAT[kk]
    if kk<ncat-1:
      hmax = ICAT[kk+1] - eps0
    else:
      hmax = 1.e3

    if hcat[kk] < hmin:
      adj_hcat = hcat[kk]-hmin
    elif hcat[kk] > hmax:
      adj_hcat = hcat[kk]-hmax

  return adj_hcat 

def redistribute_hice(hice, cice, ICAT=[], eps0=1.e-6, ck_min=1.e-3, verbose=False):
  """
    Redistribute ice by ice cats
    such that the mean grid 
    ice thickness is conserved, i.e.
    htot = sum(i_cat)(h_cat(i)*conc_cat(i)) = hice
    ctot = sum(i)(conc_cat(i)) = cice
    eps0 - close to 0 value used to check errors and "0"
    ck_min = 1.e-3  ice conc. in lower cats,  some small value >> eps0

    ice thkn in the thikest cat may not match the limits for this cat !

    Start with dumping all ice into the thickest cat 
    such that hmin(k)<hice < hmax(k)
    and then redistribute 
    into lower cats
  """
  if len(ICAT) == 0:
    ICAT = np.array([1.0e-10, 0.1, 0.3, 0.7, 1.1])
  ncat  = len(ICAT)
  ithk  = ncat-1

  err_lim = 1.e-3  # allow higher error for ITD adjustment to avoid iteration errors

  hcat = np.zeros((ncat))
  ccat = np.zeros((ncat))

  #itmp = ICAT.copy()*0.
  if hice < ICAT[0]:
  # open water
    return hcat, ccat
  if cice < eps0:
  # open water
    return hcat, ccat

  # Check that there is enough ice to be i
  # distributed over the cats for min cice and hice:
  #ck_min = 1.e-3  # some small value >> eps0
  htot_min = np.sum(ck_min*ICAT)
  if hice < htot_min or cice < ck_min*ncat:
    hcat = np.zeros((ncat))+eps0
    ccat = np.zeros((ncat))+eps0
    ccat[ithk] = cice - np.sum(ccat[:ithk]) 
    hcat[ithk] = hice/ccat[ithk]    # can be some high values due to low cice
    chcat = ccat*hcat
    # Limit max hcat  aand adjust ccat allowing to be not exact?

    return hcat, ccat

# Find ice cat. where grid cell mean hice falls in:
  icat = 1e3
  hcat_k = hice/cice     # ice thickness in a category
  for kk in range(ncat):
    hbnd = ICAT[kk]
    if kk < ncat-1:
      hbnd_up = ICAT[kk+1]
    else:
      hbnd_up = 1.e3

    if hcat_k >= hbnd and hcat_k < hbnd_up:
      icat = kk
      break
   
  assert icat < ncat, f"Could not find ice cat for {hice}"
  print(f"hice={hice:.2f}m --> cat={icat+1}:  {hbnd:.2f}/{hbnd_up:.2f}")
 
  hcat = np.zeros((ncat))
  ccat = np.zeros((ncat))
  ccat[icat]  = cice
  hcat[icat]  = hice/cice 
  chcat = ccat*hcat

  err_hice, err_cice  = check_hcice(hcat, ccat, hice, cice)
  assert not err_hice, f"1. error hice not conserved"
  assert not err_cice, f"1. error cice not conserved"
    
  for kk in range(icat):
    hmin = ICAT[kk]
    if kk<ncat-1:
      hmax = ICAT[kk+1] - eps0
    else:
      hmax = 1.e3

    ccat_k = ck_min
    hcat_k = ICAT[kk] + eps0
    dch_k = ccat_k*hcat_k
    if dch_k >= chcat[icat]:
      # cannot redistribute ice, not enough in thickest cat
      print(f"icat={kk+1} cannot redistribute ice, not enough in thickest cat")
      break

    ccat[kk] = ccat_k
    hcat[kk] = hcat_k
    chcat[kk] = dch_k

    # Update the cat with initial ice:
    cnew = ccat[icat]-ccat[kk]
    cnew = np.max([ck_min, cnew])
    ccat[icat] = cnew
    chcat[icat] = chcat[icat] - dch_k 
    hcat[icat] = chcat[icat]/ccat[icat]

    # Check if ice thkn is within the limits 
    #adj_hcat = check_hcat(hcat,ccat, ICAT, ncat, eps0, icat=icat)
 
  ctot = np.sum(ccat)
  htot = np.sum(ccat*hcat)
  print(f"redistr: ctot={ctot:.3f} cice={cice:.3f}, htot={htot:.3f} hice={hice:.3f}")
  err_hice, err_cice  = check_hcice(hcat, ccat, hice, cice)
  assert not err_hice, f"1. error hice not conserved"
  assert not err_cice, f"1. error cice not conserved"

  return hcat, ccat

def redistribute_hice_v0(hice, cice, ICAT0=[], eps0 = 1.e-6, \
    itd_method='equal', verbose=False):
  """
    The script does not work well for wierd hice/cice values, e.g.
    hice=3, cice=0.3 - thickest cat may have wrong hcat[k]
    or small hice: hice=0.1, cice=0.8 - has problems with distributing
    across all cats 

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
    htot = np.sum(hcat*ccat)
    assert abs(htot-hice) < eps0, f"Failed to keep hice Simple ITD: htot={htot:.6f}"
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
  #ccat_min = 1.e-10  # min concentration to keep in each category

  # Adjust hice for all cats except for the last one
  hcat, ccat = adjust_hcat(ICAT, hcat, ccat, ncat, hice, verbose=verbose)
   
  chcat = hcat*ccat
  htot = np.sum(chcat)
  ctot = np.sum(ccat)
  
  print(f"Up swap: ctot={ctot:.4f} htot={htot:.4f}")

  # Expected that tot ice and cice are conserved here:
  assert abs(htot-hice) < eps0, f"1. htot={htot} hice={hice} do not match"
  assert abs(np.sum(ccat)-cice) < eps0, f"1. ={htot} hice={hice} do not match"

  # Check if thickiest ice needs adjustment
  # other cats should be good
  # if hcat is too low - move from lower cats
  # if hcat is too high - increase conc, reduce hcat and adjust conc everywhere to keep cice
  hmin = ICAT[ncat-1]
  hmax = ICAT[ncat] - eps0
  hcat_k = hcat[ncat-1]
  adj_high = hcat_k > hmax # adjust excesss ice in thick. cat.
  adj_low = hcat_k < hmin  # adjust too low (and negative) thkn 
  if adj_high:
    hcat, ccat = adjust_high_thkn(hcat, ccat, ICAT, ncat, hice, eps0, verbose=True)
  elif adj_low:
    hcat, ccat = adjust_high_thkn(hcat, ccat, ICAT, ncat, hice, eps0, verbose=True)

  # Check final distribution: 
  htot = np.sum(ccat*hcat) # should be conserved 
  ctot = np.sum(ccat)
 
  print(f"Input:          hice={hice:.3f}, cice={cice:.3f}")
  print(f"After redistr:  htot={htot:.3f}, ctot={ctot:.3f}")
#  assert abs(htot-hice) < err_lim, f"3. htot={htot} hice={hice} do not match"
#  assert abs(ctot-cice) < err_lim, f"3. ctot={ctot} cice={cice} do not match"

  return hcat, ccat
  
def fcast_mo_to_cal(MMI, MF):
  """
    Find calendar month corresponding to the forecast month=MF initialized in MMI
  """
  mfcast = np.arange(MMI,MMI+12)
  mfcast = np.where(mfcast > 12, mfcast-12, mfcast)
  mcal = mfcast[MF-1]

  return mcal

def cal_mo_to_fcast(MMI, MM):
  """
    Find forecast month (lead time) corresponding to the calend. mo MM for the
    forecast initialized in MMI
  """
  mfcast = np.arange(MMI,MMI+12)
  mfcast = np.where(mfcast > 12, mfcast-12, mfcast)
  if np.any(mfcast == MM):
    imo = (mfcast == MM).argmax()
  else:
    raise Exception(f"Could not find fcast month for calend mo={MMI}")
  mf = imo+1  # f/cast month number

  return mf

def cal_months_forecast(MMI, YRI=1, nyrs=1):
  """
    Create a list of calendar months for a forecast
    initialized on month MMI
    if YRI > 0 - also return a list of years
  """
  MM = np.arange(MMI,MMI+12)
  YY = MM.copy()*0 + YRI
  YY = np.where(MM>12, YY+1, YY)
  MM = np.where(MM>12, MM-12, MM)
  MCAL = MM.copy()
  YCAL = YY.copy()
  for kk in range(2, nyrs+1):
    MCAL = np.append(MCAL,MM)
    YCAL = np.append(YCAL,YY+kk-1)

  return YCAL, MCAL

def read_SPEAR_iconc_clim_interp(YR, MMI, ens_nmb, nyrs_clim=5):
  """
    Read monthly ice conc. clim from SPEAR init = MMI
    for a given year YR
  """
  import xarray
  if nyrs_clim == 5:
    ICLIM=[[1993,1997],[1995,1999],[2000,2004],[2005,2009],[2010, 2014],[2015,2019],[2016,2020]]

  ICLIM=np.array(ICLIM)
  iclm = np.where((ICLIM[:,0] <= YR) & (ICLIM[:,1] >= YR))[0][0]
  YRS,YRE = ICLIM[iclm,:]

  pthpkl = '/work/Dmitry.Dukhovskoy/anls_output/spear_ice'
  fclim = f'spear_interpNEP_siconc_clim_{YRS}_{YRE}_MI{MMI:02d}e{ens_nmb:02d}.nc'
  dfclim = os.path.join(pthpkl,fclim)
  print(f'Opening {dfclim}')
  dset = xarray.open_dataset(dfclim)
  
  return dset
   
def read_NSIDC_iconc_clim_interp(YR, nyrs_clim=5):
  """
    Read monthly ice conc. clim from NSIDC interpolated to NEP fields
    for a given year YR
  """
  import xarray
  if nyrs_clim == 5:
    ICLIM=[[1993,1997],[1995,1999],[2000,2004],[2005,2009],[2010, 2014],[2015,2019],[2016,2020]]

  ICLIM=np.array(ICLIM)
  iclm = np.where((ICLIM[:,0] <= YR) & (ICLIM[:,1] >= YR))[0][0]
  YRS,YRE = ICLIM[iclm,:]

  # NEP grid:
  fyaml = 'paths_seasfcst.yaml'
  with open(fyaml) as ff:
    pthseas = safe_load(ff)

  pthclim = pthseas['ALL']['dirnsidc_clim']
  flclim  = f'NSIDC_NRT_interpNEP_iconc_clim_{YRS}_{YRE}.nc'
  dflclim = os.path.join(pthclim,flclim)

  print(f'Opening {dflclim}')
  dset = xarray.open_dataset(dflclim)

  return dset

def calc_iconc_ithkn_mnthmean(dnmb0,pthtest,varnm, prfx='', ndav=5, outfld='icem'):
  """
    From N-day av. (ndav) fields compute monthly mean fields of ice conc 
    ice thickness (ice volume/m2)
    From SIS2 simulations
    dnmb0 - any date in the month
  """
  import mod_time as mtime
  import mod_anls_seas as manseas
  dvR   = mtime.datevec(dnmb0)
  YYR   = dvR[0]
  MMR   = dvR[1]
  mday1 = int(mtime.datenum([YYR,MMR,1]))
  mday2 = mday1 + int(mtime.month_days(MMR,YYR))-1

  Asum = None
  icc  = 0
  for dnmbR in range(mday1+1,mday2+1,ndav):
    # Find closest output:
    YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthtest, dnmbR, fld=outfld)
    dv0  = mtime.datevec(dnmb0)
    YR0, MM0, DD0 = dv0[:3]
    jday0   = int(mtime.date2jday([YR0,MM0,DD0]))

    if MM0 != MMR:
      print(f' found month {MM0} requested {MMR}, skipping ...')
      continue

    if len(prfx) > 0:
      flice_name = f'{prfx}.icem_{YR0}_{jday0:03d}.nc'
    else:
      flice_name  = f'icem_{YR0}_{jday0:03d}.nc'
    dfsis2 = os.path.join(pthtest, flice_name)

    print(f'Reading {YR0}/{MM0:02d}/{DD0:02d}: {dfsis2}')

    dset   = xarray.open_dataset(dfsis2)

    HIce = dset['sithick'].isel(time=0).data
    CIce = dset['siconc'].isel(time=0).data
    if varnm == 'iconc':
      A2d = CIce
    elif varnm == 'ithkn':
      A2d = CIce*HIce

    if Asum is None:
      Asum = A2d.copy()
    else:
      Asum = Asum + A2d

    icc += 1

  A2d = None
  if icc>1:
    A2d = Asum / icc
  else:
    A2d = Asum.copy()

  return A2d
