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

def redistribute_hice(hice, cice, ICAT=[], eps0 = 1.e-10, itd_method='equal'):
  """
    Redistribute ice thickness (grid cell mean), 1 pnt, 
    by ice categories such that the mean grid 
    ice thickness is conserved, i.e.
    h_mean = sum(i_cat)(h_cat(i)*conc_cat(i))

    ice thickn. distribution methods:
    - simple = place all ice in the category that corresponds hice (grid cell mean value)
    - gauss = start with ice concentration aplying Gaussian filter on cice in the cat. = hice
    - equal = start with ice concentration equally distributed over the cats.

  """
  if len(ICAT) == 0:
    ICAT = np.array([1.0e-10, 0.1, 0.3, 0.7, 1.1])

  ncat  = len(ICAT)
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
  hcat = np.zeros((ncat))
  ccat = np.zeros((ncat))
  ccat[indx_cat] = cice
  
  match itd_method:
    case('simple'):
      hcat[indx_cat] = hice/cice
      return hcat, ccat
    case('gauss'):
      ccat = gauss_distr1D(chice_cat, indx_cat, sgm=1.2)
    case('equal'):
      ccat = np.zeros((ncat)) + cice/ncat

  # Adjust too high concentrations (total > cice and ccat[k]>1):
  cfctr = np.sum(ccat)/cice  # reduction factor to adjust concentration:
  ccat = ccat/cfctr

  # Initial ice distribution = max ice cat. thickness except for the thickest one
  # the thickest category is used as a loan 
  hcat = np.zeros((ncat))
  hcat[0:ncat-1] = ICAT[1:ncat]-eps0

  # Check if ice needs to be added or removed from the ice thickest cat:
  # if hice thickest > hice_lim - increase concentration in the lower categories
  # if hice thickest < ICAT[-1] - increase conce in the lower cat and try to make 0 thickets ice
  #                               if cannot increase conc in the lower cats, - add more ice to thickest
  #                               and reduce ice in the lower ice cats.

#STOP HERE

  # Perform Ice Thickness Distribution to readjust ice thicknesses within ice cat.
  # Excess hice - move up, low hice - move from upper ICAT
  hcat_adj = hcat.copy()
  for kk in range(ncat-1):
    dlt_btm = (ICAT[kk]+eps0) - hcat[kk]
    dlt_up  = ICAT[kk+1] - (hcat[kk]+eps0)
    if dlt_btm > 0.:
      hcat_adj[kk] = hcat_adj[kk] + dlt_btm
      dlt_hc = abs(dlt_btm)*ccat[kk]      # surplus fractional ice mass being redistributed 
      hcat_adj[kk+1] = hcat_adj[kk+1] - dlt_hc/ccat[kk+1]
    elif dlt_up < 0.:
      hcat_adj[kk] = hcat_adj[kk] + dlt_up
      dlt_hc = abs(dlt_up)*ccat[kk]      # surplus fractional ice mass being redistributed 
      hcat_adj[kk+1] = hcat_adj[kk+1] + dlt_hc/ccat[kk+1]



  # Distribute grid cell average ice (concentration * hice grid) 
  # by categ. then use bin-average value for hice(cat) to deduce conc(cat)
  #using linear spreading centered in the cat = hice:
#  for iflt in range(2):
#    ccat = linear_distr1D(ccat)

  chice = cice*hice
  chice_cat = np.zeros((ncat))
  chice_cat[indx_cat] = chice
  chice_cat = gauss_distr1D(chice_cat, indx_cat, sgm=1.2)

  # Conserve grid cell mean ice thickness:
  if abs(np.sum(chice_cat)-hice) > 1.e-15:
    chice_cat = hice*chice_cat/np.sum(chice_cat)   

  # Redistribute ice thickness for all categories then adjust 
  # changing the thickest one
  cat_diff = np.diff(ICAT)
  cat_diff = np.append(cat_diff,2.0)  # thickest ice
  # assume ice thicknesses = mid-values of ice cat. 
  hcat = ICAT+0.5*cat_diff[kk]
  ccat = chice_cat/hcat
  ccat = np.where(ccat<eps0, eps0, ccat)  # avoid 0 concentrations

  # Check the cell-mean hice:
  hice_new = np.sum(hcat*ccat)

  # Adjust too high concentrations (total > cice and ccat[k]>1):
  cfctr = np.sum(ccat)/cice  # reduction factor to adjust concentration:
  ccat = ccat/cfctr
  hcat = hcat*cfctr

  # Perform Ice Thickness Distribution to readjust ice thicknesses within ice cat.
  # Excess hice - move up, low hice - move from upper ICAT
  hcat_adj = hcat.copy()
  for kk in range(ncat-1):
    dlt_btm = (ICAT[kk]+eps0) - hcat[kk]
    dlt_up  = ICAT[kk+1] - (hcat[kk]+eps0)
    if dlt_btm > 0.:
      hcat_adj[kk] = hcat_adj[kk] + dlt_btm 
      dlt_hc = abs(dlt_btm)*ccat[kk]      # surplus fractional ice mass being redistributed 
      hcat_adj[kk+1] = hcat_adj[kk+1] - dlt_hc/ccat[kk+1]
    elif dlt_up < 0.:
      hcat_adj[kk] = hcat_adj[kk] + dlt_up
      dlt_hc = abs(dlt_up)*ccat[kk]      # surplus fractional ice mass being redistributed 
      hcat_adj[kk+1] = hcat_adj[kk+1] + dlt_hc/ccat[kk+1]




    
  

  # Check ice thickness: should be conserved
  dlt_hice = ICAT[0]
  hice_err = hice_new-hice
  if abs(hice_err) > dlt_hice:
    print(f"Add ice adjustment, dlt hice={hice_err:14.6f}")

  return hcat, ccat
  



