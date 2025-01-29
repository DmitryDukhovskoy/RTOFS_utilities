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


