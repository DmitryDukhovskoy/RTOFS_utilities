# To prepare MOM restart:
#  U/V fields
#  use archv - staggered U V and need to add barotropic U
#
# 1st step: interpoalte fields horizontally 
# from HYCOM archv file onto MOM grid
#
# Need to have MOM6 to reanalysis grid mapping indices
# e.g. for GOFS use hinterp_gofs2mom_gmapi.py
# to create gmapi 
# HYCOM vs MOM staggered grids:
# see: HYCOM-tools/topo/src/grid_mom6.f:
# MOM lon/lat is on a supergrid !
# for MOM: i=1, 1.5, 2, 2.5, ... j=1, 1.5, 2, 2.5, ...
# so that MOM(i=2, j=2) ==> i=1.5, j=1.5 HYCOM grid cell, i.e. in the
# middle of the grid cell
#
# Keep fields on HYCOM vertical grid
#
# Note: all below-bottom values have to be filled with above values to avoid nans
# 
# Plot NEP domain on GOFS grid
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
#import struct
import datetime
import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
import timeit
import yaml
from netCDF4 import Dataset as ncFile

#PPTHN = '/home/Dmitry.Dukhovskoy/python'
PPTHN = []
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mom6util
import mod_regmom as mrgm
#import mod_valid_utils as mvutil
importlib.reload(mcmp)

nrun       = "GOFS3.0"
nexpt      = 19.0
expt       = f"{nexpt:2.1f}"
YR         = 1993
jday       = 1
hr         = 15
fldint     = "u-vel."  # temp salin u-vel. v-vel. layer thkn: u & v - add brtrp vel!
grid_shape = 'symmetr'   # MOM grid: symmetr/nonsymmetr

if not (fldint == "u-vel." or fldint == "v-vel."):
  raise Exception ('This code for T or S, change fldint')

pthout  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_restart/'
pthgofs = '/work/Dmitry.Dukhovskoy/data/GOFS3.0/expt_19.0/'

if fldint == "temp":
  grid_var = 'hgrid'
  fldiout  = 'temp'
elif fldint == "salin":
  grid_var = 'sgrid'
  fldiout  = 'salin'
elif fldint == "u-vel.":
  grid_var = 'ugrid'
  fldiout  = 'uvel'
elif fldint == "v-vel.":
  grid_var = 'vgrid'
  fldiout  = 'vvel'
elif fldint == "thknss":
  grid_var = 'hgrid'
  fldiout  = 'lrthknss'


with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthrun  = dct[nrun][expt]["pthrun"]
pthgrid = dct[nrun][expt]["pthgrid"]
ftopo   = dct[nrun][expt]["ftopo"]
fgrid   = dct[nrun][expt]["fgrid"]

LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
HH = np.where(HH>=0, np.nan, HH)

huge   = 1.e25
rg     = 9806.

# 190_archv.1993_001_15.a
fhycom = f"{nexpt*10:3.0f}_archv.{YR}_{jday:03d}_{hr:02d}"
fina   = os.path.join(pthrun,fhycom) + '.a'
finb   = os.path.join(pthrun,fhycom) + '.b'
A3d, idmh, jdmh, kdmh = mhycom.read_hycom(fina,finb,fldint)
A3d[np.where(A3d>huge)] = np.nan

# Read layer pressures:
dH, _, _, _ = mhycom.read_hycom(fina,finb,'thknss')
dH          = dH/rg
dH          = np.where(dH>huge, np.nan, dH)
dH          = np.where(dH<1.e-12, 1.e-12, dH)  # avoid 0 dH for interpolation

# Saved gmapi:
dnmb    = mtime.jday2dnmb(YR, jday)
dv      = mtime.datevec(dnmb)
ds      = mtime.datestr(dnmb)
MM      = dv[1]
DM      = dv[2]
rdate   = f"{YR}{MM:02d}{DM:02d}"
flgmap  = f"gofs2mom_nep_gmapi-{grid_var}_{grid_shape}.pkl"
dflgmap = os.path.join(pthout,flgmap)

# Output interpolated fields on hycom layers:
flintrp  = f"gofs2mom_nep_hrzi-{fldiout}_{grid_shape}_{rdate}.pkl"
dflintrp = os.path.join(pthout,flgmap)


#
# Read MOM6 NEP domain nx=342, ny=816
pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom) 
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom) 
hlon, hlat  = mom6util.read_mom6grid(dfgrid_mom, grid=grid_shape, grdpnt=grid_var)
HHM         = mom6util.read_mom6depth(dftopo_mom) 
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]
# Convert to -180/180:
hlon   = np.where(hlon > 180.0, hlon-360., hlon)

Iall = np.where(HHM.flatten()<=0.)[0]
nall = Iall.shape[0]

# gmapi mapping indices MOM6 <---> GOFS:
with open(dflgmap,'rb') as fid:
  INDX, JNDX = pickle.load(fid)


# Use bilinear interpolation 
# Points are mapped onto a reference square for interpolation then 
# Perform interpolation on reference rectangle, that is 
# similar to interp on actual rectangle
# The advantage of using a reference quadrilateral is that basis
# functions computed once and used for all grid points - saved lots of 
# computational time
#     See p. 83-90, Gockenbach, "Understanding Finite-Element Methods" textbook
import mod_bilinear as mblnr
importlib.reload(mblnr)

# Find basis functions for a reference rectangle:
phi1,phi2,phi3,phi4 = mblnr.basisFn_RectRef()
phi_basis           = np.array([phi1, phi2, phi3, phi4]).transpose() # basis funs in columns

A3di = np.zeros((kdmh,jdm,idm))*np.nan
kcc  = 0
tic  = timeit.default_timer()
ticR = timeit.default_timer()
print(f'Interpolating HYCOM GOFS --> MOM6 {fldint} {grid_shape} {grid_var}')

for ikk in range(nall):
  I1 = Iall[ikk]
  jj, ii = np.unravel_index(I1,HHM.shape)
  kcc += 1

# If continue from saved:
  if len(Irsb) > 0:
    if np.min(abs(Irsb-I1)) == 0:
      if (kcc % 1000) == 0 or I1 == Irsb[-1]:
        print(f" ---> Skipping {kcc/nall*100.:3.1f}%...")
      continue

  if (kcc % 2000) == 0:
    toc    = timeit.default_timer()
    pprc   = kcc/nall*100.
    dtic   = (toc-ticR)/60.
    tictot = (toc-tic)/60.
    print(f"  {pprc:4.1f}% done dt={dtic:6.2f} min, ttot={tictot:7.1f} min")
    ticR = timeit.default_timer()

  II = np.squeeze(INDX[jj,ii,:])
  JJ = np.squeeze(JNDX[jj,ii,:])

  x0  = hlon[jj,ii]
  y0  = hlat[jj,ii]

  ii1, jj1 = II[0], JJ[0]
  ii2, jj2 = II[1], JJ[1]
  ii3, jj3 = II[2], JJ[2]
  ii4, jj4 = II[3], JJ[3]

  xx  = np.array([LON[jj1,ii1], LON[jj2,ii2], LON[jj3,ii3], LON[jj4,ii4]])
  yy  = np.array([LAT[jj1,ii1], LAT[jj2,ii2], LAT[jj3,ii3], LAT[jj4,ii4]])
  aa1 = np.squeeze(A3d[:,jj1,ii1])    
  aa2 = np.squeeze(A3d[:,jj2,ii2])
  aa3 = np.squeeze(A3d[:,jj3,ii3])
  aa4 = np.squeeze(A3d[:,jj4,ii4])
  qq1 = np.squeeze(dH[:,jj1,ii1])    
  qq2 = np.squeeze(dH[:,jj2,ii2])
  qq3 = np.squeeze(dH[:,jj3,ii3])
  qq4 = np.squeeze(dH[:,jj4,ii4])
#  kk = 0
  if len(np.where(aa1 == np.isnan)[0]) > 0:
    aa1 = mrgm.check_bottom(aa1) 
  if len(np.where(aa2 == np.isnan)[0]) > 0:
    aa2 = mrgm.check_bottom(aa2) 
  if len(np.where(aa3 == np.isnan)[0]) > 0:
    aa3 = mrgm.check_bottom(aa3) 
  if len(np.where(aa4 == np.isnan)[0]) > 0:
    aa4 = mrgm.check_bottom(aa4) 

# Interpolate with dH taken into account:

  HQT = np.array([aa1*qq1, aa2*qq2, aa3*qq3, aa4*qq4]).transpose()
  QT  = np.array([qq1, qq2, qq3, qq4]).transpose()
  HT  = np.array([aa1, aa2, aa3, aa4]).transpose()

# Map X,Y ---> Xhat, Yhat on reference quadrialteral 
# i.e. map GOFS grid coordinate to a reference quadrilateral 
# to do bilinear interpolation 
  xht, yht = mblnr.map_x2xhat(xx, yy, x0, y0)
# Perform interpolation on reference rectangle, that is 
# similar to interp on actual rectangle
  hqintp = mblnr.bilin_interp1D(phi1, phi2, phi3, phi4, xht, yht, HQT)
  qintp  = mblnr.bilin_interp1D(phi1, phi2, phi3, phi4, xht, yht, QT)
#  hintp = mblnr.bilin_interp1D(phi1, phi2, phi3, phi4, xht, yht, HT)
  hintp  = hqintp/qintp
#
# Check:
  if np.max(abs(hintp)) > np.max(abs(HT)):
    raise Exception("Max test violated: Check interpolation")

  A3di[:,jj1,ii1] = hintp

  if kcc%10000 == 0:
    print('Saving to ' + dflintrp)
    with open(dflintrp,'wb') as fid:
      pickle.dump(A3di,fid)

print('END: Saving to '+dflgmap)
with open(dflintrp,'wb') as fid:
  pickle.dump(A3di,fid)


f_plt = False
if f_plt:
  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])

  ax1.plot(x0,y0,'r*')
  ax1.plot(LON[JJ,II],LAT[JJ,II],'o')
 
