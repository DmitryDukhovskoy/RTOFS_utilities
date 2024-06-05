# To prepare MOM restart:
#  T/S fields
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
importlib.reload(mrgm)

nrun       = "GOFS3.0"
nexpt      = 19.0
expt       = f"{nexpt:2.1f}"
YR         = 1993
jday       = 1
hr         = 15
fldint     = "salin"  # temp salin thknss
grid_shape = 'symmetr'   # MOM grid: symmetr/nonsymmetr
f_thknss   = False     # interpolate sum(A*dH)/dH(j,i) <- errors due to bottom val
if not (fldint == "temp" or fldint == "salin" or fldint == 'thknss'):
  raise Exception ('This code for T or S, change fldint')

pthout  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_restart/'
pthtmp  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_tmp/'
pthgofs = '/work/Dmitry.Dukhovskoy/GOFS3.0/expt_19.0/'

if fldint == "temp":
  grid_var = 'hgrid'
  fldiout  = 'temp'
elif fldint == "salin":
  grid_var = 'hgrid'
  fldiout  = 'salin'
elif fldint == "thknss":
  grid_var = 'hgrid'
  fldiout  = 'lrthknss'

print(f" INTERPOLATING {nrun}-{expt} --> MOM {fldint}")

with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthrun  = dct[nrun][expt]["pthrun"]
pthgrid = dct[nrun][expt]["pthgrid"]
ftopo   = dct[nrun][expt]["ftopo"]
fgrid   = dct[nrun][expt]["fgrid"]

LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
#HH = np.where(HH>=0, np.nan, HH)
Jocn_hycom, Iocn_hycom  = np.where(HH<0)

huge   = 1.e25
rg     = 9806.

# 190_archv.1993_001_15.a
fhycom = f"{nexpt*10:3.0f}_archv.{YR}_{jday:03d}_{hr:02d}"
fina   = os.path.join(pthrun,fhycom) + '.a'
finb   = os.path.join(pthrun,fhycom) + '.b'
A3d, idmh, jdmh, kdmh = mhycom.read_hycom(fina,finb,fldint)
A3d[np.where(A3d>huge)] = np.nan
if fldint == 'thknss':
  A3d = A3d/rg

# Read layer pressures:
if f_thknss:
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
dflintrp = os.path.join(pthtmp,flintrp)


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

Irsb = []
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
  II, JJ = mblnr.sort_gridcell_indx(II,JJ)


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
  if f_thknss:
    qq1 = np.squeeze(dH[:,jj1,ii1])    
    qq2 = np.squeeze(dH[:,jj2,ii2])
    qq3 = np.squeeze(dH[:,jj3,ii3])
    qq4 = np.squeeze(dH[:,jj4,ii4])
#  kk = 0
  la1 = len(np.where(np.isnan(aa1))[0])
  la2 = len(np.where(np.isnan(aa2))[0])
  la3 = len(np.where(np.isnan(aa3))[0])
  la4 = len(np.where(np.isnan(aa4))[0])
# Fill land:
  if la1 == kdmh or la2 == kdmh or la3 == kdmh or la4 == kdmh:
    aa1, aa2, aa3, aa4 = mrgm.fill_land(aa1,aa2,aa3,aa4,HH,A3d,\
                                       JJ,II,Jocn_hycom,Iocn_hycom)
    la1 = len(np.where(np.isnan(aa1))[0])
    la2 = len(np.where(np.isnan(aa2))[0])
    la3 = len(np.where(np.isnan(aa3))[0])
    la4 = len(np.where(np.isnan(aa4))[0])
# Fill bottom
  if la1 > 0:
    aa1 = mrgm.check_bottom(aa1) 
  if la2 > 0:
    aa2 = mrgm.check_bottom(aa2) 
  if la3 > 0: 
    aa3 = mrgm.check_bottom(aa3) 
  if la4 > 0:
    aa4 = mrgm.check_bottom(aa4) 


  if x0 < 0.:
    x0 = x0+360.
  xx = np.where(xx<0., xx+360., xx)
#
# Use cartesian coordinates for mapping
  XV, YV, x0c, y0c = mblnr.lonlat2xy_wrtX0(xx, yy, x0, y0)

# Interpolate with dH taken into account: <--- errors in shallow regions
# due to very thin layer thicknessese, use simple interpolation
# interpolation may not be accurate in shallow regions
# due to potential contamination of below-bottom values 
# mixed with true values in above-bottom
# values on quadrilateral box used for interpolation
  if f_thknss:
    HQT = np.array([aa1*qq1, aa2*qq2, aa3*qq3, aa4*qq4]).transpose()
    QT  = np.array([qq1, qq2, qq3, qq4]).transpose()

  HT  = np.array([aa1, aa2, aa3, aa4]).transpose()

# Map X,Y ---> Xhat, Yhat on reference quadrialteral 
# i.e. map GOFS grid coordinate to a reference quadrilateral 
# to do bilinear interpolation 
#  xht, yht = mblnr.map_x2xhat(xx, yy, x0, y0)  # geogr coord
  xht, yht = mblnr.map_x2xhat(XV, YV, x0c, y0c)   # cartesian coord


# Perform interpolation on reference rectangle, that is 
# similar to interp on actual rectangle
  if f_thknss:
    hqintp = mblnr.bilin_interp1D(phi1, phi2, phi3, phi4, xht, yht, HQT)
    qintp  = mblnr.bilin_interp1D(phi1, phi2, phi3, phi4, xht, yht, QT)
    hintp  = hqintp/qintp
  else:
    hintp  = mblnr.bilin_interp1D(phi1, phi2, phi3, phi4, xht, yht, HT)
#
# Check: abs. values of interpolated values <= original data
  mxHT = np.max(abs(HT), axis=1)
  dmx  = abs(hintp)/mxHT
#  idmx = np.argwhere(dmx<0)
#  idmn = np.argwhere(dmn<0)
  if max(dmx-1.) > 0.1:
    print(f"!!! Min/Max test violated: i/j={ii}/{jj} dlt: {np.max(dmx)}") 
    if max(dmx-1.) > 1.:
      raise Exception("MinMax test: Interp error Check interpolation")

  A3di[:,jj,ii] = hintp

  if kcc%20000 == 0:
    print('Saving to ' + dflintrp)
    with open(dflintrp,'wb') as fid:
      pickle.dump(A3di,fid)

print('END: Saving to '+dflintrp)
with open(dflintrp,'wb') as fid:
  pickle.dump(A3di,fid)


f_plt = False
if f_plt:
  plt.ion()

  if x0 < 0.:
    x0 = x0+360.
  xx = np.where(xx<0., xx+360., xx)

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])
  ax1.plot(x0,y0,'r*')
  ax1.plot(xx,yy,'o')
#  ax1.plot(LON[JJ,II],LAT[JJ,II],'o')
 
