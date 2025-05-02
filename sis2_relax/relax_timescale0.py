"""
  Create inverse relaxation e-folding time scale
  Time scale can vary spatially allowing different relaxation
  rates for different parts of the domain

  Relaxation field is created using complex transformation/mapping technique
 
"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
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
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc
from mod_utils_fig import bottom_text

# Select max and min relaxation time scalres, hrs
# max relaxation - strongest, typically along the OBs
# min relaxation - somewhere in the domain where sea ice presents
#                  e.g. Bering strait
# rx_min, ry_min - approximate # of i, j pnts from the ice OBs
#                  i.e., from i=imax to Ber. Str. (342-200)
# relaxation time scales will be going to 0 away from the ice OBs
rate_max_hrs  = 1.                     # max relaxation time
Irate_max_sec = 1./(rate_max_hrs*3600.)  # relaxation rate, s-1

f_save = False           # Save netcdf relax file
check_rlx = True         # Plot relaxation field
check_ref_domain = True  # Plot transformations of the reference domain

rlx_name = 'relax_rate' # name of the variable, should be the same in the SIS_input

btx  = 'relax_timescale.py'

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

# MOM6 NEP topo/grid:
run_name   = 'seasonal_fcst_daily'
pthtopo    = gridfls['MOM6_NEP'][run_name]['pthgrid']
fgrid      = gridfls['MOM6_NEP'][run_name]['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"][run_name]["ftopo"]
outdir     = gridfls['MOM6_NEP'][run_name]['pthoutp']
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
# Hgrid lon. lat:
hlon, hlat  = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

# Define domain in a complex plane with 
# North-East boundaries along the Im axis (strongest relaxation)
# and slowest relax. rate at a distance dist_rlx 
# rpolar  = approximate # of grid points from the NE corner  of the NEP domain
# to ~55N
# Define a rectangular domain D, transformed relax zone is a square = rpolar
rpolarN = 210   # this defines the size of the relaxation zone along N. bndry axis (index space)
rpolarE = 260   # -"- -"- along the E. bndry axis
#rpolar = 240
Y = np.arange(-rpolarN,rpolarE+1)
X = np.arange(0,idm)
XR = np.arange(-X[-1],1)     # for rotated domain that is in the X<0 plane
YR = np.arange(-Y[-1],-Y[0]+1)  # rotated domain 
idmD = len(X)
jdmD = len(Y)

# Specify the number of grid points from the boundary to keep max relaxation
# Define a reference rectangular domain with relaxation decaying from the left bndry to the right
# Decrease the width of the region of strong relaxation along the N. bndry going south
# and increase the width along the E. bndry going south to have >0 relax over Eastern Bering shelf
irmax = 90    # where to begin exponential decay of the relaxation
ieN = 135     # start pnt (in domain D indices) where the width of max rlx starts changing
isN = 10 
irmaxN = 1    # width of max rlx at the S. end of N. boundary
isE = jdmD-230
ieE = jdmD
irmaxE = 120

# Adjust exponential decay of the relaxation off the Y=0 axis:
f_adj = True
sgmx = idm/4  # controls exponential decay in the Gaussian, decrease the denom. to slow the decay
RLX = np.zeros((jdmD, idmD))
for jj in range(jdmD):
  irmax0 = irmax
  if f_adj:
    if jj>=isN and jj<=ieN:
      irmax0 = irmaxN + int((irmax-irmaxN)/(ieN-isN-1)*(jj-isN))
    elif jj<isN:
      irmax0 = irmaxN
    elif jj>=isE and jj<=ieE:
      irmax0 = irmax + int((irmaxE-irmax)/(ieE-isE-1)*(jj-isE))

  aa = np.exp(-((X-irmax0)**2/sgmx**2))
  aa = aa/np.max(aa)
  aa[:irmax0] = aa[irmax0]
  RLX[jj,:] = aa*Irate_max_sec

# Perform 1st mapping using z**1/2
# Use the fact that mapped domain is symmetric wrt real axis X
RMAP1 = np.zeros((jdmD, idmD))*np.nan
RMAP2 = RMAP1.copy()*np.nan
jD0 = np.argmin(np.abs(Y))
for ii in range(idmD):
  for jj in range(jdmD):
    # Y Distance wrt to jD0:
    xx = X[ii]
    yy = Y[jj]
    RR  = np.sqrt(xx**2+yy**2)
    phi = np.arctan2(yy,xx)
    # Note that under the 1st mapping, RR should be sqrt(RR) 
    # To keep RR unchanged during the transformation and keep the same grid dim
    # apply another transformation z=z*|z|, i.e. sqrt(RR) --> RR
    x_map = RR*np.cos(phi/2)
    y_map = RR*np.sin(phi/2)
    imap = np.max(np.where(X<=x_map)[0])
    jmap = np.max(np.where(Y<=y_map)[0])
    RMAP1[jmap, imap] = RLX[jj,ii]

    # Perform 2nd mapping - rotation by 225 degree angle:
    theta = 225*np.pi/180.
    RRmap1 = np.sqrt((x_map)**2 + (y_map)**2)
    phi1 = np.arctan2(y_map,x_map)
    rmap0 = RMAP2[jD0,jD0]    # debugging
    x_rot = RRmap1*np.cos(phi1+theta)
    y_rot = RRmap1*np.sin(phi1+theta)
    irot = np.min(np.where(XR>=x_rot)[0])
    jrot = np.min(np.where(YR>=y_rot)[0])
    #RMAP2[jrot,irot] = RLX[jj,ii]
    # To avoid gaps: fill neighboring grid cells
    ip1 = np.min([irot+1, idmD])
    im1 = np.max([irot-1, 0])
    jp1 = np.min([jrot+1, jdmD])
    jm1 = np.max([jrot-1, 0])
    RMAP2[jm1:jp1,im1:ip1] = RLX[jj,ii]

# Imbed transformed relax. field into relaxation array PSI
jdmR, idmR = RMAP2.shape
# subset part of the domain that contains transformed field
# domain: Y<=0 and x<=0
iyax0 = max(np.where(YR<=0)[0])
ixax0 = max(np.where(XR<=0)[0])
AA = RMAP2[:iyax0+1,:ixax0]    
AA = np.where(np.isnan(AA), 0., AA)
jdmA, idmA = AA.shape
isD = idm-idmA
jsD = jdm-jdmA
RLXIS = np.zeros((jdm,idm))
RLXIS[jsD:,isD:] = AA    # relaxation rate, s-1
# No relaxation in the Gulf of Alaska:
RLXIS[:575,170:] = 0.0


# Add land mask and make southern domains = 0
lat_cut = 60.
RLXIS = np.where(HH>=0, 0.0, RLXIS)
RLXIS = np.where(hlat<lat_cut, 0.0, RLXIS)

# For checking, relaxation time, hrs:
RLXHR = RLXIS.copy()
RLXHR = np.where(RLXHR==0., np.nan, RLXHR)
RLXHR = 1./RLXHR * 1/3600.

# Write relax time scale:
dflstat = os.path.join(pthtopo, 'ocean_static.nc')
ds_stat = xarray.open_dataset(dflstat)
ds_rlx  = ds_stat['wet'].copy()
ds_rlx.name = rlx_name
ds_rlx *= RLXIS
ds_rlx = ds_rlx.to_dataset()
ds_rlx[rlx_name].attrs['units'] = 's-1'
ds_rlx[rlx_name].attrs['cell_method'] = 'time: point'

cwd = os.getcwd()
ds_rlx.attrs["history"] = f"Created {cwd}/{btx}"


if f_save:
  encoding = {rlx_name: {'_FillValue': None}}
  flout = f'relax_rate_{int(rate_max_hrs):03d}hrs.nc'
  pthsis = gridfls['MOM6_NEP'][run_name]['pthsis']
  dflout = os.path.join(pthsis, flout)

  print(f'Saving SIS2 relaxation time scale --> {dflout}')
  ds_rlx.to_netcdf(
      dflout,
      format='NETCDF3_64BIT',
      engine='netcdf4',
      encoding=encoding
  )


