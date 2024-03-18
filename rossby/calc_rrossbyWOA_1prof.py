# The code is similar to calc_rrossby_etopoWOA.py
# but computes Rossby radius for 1 profile from the WOA
# specified by lon/lat (or grid indices) 
# 
# 1st baroclinic Rossby radius 
# 2 approaches: WKB-theory and numerically solve eig/value problem
#
# Use bottom topography for depths deeper than WOA last depth level (-5500m)
#  Topography interpolated from ETOPO2
#
#
# Calculate N2 following Chelton 1996
# N2 is calculated at the middle of the depth interval
# density is adiabatically adjusted to the mid-grid depth
#
import os
import numpy as np
import matplotlib.pyplot as plt
#import torch
import sys
import pdb
import netCDF4
import importlib
from netCDF4 import Dataset as ncFile
import timeit

sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
import mod_utils_fig
# importlib.reload(mod_utils_fig)
from mod_utils_fig import bottom_text
from mod_utils_fig import plot_fld2D

lon_plt = 137.
lat_plt = 0.0

grd=0.25
if grd==0.25:
  cgrd=4
woa='woa18'
seas=13    # season: 1-12 - monthly, 13-16 - spring, summer, ...

pthout = '/data/ncei2/work/dmitry.dukhovskoy/data_output/'
#ftmp='/data/ncei2/w18C/analysis/all_0/{0}/mean/M02001'.format(grd)
#fsal='/data/ncei2/w18C/analysis/all_0/{0}/mean/s013'.format(grd)

#nlat=720
#nlon=1440
#nz=100
#dt=np.dtype((np.float32,(nlat,nlon)))

btx = 'calc_rrossbyWOA_1prof.py'
urlT='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/temperature/decav/{0}/'.\
      format(grd)
urlS='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/salinity/decav/{0}/'.\
      format(grd)

tfnm='{0}_decav_t{1:02d}_{2:02d}.nc'.format(woa,seas,cgrd)
sfnm='{0}_decav_s{1:02d}_{2:02d}.nc'.format(woa,seas,cgrd)

#

def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].squeeze()
  dmm = np.copy(dmm0)
  return dmm

def lookup_ncvar(nc):
  ii=0
  for var in nc.variables.values():
    ii+=1
    print('--------\n')
    print('Var # {0}'.format(ii))
    print(var)

iz=0 
tvar='t_an'
svar='s_an'

furl=urlT+tfnm
T=read_field(furl,tvar)
T[T>1.e10]=np.nan
kdm=T.shape[0]
jdm=T.shape[1]
idm=T.shape[2]

ZZ=read_field(furl,'depth')
ZZ=-abs(ZZ)
lat=read_field(furl,'lat')
lon=read_field(furl,'lon')

LON=np.zeros((jdm,idm))
LAT=np.zeros((jdm,idm))
for ii in range(idm):
  LAT[:,ii]=lat

for jj in range(jdm):
  LON[jj,:]=lon

furlS=urlS+sfnm
S=read_field(furlS,svar)
S[S>1.e10]=np.nan


# Land-sea mask
aa=T[iz,:,:].squeeze()
ss=S[iz,:,:].squeeze()
LMsk=aa.copy();
LMsk[np.isfinite(aa)]=1.
LMsk[~np.isfinite(aa)]=0.

zz=ZZ

# Find index to plot:
dmm = abs(lat-lat_plt)
jplt = np.argmin(dmm)
dmm = abs(lon-lon_plt)
iplt = np.argmin(dmm)


import mod_swstate
#importlib.reload(mod_swstate)
from mod_swstate import sw_press
from mod_swstate import adiab_Tgrad
from mod_swstate import sw_ptmp
from mod_swstate import sw_dens0
from mod_swstate import sw_smow
from mod_swstate import sw_seck
from mod_swstate import sw_dens

f_test = 0
if f_test > 0:
  from test_sw_state import test_sw
  test_sw(t=20., s=34., pr1=100., pref=3000., p=6000.)
# --------------------------------------
#  Calculate N2 = -g/rho*d2(rho)/dz2 
#  Using in situ T and S - calculate rho relative to the
#  mid-grid depth - following Chelton, 1996
# Use the neutral density gradient method, Chelton et al., 1996

print('Calculating pressure at midpoints ...')
grav = 9.81
rho0 = 1025.0  
Z_phi=np.zeros(kdm-1)  # depths for e/function Phi
N2 = np.zeros((kdm-1,jdm,idm))
n2fill=1.e-8     # missing values
for kk in range(kdm-1):
  print(' Layer {0}'.format(kk))
  z1 = ZZ[kk]
  z2 = ZZ[kk+1]
  t1 = T[kk,:,:].squeeze()
  t2 = T[kk+1,:,:].squeeze()
  s1 = S[kk,:,:].squeeze()
  s2 = S[kk+1,:,:].squeeze()
  Z1=z1*np.ones((jdm,idm))
  Z2=z2*np.ones((jdm,idm))
  p1_db, p1_pa = sw_press(Z1,LAT)  # pressure upper interface
  p2_db, p2_pa = sw_press(Z2,LAT)  # pressure bottom interface
#
# Find mid-point of the layers
  p0_db = 0.5*(p1_db+p2_db)
  z0    = 0.5*(z1+z2)
  Z_phi[kk] = z0
#
# Calculate rho(z1--->z0) with depth reference at midpoint
# and rho(z2--->z0) at midpoint 
  t_z1z0 = sw_ptmp(s1,t1,p1_db,p0_db)
  t_z2z0 = sw_ptmp(s2,t2,p2_db,p0_db)  
  rho_z1z0 = sw_dens(s1,t_z1z0,p0_db)
  rho_z2z0 = sw_dens(s2,t_z2z0,p0_db)

# Calculate d(rho)/dz for z0 - center-difference
  drho_dz = (rho_z1z0 - rho_z2z0)/(z1 - z2)
  N2z0 = -grav/rho0*drho_dz 
  N2[kk,:,:]=N2z0
  print('k={0}, z={1}, min/max N2: {2}, {3}'.
         format(kk,z0,np.nanmin(N2z0),np.nanmax(N2z0)))

#
# If N2 < 0 - density inversion happens in some ~homogeneous layers
# when parcels are brought down from z1 and up from z2
# replace with above N2 or below of surface layer
print('Fixing N2<0')
for kk in range(kdm-1):
  N2z0=N2[kk,:,:].squeeze()
  if kk > 0:
    dmm = N2[kk-1,:,:].squeeze()
    N2z0[np.where(N2z0<0.)] = dmm[np.where(N2z0<0.)]
  else:
    dmm = N2[kk+1,:,:].squeeze()
    N2z0[np.where(N2z0<0.)] = dmm[np.where(N2z0<0.)]

  N2z0[np.where(N2z0<0.)] = n2fill
  N2[kk,:,:]=N2z0

  print('k={0}, z={1}, min/max N2: {2}, {3}'.
         format(kk,z0,np.nanmin(N2z0),np.nanmax(N2z0)))

# ------------
# Read TOPO
# -------------
pthout = '/data/ncei2/work/dmitry.dukhovskoy/topo/woa/'
fout = pthout+'etopo2woa025_bilinear.dat'
ijdm = idm*jdm
fid = open(fout,'rb')
try:
  fid.seek(0)
  HW = np.fromfile(fid,dtype='<f4', count=ijdm)
finally:
  fid.close()
HW = HW.reshape(jdm,idm)

#
# =======================================================
#
# Numerically Solve Sturm-Liouville e/value problem 
#
# =======================================================
import mod_solver
importlib.reload(mod_solver)
from mod_solver import runmn
from mod_solver import form_mtrxA
from mod_solver import eig_unfoldA
from mod_solver import eig_wkb

#RsbNum = LMsk.copy()
#RsbNum = np.where(LMsk==0,np.nan,0.0)
#RsbWkb = RsbNum.copy()

omg = 7.29e-5 # Earth angular velocity
Iocn = np.where(LMsk.flatten()>0)[0]
nocn = Iocn.shape[0]
cc = 0
#tic = timeit.default_timer()
#ticR = timeit.default_timer()
print('Solving Strum-Liouville ...')
#for iocn in range(nocn):
#  I1 = Iocn[iocn]
#  jj, ii = np.unravel_index(I1,LMsk.shape)
#  cc += 1

ii = iplt
jj = jplt

N2z = N2[:,jj,ii].squeeze()
#
# Find bottom:
# Escape coastal regions < 20 m
k1 = np.where(np.isnan(N2z))[0]
if k1.size:
  kbtm = k1.min()-1
#
# Fill below bottom:
  N2z[kbtm+1:]=N2z[kbtm]
  zbtm = ZZ[kbtm+1]  # bottom depth = bottom interface of the near-bottom gr.cell
else:
# Case when Bottom is deeper than WOA last depth level
# extend last grid cell to the actual bottom
#    print('Deep site: ii={0}, jj={1}, Hbtm={2:6.1f}'.\
#          format(ii,jj,HW[jj,ii]))
  kbtm = N2z.shape[0]-1
  zbtm = HW[jj,ii]

# Skip shallow regions
if kbtm < 5:
  print('Shallow location: zbtm={0:.2f} m'.format(zbtm))
  raise NameError('Exit shallow depth')

# Filter N-point running mean
N2zF = runmn(N2z,ZZ,mnwnd=7)

# Create Matrix A with Dk, Dk+1 for 2nd derivative of N2
AA = form_mtrxA(Z_phi,kbtm,zbtm)
#
# Form 1/N2*AA:
# The matrix is cutoff by the bottom depth
ka = AA.shape[0]
na = AA.shape[1]
N2A = np.zeros((ka,na))
for kk in range(ka):
#  n2 = N2z[kk]
  n2 = N2zF[kk]  # filtered profile
  if n2 == 0:
    n2=1.e-12
  N2A[kk,:] = 1./n2*AA[kk,:]

#
# For now: use python eig function, for this
# form Matrix A unfolding
# the 3 -elemnts form and find eigenvalues/vectors
# W - eigenvalues, not sorted our by magnitude!
# V- corresponding eigenvectors (in columns)
W, V = eig_unfoldA(kbtm,N2A)

# Calculate
# Choose lmbd1 = as min(abs(W)) and lmbd1<0
im = np.where(np.abs(W) == np.min(np.abs(W)))[0][0]
#  im = np.where((W < 0.) & (np.abs(W) == np.min(np.abs(W))))[0][0]
if W[im] > 0.:
  print('ERROR: W[im] >0: im={0}, W[im]={1}, zbtm={4:6.1f}, ii={2}, jj={3}'.\
         format(im, W[im], ii, jj, zbtm))

latj = lat[jj]
Rearth = 6371.e3  # Earth R
fcor = 2.*omg*np.sin(np.deg2rad(latj))
betta = (2.*omg*np.cos(np.deg2rad(latj)))/Rearth

if abs(latj) >= 5.0: 
  lmbd = np.sqrt(-1./(W[im]*fcor**2))
else:
  lmbd = (-1./(4.*(betta**2)*(W[im])))**(1./4.)

RsbNum = lmbd*1.e-3  # km
Phi = V[:,im]   # eigenfunction

print(" All Done R = {0:.4f} km".format(RsbNum))

# Test - plot N profiles
# Convert to cycles/hr
f_testN = 0
if f_testN == 1:
  N_hrz = np.sqrt(N2z)
  N_chr = 3600.*N_hrz   # cycles per hr
  NF_chr = 3600.*np.sqrt(N2zF)
  plt.ion()
  fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
  fig1.clf()
  ax1 = plt.axes([0.1, 0.1, 0.35, 0.8]) 
  plt.plot(N_chr,Z_phi)  
  plt.plot(NF_chr,Z_phi)
  ax1.plot([0,max(N_chr)],[zbtm,zbtm],'k--')
  ax1.grid(True)
  ax1.set_xlabel('cycle/hr')
  ctl = 'N^2, Zbtm={0:.1f}, {1}'.format(zbtm,tfnm)
  ax1.set_title(ctl)

  Rrsb = RsbNum
  x0 = lon[ii]
  y0 = lat[jj]
  ww = W[im]
  vvk = V[:,im]
  zzk = Z_phi[0:kbtm+1]
  nzk = zzk.shape[0]
# 
# Add surface and bottom to eig/functions
  zzV=np.zeros(nzk+2)
  zzV[0] = 0.
  zzV[1:nzk+1]=zzk
  zzV[nzk+1]=zbtm

# Add 0 at the ends for eig/functions:
  vvV = np.zeros(zzV.shape[0])
  vvV[1:nzk+1] = vvk

  ax2 = plt.axes([0.55, 0.1, 0.35, 0.8])
  ax2.plot(vvV,zzV,'.-')
  ax2.grid(True)
  ctl2 = 'Eig/vect {1:5.2f}E, {2:5.2f}N, Rr={0:6.0f} km'.\
         format(Rrsb,x0,y0)
  ax2.set_title(ctl2)

  bottom_text(btx,pos=[0.01,0.02])



