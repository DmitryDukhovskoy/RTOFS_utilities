# Calculate N2 following Chelton 1996
# N2 is calculated at the middle of the depth interval
# density is adiabatically adjusted to the mid-grid depth
#
# Using WKB approach described in Chelton 1996
##reset -f 
from IPython import get_ipython
get_ipython().magic('reset -sf')

import os
import numpy as np
import matplotlib.pyplot as plt
import torch
import sys
import pdb
import netCDF4
import importlib
from netCDF4 import Dataset as ncFile

sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
import mod_utils_fig
from mod_utils_fig import bottom_text
from mod_utils_fig import plot_fld2D

grd=0.25
#ftmp='/data/ncei2/w18C/analysis/all_0/{0}/mean/M02001'.format(grd)
#fsal='/data/ncei2/w18C/analysis/all_0/{0}/mean/s013'.format(grd)

#nlat=720
#nlon=1440
#nz=100
#dt=np.dtype((np.float32,(nlat,nlon)))

btx = 'calcN2.py'

urlT='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/temperature/decav/{0}/'.\
      format(grd)
urlS='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/salinity/decav/{0}/'.\
      format(grd)

tfnm='woa18_decav_t13_04.nc'
sfnm='woa18_decav_s16_04.nc'

#
# Reloading modules:
def reload_mod(__name__):
  try:
      del sys.modules[__name__]
  except AttributeError:
      pass
  import __name__



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


def plot_fld(aa,tvar,z0,cl1,cl2):
  jdm=aa.shape[0]
  idm=aa.shape[1]

  plt.ion()
  cmpS = plt.cm.get_cmap('turbo')
  fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im = ax1.pcolormesh(aa,cmap=cmpS,shading='flat')
  ax1.axis('scaled')
  ax1.set(xlim=(0,idm),ylim=(0,jdm))
  ax1.autoscale(enable=True, axis='both', tight=True)

  ax1.set_title('Var {0}, z={1}'.format(tvar,z0))
  fig1.colorbar(im,ax=ax1,orientation='horizontal')

  im.set_clim(cl1,cl2)


  plt.show()
  figout='plot_field.jpg'
  plt.savefig(figout)
  plt.close(fig1)

  wwwpth='/net/www-dev/www/ncei_wod/'
  wwwfig='pthn_fig.jpg'
  cmd = ('mv {0} {1}{2}'.format(figout,wwwpth,wwwfig)) 
  os.system(cmd)  


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


# Pick some location:
aa=T[iz,:,:].squeeze()
ss=S[iz,:,:].squeeze()
LMsk=aa.copy();
LMsk[np.isfinite(aa)]=1.
LMsk[~np.isfinite(aa)]=0.

ii=600
jj=400
zz=ZZ
tt=T[:,jj,ii].squeeze()
ss=S[:,jj,ii].squeeze()

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
# Test: adiabatic T gradient
# t=20.0, s=34.0, pres=10.0 db --> adiab T grad = 1.83342653e-4
  t=20.0
  s=34.0
  pres=10.0
  adtg = adiab_Tgrad(s,t,pres)
  print('adiab T grad(s={0},t={1},p={2}db) = {3}'.
         format(s,t,pres,adtg))

  # Test: potential T from pressure P1 to Pref
  # for t=20, s=34, pr1=1000db, pref=2000db
  # pot. temp = 20.1958294083943 ...
  t=20.0
  s=34.0
  pr1=1000.0
  pref=2000.0
  ptmp = sw_ptmp(s,t,pr1,pref)
  print('Potential T(s={0},t={1},p1={2}db, pref={3}) = {4}'.
         format(s,t,pr1,pref,ptmp))


# Test rho pure water at P=0
# at t=24.4 rho_S0 = 997.198518670105
  t=24.4
  rho_S0 = sw_smow(t)
  print('Pure water dens(T={0}) = {1}'.format(t,rho_S0)) 

  # Test depth ---> pressure conversion
  # should get 7500.00 db
  #pres_db, pres_pa =sw_press(7321.45,30)

#
# Test sea water dens at the surface P=0 
# t=11.5, s=45.2, rho=1034.63853438209
  t=11.5
  s=45.2
  rho_P0 = sw_dens0(s,t)
  print('Surface Water dens(S={0}, T={1}) = {2}'.format(s,t,rho_P0)) 


# Secant bulk formula
# t=11.5, s=45.2, P=5000 db
# K secnat = 25057.5763828814
  t=11.5
  s=45.2
  p=5000.0
  Kscnt = sw_seck(s,t,p)
  print('Secant bulk modulus(S={0}, T={1}, P={2}) = {3} bar'.format(s,t,p,Kscnt))
  
#
# in situ water dens
# t=11.5, s=45.2, P=5000 db
# 1055.7041012411848 kg/m3
  t=11.5
  s=45.2
  p=5000.0
  rho = sw_dens(s,t,p)
  print('rho(S={0}, T={1}, P={2}) = {3} kg/m3'.format(s,t,p,rho))

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
# If N2 < 0 - density inversion happnes in some ~homogeneous layers
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


#

# Plot T
#cl1=-2.
#cl2=10
#ctitle='T surf'
#plot_fld2D(t1,ctitle,cl1,cl2,btx=btx) # quick plot on the web


kk=0
aa=N2[kk,:,:].squeeze()
zz0=Z_phi[kk]
ctitle='N2 zz0={0}'.format(zz0)

cl1=0.
cl2=0.0005
plot_fld2D(aa,ctitle,cl1,cl2,btx=btx) # quick plot on the web
#

#
 

# Solve Storm-Liouville e/value problem for 1 profile
ii=600
jj=400
omg=7.29e-5
phi=lat[jj]
fcor=2.*omg*np.sin(np.deg2rad(phi))
N2z = N2[:,jj,ii].squeeze()
#
# Find bottom:
kbtm = np.where(np.isnan(N2z))[0].min()-1
zbtm = ZZ[kbtm]
#
# Fill below bottom:
N2z[kbtm+1:]=N2z[kbtm]

#
# Filter 5-point running mean
import mod_solver
importlib.reload(mod_solver)
from mod_solver import runmn
N2zF = runmn(N2z,ZZ,mnwnd=9)


#
# Convert to cycles/hr
N_hrz = np.sqrt(N2z)
N_chr = 3600.*N_hrz   # cycles per hr
NF_chr = 3600.*np.sqrt(N2zF)
f_plt=0
if f_plt==1:
  plt.ion()
  fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
  plt.plot(N_chr,Z_phi)  
  plt.plot(NF_chr,Z_phi)


# Create Matrix A with Dk, Dk+1 for 2nd derivative of N2
from mod_solver import form_mtrxA
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


#K=N2A[::5,:]

#
# For testing - form Matrix A unfolding
# the 3 -elemnts form and find eigenvalues/vectors
from mod_solver import eig_unfoldA
W, V = eig_unfoldA(kbtm,N2A)

#
# Obtain eig.value from WKB approximation solution
from mod_solver import eig_wkb
wm = eig_wkb(kbtm,N2z,ZZ,fcor,mode=1)

# Calculate
# Choose lmbd1 = as min(abs(W))
im=np.where(np.abs(W) == np.min(np.abs(W)))[0]
lmbd=np.sqrt(-1./(W[im]*fcor**2))

f_plt=0
if f_plt>0:
  plt.ion()
  fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)






