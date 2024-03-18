# 1st baroclinic Rossby radius 
#
# Use bottom topography for depths deeper than WOA last depth level (-5500m)
#  Topography interpolated from ETOPO2
# Calculate N2 following Chelton 1996
# N2 is calculated at the middle of the depth interval
# density is adiabatically adjusted to the mid-grid depth
#
import numpy as np
#import importlib
#import mod_swstate
#importlib.reload(mod_swstate)
from mod_swstate import sw_press
from mod_swstate import adiab_Tgrad
from mod_swstate import sw_ptmp
from mod_swstate import sw_dens0
from mod_swstate import sw_smow
from mod_swstate import sw_seck
from mod_swstate import sw_dens

def calcN2_1D(T,S,ZZ,lat0,info=1): 
  """
    #  Calculate N2 = -g/rho*d2(rho)/dz2 
    #  Using in situ T and S - calculate rho relative to the
    #  mid-grid depth - following Chelton, 1996
    # Use the neutral density gradient method, Chelton et al., 1996
    Input: T,S,ZZ - observed T/S and depths
           lat0 - latitude
           info = 1 - print depth, computed N2
                = 0 - no print out 
  """
  print('Computing N2 at mid-grid depths, {0} levels'.format(ZZ.shape[0]))
  kdm = T.shape[0]
#print('Calculating pressure at midpoints ...')
  grav  = 9.81
  rho0  = 1025.0  
  Z_phi = np.zeros(kdm-1)  # depths for e/function Phi
  N2    = np.zeros(kdm-1)
  n2fill= 1.e-8     # missing values
  for kk in range(kdm-1):
#    print(' Layer {0}'.format(kk))
    z1 = ZZ[kk]
    z2 = ZZ[kk+1]
    t1 = T[kk]
    t2 = T[kk+1]
    s1 = S[kk]
    s2 = S[kk+1]
#    Z1=z1*np.ones((jdm,idm))
#    Z2=z2*np.ones((jdm,idm))
    p1_db, p1_pa = sw_press(z1,lat0)  # pressure upper interface
    p2_db, p2_pa = sw_press(z2,lat0)  # pressure bottom interface
#
# Find mid-point of the layers
    p0_db     = 0.5*(p1_db+p2_db)
    z0        = 0.5*(z1+z2)
    Z_phi[kk] = z0
#
# Calculate rho(z1--->z0) with depth reference at midpoint
# and rho(z2--->z0) at midpoint 
    t_z1z0   = sw_ptmp(s1,t1,p1_db,p0_db)
    t_z2z0   = sw_ptmp(s2,t2,p2_db,p0_db)  
    rho_z1z0 = sw_dens(s1,t_z1z0,p0_db)
    rho_z2z0 = sw_dens(s2,t_z2z0,p0_db)

# Calculate d(rho)/dz for z0 - center-difference
    drho_dz = (rho_z1z0 - rho_z2z0)/(z1 - z2)
    N2z0    = -grav/rho0*drho_dz 
    N2[kk]  = N2z0
    if info == 1:
      print('k={0}, z={1}, N2: {2}'.format(kk,z0,N2z0))

# If N2 < 0 - density inversion happens in some ~homogeneous layers
# when parcels are brought down from z1 and up from z2
# replace with above N2 or below of surface layer
  N2[np.where(N2==0.)] = n2fill
  II = np.where(N2<0.)[0]
  JJ = np.where(N2>0.)[0]
  if II.size>0: 
    print('Fixing N2<0, {0} occurrences'.format(II.shape[0]))
    for i in range(II.shape[0]):
      kk = II[i]
      if kk > min(JJ):
        j = max(np.where(JJ<kk)[0])
        krpl = JJ[j]
        N2[kk] = N2[krpl]
      else:
        j = min(np.where(JJ>kk)[0])
        krpl = JJ[j]
        N2[kk] = N2[krpl]

#      N2z0 = n2fill
#      N2[kk]=N2z0


  return N2, Z_phi

def intrpN2_bottom(N2z,Z_phi,ZZ,zbtm):
# most of the profiles are much shallower than the bottom
# Following Chelton et al., interpolate N2 to the bottom
# assuming N2 ~0 at the bottom
#
# ZZ - layer interface depths - need to update adding new depths
# Z_phi - depths where N2 is computed
# N2z - N2 at z-phi points
#
  kN = Z_phi.shape[0]
  zE = Z_phi[-1]   # last z-phi point
  zzE = ZZ[-1]     # last interface depth

# Depths array:
  dZ = 200.
  ZI = np.arange(0.,11000.,dZ)
  ZI = -ZI

  dmm = abs(ZI-zzE)
  kS = np.argmin(dmm)
  if abs(ZI[kS]) < abs(zzE):
    kS = kS+1

  if abs(zE-zzE) < dZ/2. or abs(ZI[kS]) >= abs(zbtm):
    print(' No need interp N2 to bottom, last intrf Z {0:.1f} btm={1:.1f}'.\
          format(zzE, zbtm))
    return N2z, Z_phi, ZZ

  dmm = abs(ZI-zbtm)
  kE = np.argmin(dmm)
  if abs(ZI[kE]) >= abs(zbtm):
    kE = kE-1                  # keep last value above the bottom


  ZZadd = ZI[kS:kE+1]
  ZZadd = np.append(ZZadd,zbtm)
  ZZadd = np.insert(ZZadd,0,zzE)
  nZ = ZZadd.shape[0]

# Linear interpolation:
  N0 = N2z[-1]
  N1 = 1.e-18   
  Nint = np.zeros((nZ-1))
  Zint = np.zeros((nZ-1))
  cc = -1
  for kk in range(nZ-1):
    zz1 = ZZadd[kk]
    zz2 = ZZadd[kk+1]
    zz = 0.5*(zz2+zz1)
    cc += 1
    Nint[cc] = N0*(zz-zbtm)/(zE-zbtm)+N1*(zz-zE)/(zbtm-zE)
    Zint[cc] = zz

# Add to original N2 profile:
  N_new = np.concatenate((N2z,Nint))
  Z_phi_new = np.concatenate((Z_phi,Zint))
  ZZ_new = np.concatenate((ZZ,ZZadd[1:])) 

  return N_new, Z_phi_new, ZZ_new

#
# =======================================================
#
# Numerically Solve Sturm-Liouville e/value problem 
# The e/value problem is solved for A*Phi = mu*Phi
# to find e/values mu(i), i=1, ..., m
# No depth control/ QC of N2 is done here
# 
# =======================================================
def solve_SturmLiouv_1D(N2z,Z_phi,ZZ,zbtm_ETOPO,lat0):
# Solve A*phi = mu*phi
# 
# N2z - N2 at Z_phi points- mid depths of grid cells
# ZZ - depths of the grid-cell interfaces
#
#  import mod_solver
#  importlib.reload(mod_solver)
  from mod_solver import runmn
  from mod_solver import form_mtrxA
  from mod_solver import eig_unfoldA
  from mod_solver import eig_wkb

  print('Solving Sturm-Liouville ...')
  omg = 7.29e-5 # Earth angular velocity

# Find "bottom": here profiles are not to the actualy bottom
# Depending to observ. casts, these typically cover only 
# some fraction of the water column
# Assumed that last depth = last observed T/S  
# Escape coastal regions < 20 m
#    k1 = np.where(np.isnan(N2z))[0]
#
# Fill below bottom:
# Assume that actual bottom is deeper than deepest observed value
# extend last grid cell to the actual bottom
  kbtm = N2z.shape[0]-1
  zbtm = zbtm_ETOPO
#
# Filter N-point running mean
  N2zF = runmn(N2z,ZZ,mnwnd=3)

# Test - plot N profiles
# Convert to cycles/hr
  f_testN = 0
  if f_testN == 1:
    N_hrz = np.sqrt(N2z)
    N_chr = 3600.*N_hrz   # cycles per hr
    NF_chr = 3600.*np.sqrt(N2zF)
    plt.ion()
    fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
    plt.plot(N_chr,Z_phi)  
    plt.plot(NF_chr,Z_phi)


# Create Matrix A with Dk, Dk+1 for 2nd derivative of N2
  AA = form_mtrxA(Z_phi,kbtm,zbtm)
#
# Form 1/N2*AA:
# The matrix is cutoff by the bottom depth
  ka = AA.shape[0]
  na = AA.shape[1]
  N2A = np.zeros((ka,na))
  for kk in range(ka):
    n2 = N2zF[kk]  # filtered profile
    if n2 == 0:
      n2 = 1.e-12
    N2A[kk,:] = 1./n2*AA[kk,:]

#
# Use python eig function, for this
# form Matrix A unfolding
# the 3 -elemnts form and find eigenvalues/vectors
# W - eigenvalues, not sorted our by magnitude!
# V- corresponding eigenvectors (in columns)
  W, V = eig_unfoldA(kbtm,N2A)

# Calculate
# Choose lmbd1 = as min(abs(W)) and lmbd1<0
  W[np.where(W>0.)] = np.nan  # ignore e/values that result complex Radius
  im = np.where(np.abs(W) == np.nanmin(np.abs(W)))[0][0]
#  im = np.where((W < 0.) & (np.abs(W) == np.min(np.abs(W))))[0][0]
  if W[im] > 0.:
    print('WARNING: W[im] >0: im={0}, W[im]={1}, zbtm={4:6.1f}, ii={2}, jj={3}'.\
           format(im, W[im], ii, jj, zbtm))

  Rearth = 6371.e3  # Earth R
  fcor = 2.*omg*np.sin(np.deg2rad(lat0))
  betta = (2.*omg*np.cos(np.deg2rad(lat0)))/Rearth

  if abs(lat0) >= 5.0: 
    lmbd = np.sqrt(-1./(W[im]*fcor**2))
  else:
    lmbd = (-1./(4.*(betta**2)*(W[im])))**(1./4.)

# Baroclinic gravity wave Phase speed c(mode m):
  Cph1 = np.sqrt(-1/W[im])
  Phi1 = V[:,im]          # eigenvector

  RsbNum = lmbd*1.e-3  # km

  return RsbNum, Cph1, Phi1 



