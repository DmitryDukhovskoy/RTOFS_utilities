# Convert binary cice4 to binary netcdf cice6 restart
# cice restart file is from GOFS3.1 GLBb0.08 - expt 93.0
#
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime

# CICE4 restart date:
YRc4     = 2020
MMc4     = 1
MDc4     = 1
HRc4     = 9
pthrst4  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/GLBb0.08_expt93.0/'
cicerst4 = 'cice.restart.{0}{1:02d}{2:02d}{3:02d}'.format(YRc4,MMc4,MDc4,HRc4)

# Restart CICE6 template:
YRtmp       = 2000
MMtmp       = 12
MDtmp       = 16
HRtmp       = 0
pthrstT     = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/MOM6_CICE6/' +\
              'cice_restart/'
cicerstT    = 'iced.{0}-{1:02d}-{2:02d}-{3:02d}tmp.nc'.\
               format(YRtmp,MMtmp,MDtmp,HRtmp)

# Output restart:
# Desired date/time
YRnew   = 2020
MMnew   = 1
MDnew   = 1
HRnew   = 0
pthrst6     = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/MOM6_CICE6/' +\
              'cice_restart/'
cicerst6    = 'iced.{0}-{1:02d}-{2:02d}-{3:02d}new.nc'.\
               format(YRnew,MMnew,MDnew,HRnew)


# File names with paths:
fl_restart4 = pthrst4 + cicerst4
fl_restartT = pthrstT + cicerstT
fl_restart6 = pthrst6 + cicerst6

print('Input CICE4 restart:    ' + fl_restart4)
print('Template CICE6 restart: ' + fl_restartT)
print('Output CICE6 restart:   ' + fl_restart6)


import mod_datm_utils as mdatm
importlib.reload(mdatm)

# Create new restart from template for writing CICE4 fields
mdatm.cice6_newfile(fl_restartT, fl_restart6)

# Grid CICE4 - unformatted binary file
pthgrd4 = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_paraCd/fix/'
grdfl4  = 'regional.cice.r'
fgrdin4 = pthgrd4 + grdfl4

# Grid, use CICE6 grid - should be similar to CICE4
pthgrd  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6/'
grdfl   = 'grid_cice_NEMS_mx008.nc'
fgrdin  = pthgrd + grdfl

pthdpth = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/008mom6cice6/INPUT/'
dpthfl  = 'depth_GLBb0.08_09m11ob2_mom6.nc'
fdpthin = pthdpth + dpthfl

print('Creating CICE6 restart for {0}/{1:02d}/{2:02d} {3:02d}hr UTC'.\
      format(YRnew, MMnew, MDnew, int(HRnew)))
print('CICE4 restart:  ' + pthrst4 + cicerst4)
print('CICE6 template: ' + pthrst6 + cicerst6)
print('New restart will modify the template ' + cicerst6)
print(' =================================== \n')


# CICE parameters:
puny      = 1.e-11
c0        = 0.0
c1        = 1.0
c2        = 2.0
p5        = 0.5
Lsub      = 2.835e6    # latent heat sublimation fw (J/kg)
Lvap      = 2.501e6    # latent heat vaporization fw (J/kg)
Lfresh    = Lsub - Lvap # latent heat of melting of fresh ice (J/kg)
cp_ice    = 2106.       # specific heat of fresh ice (J/ kg/K)
rhos      = 330.        # density of snow (kg/m3)
hs_min    = 1.e-4       # min snow thickness for computing Tsno (m)
nsal      = 0.407
msal      = 0.573
min_salin = 0.1      # threshold for brine pocket treatment
saltmax   = 3.2        # max S at ice base
hg        = 1.e20    # bad values, land mask, etc.
rdn2dgr   = 180./np.pi
dgr2rdn   = np.pi/180.


# ====================================
# 
#  Input parameters - check with ice_in
#  Edit mod_cice6_utils
# =====================================
# CICE4 params and restart fields
# Note in GOFS3.1 CICE has 1 row less than HYCOM
# ny = 3297 (and it is 3298 in HYCOM)
#cice4 = mc6util.param_cice4()
cice4 = mc6util.cice4(nx=4500, ny=3297)


# In CICE4, can be isosalaine = const S or
# constant shape S
salin = np.zeros((cice4.nilyr+1))
if saltmax > min_salin:
  l_brine = True

  for k in range(cice4.nilyr):
    zn = (float(k+1)-0.5)/(float(cice4.nilyr))
    salin[k] = (saltmax/2.)*(1.-np.cos(np.pi*zn**(nsal/(msal+zn))))
  salin[k+1] = saltmax

else:
  l_brine = False
  salin = 0.0

# Read restart  
#  Read CICE4 restart fields
#  unformatted binary big endian

# ----------------------------
def read_rest_cice4(fid):
  """
    Read 2D record from open file
    secuential binary file
  """
  fdump = np.fromfile(fid, dtype='>i4', count=1)[0]
  A     = np.fromfile(fid, dtype='>f8', count=nx*ny)
  A     = np.reshape(A,(ny,nx), order='C')
  fdump = np.fromfile(fid, dtype='>i4', count=1)[0]
  return A

def print_minmax(sfld,A):
  print('   {2} min/max:  {0}/{1}'.format(np.min(A),np.max(A),sfld))
  return

# ----------------------------
spval = 1.e30

try:
  fid = open(fl_restart4, 'rb')
except:
  print('Could not open '+fl_restart4)

print('Reading restart: ' + fl_restart4)

fid.seek(0)
# Read Fortran binary
fdump   = np.fromfile(fid, dtype='>i4', count=1)[0]
istep   = np.fromfile(fid, dtype='>i4', count=1)[0]
runtime = np.fromfile(fid, dtype='>f8', count=1)[0]  # total elapsed time, sec
frtime  = np.fromfile(fid, dtype='>f8', count=1)[0]  # forcing time, sec
fdump   = np.fromfile(fid, dtype='>i4', count=1)[0]

print('Restart: step={0}, total time(yrs)={1}, forcing last update(hrs)={2}'\
       .format(istep,runtime/(3600*24*365.25),frtime/3600.))

# Read state variables:
# Tsfc is the only tracer read in this file
nx     = cice4.nx
ny     = cice4.ny
ncat   = cice4.ncat
ntilyr = cice4.ntilyr  # total # of icelrs * cat 
ntslyr = cice4.ntslyr

aicen = np.zeros((ncat,ny,nx), dtype='float64')
vicen = np.zeros((ncat,ny,nx), dtype='float64')
vsnon = np.zeros((ncat,ny,nx), dtype='float64')
trcrn = np.zeros((ncat,ny,nx), dtype='float64')

for n in range(ncat):
  print(' Category {0}'.format(n+1))
# Read ice area for category n
  A = read_rest_cice4(fid)
  aicen[n,:,:] = A
  print_minmax('ice area',A)

# Read ice volume/per m2/ for category n
  A = read_rest_cice4(fid)
  vicen[n,:,:] = A
  print_minmax('ice vol',A)

# Read snow volume/m2 for category n
  A = read_rest_cice4(fid)
  vsnon[n,:,:] = A
  print_minmax('snow vol',A)

# Read tracer 1 (surf T = Tsfcn by categories) - only 1 tracer in CICE4
  A = read_rest_cice4(fid)
  trcrn[n,:,:] = A
  print_minmax('surf T',A)

# Ice energy (J/m2):
eicen = np.zeros((ntilyr,ny,nx), dtype='float64')
print('\n Ice energy:')

for k in range(ntilyr):
  A = read_rest_cice4(fid)
  eicen[k,:,:] = A
  print_minmax('{0} eicen'.format(k+1),A)

# Snow energy:
esnon = np.zeros((ntslyr,ny,nx), dtype='float64')
print('\n Snow energy:')

for k in range(ntslyr):
  A = read_rest_cice4(fid)
  esnon[k,:,:] = A
  print_minmax('{0} esnon'.format(k+1),A)

# Velocities:
print('\n Velocity components:')
A = read_rest_cice4(fid)
uvel = A.copy()
print_minmax('U vel',A)

A = read_rest_cice4(fid)
vvel = A.copy()
print_minmax('V vel',A)

# Radiation fields
# 4 radiative categories
# for calculating albedo for visible and IR wavelengths
# and penetrating sh/wave

# Scale factor to change MKS units
# for shortwave components
# default = 1
print('\n Radiation fields, W/m2: ')
A = read_rest_cice4(fid)
scale_factor = A.copy()
print_minmax('Scale Factor',A)

A = read_rest_cice4(fid)
swvdr = A.copy()
print_minmax('Sh/wave down vis. direct',A)

A = read_rest_cice4(fid)
swvdf = A.copy()
print_minmax('Sh/wave down vis. diff',A)

A = read_rest_cice4(fid)
swidr = A.copy()
print_minmax('Sh/wave down near IR dir',A)

A = read_rest_cice4(fid)
swidf = A.copy()
print_minmax('Sh/wave down near IR diff',A)

# Ocean stress, N/m2
print('\n Ocean stress components, N/m2:')

A = read_rest_cice4(fid)
strocnxT = A.copy()
print_minmax('ocean stress x-comp',A)

A = read_rest_cice4(fid)
strocnyT = A.copy()
print_minmax('ocean stress y-comp',A)

# Internal stress, stress tensor kg/s2
# (1) northeast, (2) northwest, (3) southwest, (4) southeast
print('\n Internal stress, kg/s2: ')

A = read_rest_cice4(fid)
stressp_1 = A.copy()
print_minmax('stressp_1',A)

A = read_rest_cice4(fid)
stressp_3 = A.copy()
print_minmax('stressp_3',A)

A = read_rest_cice4(fid)
stressp_2 = A.copy()
print_minmax('stressp_2',A)

A = read_rest_cice4(fid)
stressp_4 = A.copy()
print_minmax('stressp_4',A)

A = read_rest_cice4(fid)
stressm_1 = A.copy()
print_minmax('stressm_1',A)

A = read_rest_cice4(fid)
stressm_3 = A.copy()
print_minmax('stressm_3',A)

A = read_rest_cice4(fid)
stressm_2 = A.copy()
print_minmax('stressm_2',A)

A = read_rest_cice4(fid)
stressm_4 = A.copy()
print_minmax('stressm_4',A)

A = read_rest_cice4(fid)
stress12_1 = A.copy()
print_minmax('stress12_1',A)

A = read_rest_cice4(fid)
stress12_3 = A.copy()
print_minmax('stress12_3',A)

A = read_rest_cice4(fid)
stress12_2 = A.copy()
print_minmax('stress12_2',A)

A = read_rest_cice4(fid)
stress12_4 = A.copy()
print_minmax('stress12_4',A)

# Ice mask for dynamics
print('\n Ice Mask for Dynamics: ')

A = read_rest_cice4(fid)
iceumask = A.copy()
print_minmax('iceumask',A)

# For trully coupled HYCOM-CICE these fields
# are not needed
# if defined ocean mixed layer in CICE
# This is for ocean mixed layer defined (for GOFS )
print('\n Ocean mixed layer: \n');

A = read_rest_cice4(fid)
sst = A.copy()
print_minmax('sst',A)

A = read_rest_cice4(fid)
frzmlt = A.copy()
print_minmax('frzmlt',A)

fid.close()

# Mask out land points:
print(' Masking out fields ')
aicen        = np.where(aicen > 0.5*spval, 0., aicen)
vicen        = np.where(vicen > 0.5*spval, 0., vicen)
vsnon        = np.where(vsnon > 0.5*spval, 0., vsnon)
trcrn        = np.where(trcrn > 0.5*spval, 0., trcrn)
eicen        = np.where(eicen > 0.5*spval, 0., eicen)
esnon        = np.where(esnon > 0.5*spval, 0., esnon)
uvel         = np.where(uvel > 0.5*spval, 0., uvel)
vvel         = np.where(vvel > 0.5*spval, 0., vvel)
scale_factor = np.where(scale_factor > 0.5*spval, 0., scale_factor)
swvdr        = np.where(swvdr > 0.5*spval, 0., swvdr)
swvdf        = np.where(swvdf > 0.5*spval, 0., swvdf)
swidr        = np.where(swidr > 0.5*spval, 0., swidr)
swidf        = np.where(swidf > 0.5*spval, 0., swidf)
strocnxT     = np.where(strocnxT > 0.5*spval, 0., strocnxT)
strocnyT     = np.where(strocnyT > 0.5*spval, 0., strocnyT)
stressp_1    = np.where(stressp_1 > 0.5*spval, 0., stressp_1)
stressp_3    = np.where(stressp_3 > 0.5*spval, 0., stressp_3)
stressp_2    = np.where(stressp_2 > 0.5*spval, 0., stressp_2)
stressp_4    = np.where(stressp_4 > 0.5*spval, 0., stressp_4)
stressm_1    = np.where(stressm_1 > 0.5*spval, 0., stressm_1)
stressm_3    = np.where(stressm_3 > 0.5*spval, 0., stressm_3)
stressm_2    = np.where(stressm_2 > 0.5*spval, 0., stressm_2)
stressm_4    = np.where(stressm_4 > 0.5*spval, 0., stressm_4)
stress12_1   = np.where(stress12_1 > 0.5*spval, 0., stress12_1)
stress12_3   = np.where(stress12_3 > 0.5*spval, 0., stress12_3)
stress12_2   = np.where(stress12_2 > 0.5*spval, 0., stress12_2)
stress12_4   = np.where(stress12_4 > 0.5*spval, 0., stress12_4)
sst          = np.where(sst > 0.5*spval, 0., sst)
frzmlt       = np.where(frzmlt > 0.5*spval, 0., frzmlt)


# ----------------------
# See: ufs-weather-model/CICE-interface/CICE/cicecore/cicedynB/infrastructure/io/io_netcdf
# ice_restart.F90
#
Lmsk  = mc6util.read_ncfile(fdpthin,'wet')
cice6 = mc6util.cice6()

f_chck = False
if f_chck:
  ncdata = ncFile(fl_restart6,'r')
  coszen = ncdata['coszen'][:]

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im1 = ax1.pcolormesh(coszen)
  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2)
  im1.set_clim(0,1)

  stl = 'coszen, {0}/{1}/{2} {3}:{4}'.format(YRnew,MMnew,MDnew,HRnew,0)
  ax1.set_title(stl)

btx = 'cice4_cice6_restart.py'

f_plotPolar = False
if f_plotPolar:
  ulati6 = mc6util.read_ncfile(fgrdin, 'ulat')
  uloni6 = mc6util.read_ncfile(fgrdin, 'ulon')
  LAT, LON = mc6util.grid_rad2dgr(ulati6, uloni6)
#  Umsk = mc6util.read_ncfile(fl_restart6,'iceumask')
#  fldplt = 'fsnow'  # snow fall rate ? kg/m2/sec
#  fsnow = mc6util.read_ncfile(fl_restart6,'fsnow')
  fldplt = 'sice001'
  icat = 5
  A3 = mc6util.read_ncfile(fl_restart6,fldplt)
  A2D = A3[icat-1,:,:]
#  A2D = mc6util.read_ncfile(fl_restart6,fldplt)
  A2D = np.where(Lmsk==0, np.nan, A2D)
  rmin = 0.0
  rmax = 3.
  stl = 'CICE6, {0}/{1}/{2}, {3}, icat={4} '.format(YRtmp,MMtmp,MDtmp,fldplt,icat)
  mc6util.plot_polar_2D(LON, LAT, A2D, region='Arctic',  \
                  rmin=rmin, rmax=rmax, stl=stl)
  bottom_text(btx)

# Plot area-mean snow thickness:
#  fldplt = 'hsnow_mean, m'
#  hsmn = np.sum(vsnon, axis=0)  # aggregated snow vol (m3/m2) = mean snow thickness
#  A2D  = hsmn
#  A2D = np.where(Lmsk==0, np.nan, A2D)
#  rmin = 0.0
#  rmax = 0.5
#  stl = 'CICE6, {0}/{1}/{2}, {3} '.format(YRtmp,MMtmp,MDtmp,fldplt)
#  rmin = -1.8
#  rmax = 10.
#  mc6util.plot_polar_2D(LON, LAT, A2D, region='Antarctic',  \
#                  rmin=rmin, rmax=rmax, stl=stl)
#  bottom_text(btx)
#------------------------------
  
# Edit/ create missing fields:
# coszen - cosine of solar zenith angle
#          negative for sun below horizon
# zen angle = 90 - solar inclination(time, location)
# Here, a simple calculation of zenith angle is performed
# To match CICE6 restart fields, CICE4 does not have coszen
# For accurate calculation: see icepack/icepack_orbital.F90
#
# Check grid info
# lon/lat of U-points, here just check to make sure
# that both grids are Arakawa B grid, if not - need to add
# module to relocate and remap U/V components
#
# Lon/lat in CICE4:
ulati4 = mc6util.read_cice4_grid(fgrdin4, 'ulati', IDM=nx, JDM=ny)
uloni4 = mc6util.read_cice4_grid(fgrdin4, 'uloni', IDM=nx, JDM=ny)

# Read lon/lat from CICE6 restart template:
ulati6 = mc6util.read_ncfile(fgrdin, 'ulat')
uloni6 = mc6util.read_ncfile(fgrdin, 'ulon')

# It seems to be an error in CICE4 regional grid:
# last column (in the top ~100 rows) in lat is repeated (end-1) column
# in CICE6, these are different columns
mc6util.check_cice_grids(ulati4,uloni4,ulati6,uloni6)

# Compute coszen for CICE4 restart data
# adjusting for new restart date
# UTC --> time_zone=0
# ==> TO DO: Add land mask = 0 to match CICE6 restart 
DVnew      = [int(YRnew),int(MMnew),int(MDnew),int(HRnew)]
dnmb_new   = mtime.datenum(DVnew)
coszen_new = mc6util.compute_coszen(ulati6, uloni6, dnmb_new, time_zone=0)
#mc6util.add_fld2D_ncfile(fl_restart6, 'coszen')

# Snow fall rate - unknown, make 0:
fsnow = np.zeros((ny,nx))

# Volume-weighted ice age, iage:
# see: icepack_age.F90, icepack_therm_vertical.F90
# cannot be reconstructed from CICE4:
iage = aicen*0.0

# alvl - level (undeformed) ice area fraction by categories
# alvl cannot be recovered from CICE4 restart
# assume no ridges only level ice
alvl = np.where(aicen > 0.0, 1., 0.)

# vlvl - level ice volume fraction of vicen
# cannot be recovered from CICE4
vlvl = np.where(vicen > 0.0, 1., 0.)

# apnd - pond area tracer, ponded level-ice fraction
# unponded ice, cat=i: (1-apnd)*alvl*aice
# ponded ice: apnd*alvl*aice, where aice - ice area fraction
# cannot be recovered from CICE4
# assume no ponds
apnd = iage

# hpnd - melt pond depth, m
# cannot be recovered from CICE4
# hpnd =0
hpnd = iage

# ipnd - melt pond refrozen lid thickness
# ipnd cannot be rcovered from CICE4
ipnd = iage 

# dhs - change in snow thickness
# icepack_therm_vertical.F90: 
# dhs = hsn - works - hsn_new, where hsn_new - thickness of new snow (snowfall)
# Cannot be recovered for CICE4
dhs = iage

# ffrac - fraction of fsurfn (atm-ice surface heat flux  (W/m2)) 
#         over pond used to melt ipond (pond ice)
# Cannot be recovered from CICE4
ffrac = iage

# Internal ice energy in CICE 4 = 5 cat x 4 lrs, J/m2
# in CICE6 - ice enthalpy by layers J/m3
# ncat, ny, nx
# convert eicen, esnon in CICE4 to qicen, qsnon in CICE6  - energy by volume J/m3
#
# eicen, the internal ice energy in layer k
# eicen(k) =  vicen(k)/Ni * qicen (ice layer enthalpy), Ni - # of ice layers
# CICE4 to CICE6 conversion:
# qicen [J/m3] = eicen[J/m2] * Ni/vicen [m]
#
# In CICE6 enthalpies are carried in tracer array:
#  trcrn(i,j,nt_qsno+k-1,n) , k = 1,nslyr, n=1,..,5 cat - snow enth., J/m3
#
# see: https://github.com/NOAA-EMC/CICE/blob/5840cd1931e2e32b9dfded0c19049d0f1ec3d04c/configuration/tools/cice4_restart_conversion/convert_restarts.f90
#
# Get Combined ice layers by categories for 
# eicen([1,...,ncat*nilyr],j,i) ilyr1=[4, 8,..,20] - start-end indices of 
# eicen for each ice category
# and esnon(slyr1,j,i) slyr1=[1,...,5]
#
# The process has 2 steps:
# (1) convert esnon (snow/ice int. energy)(J/m2) --> qsnon 
#       (snow/ice layer enthalpy, J/m3)
#     and eice (J/m2) --> qice (J/m3)
#     keep snow/ice on cice4 vertical layers
# (2) if # of layers in CICE6 is different, remap 
#     qsnon and qice onto N layers of CICE6
#
nilyr = cice4.nilyr
nslyr = cice4.nslyr
ilyr1 = np.zeros((ncat), dtype='int')
slyr1 = np.zeros((ncat), dtype='int')
ilyr1[0] = 0   # 1st index, python
slyr1[0] = 0
for n in range(1,ncat):
  ilyr1[n] = ilyr1[n-1] + nilyr
  slyr1[n] = slyr1[n-1] + nslyr

# Surface T, enthalpy and maxT:
# see code ufs-weather-model/CICE-interface/CICE/icepack/columnphysics/
# icepack_therm_vertical.F90: L 752
#       !-----------------------------------------------------------------
#      ! Tmax based on the idea that dT ~ dq / (rhos*cp_ice)
#      !                             dq ~ q dv / v
#      !                             dv ~ puny = eps11
#      ! where 'd' denotes an error due to roundoff.
#      !-----------------------------------------------------------------
#
# Note qsnon, qice < 0 !
Tmin = -100.   # minimum snow T
qsnon = np.zeros((nslyr,ncat,ny,nx))
qicen = np.zeros((nilyr,ncat,ny,nx))
for n in range(ncat):                 # ice categories
  for k in range(nslyr):              # snow layers
    vsn    = vsnon[n,:,:]               # snow volume per m2 in cat=n
    ai_cat = aicen[n,:,:]
    ai_cat = np.where(ai_cat < 1.e-20, 1.e-20, ai_cat)
    hsn    = vsn/ai_cat
    vsn    = np.where(vsn < 1.e-20, 1.e-20, vsn)

    iesnon = slyr1[n]+k
    qsn    = esnon[iesnon,:,:]*float(nslyr)/vsn  # energy (J/m2) --> enthalpy (J/m3)
    qsn    = np.where(qsn > -rhos*Lfresh, -rhos*Lfresh, qsn)
    qsn    = np.where(ai_cat <= puny, 0., qsn)  # no ice
    qsn    = np.where(hsn <= hs_min/(float(nslyr)), 0., qsn)   # no snow
# Units of qsn = J/m3, puny = tiny vol (m3/m2):
    Tmax   = -qsn*puny*float(nslyr)/(rhos*cp_ice*vsn)  # Large T for small vsnon

    vsn    = np.where(hsn <= hs_min/(float(nslyr)), 0., vsn)
    Tmax   = np.where(vsn < 1.e-11, 0., Tmax)

#
# Snow Temperature for checking qsnon:
# in icepack_itd code trcrn is slice from 5D trcrn:
# see cicecore/cicedynB/infrastructure/ice_restart_driver.F90
# 659       do j = 1, ny_block
# 660       do i = 1, nx_block
# 661          if (tmask(i,j,iblk)) &
# 662             call icepack_aggregate(ncat  = ncat,                  &
# 663                                    aicen = aicen(i,j,:,iblk),     &
# 664                                    trcrn = trcrn(i,j,:,:,iblk),   &
    qT0 = -Lfresh*rhos
    qsn = np.where(hsn <= hs_min/(float(nslyr)), qT0, qsn)
    zTsn = (Lfresh + qsn/rhos)/cp_ice
# Check max snow T from max enthalpy, should be ~0:
    JJ,II = np.where(zTsn > 1.e-11)
    if len(JJ) > 0:
      qsn[JJ,II] = qT0 

# Zap entire snow volume if T is out of bounds:
    print('Snow enthalpy: lr = {0}, Tmax = {1}'.format(k+1,np.max(Tmax)))
    print('cat={0} slyr={1} min/max snow T: {2}/{3}'.\
         format(n, k+1, np.min(zTsn), np.max(zTsn)))
    J1,I1 = np.where(zTsn < Tmin)
    J2,I2 = np.where(zTsn > Tmax)
    if len(J1) > 0 or len(J2) > 0:
      print('Tsnow is out of bound zeroing snow volume ')
      vsn = np.where(zTsn < Tmin, 0., vsn)
      qsn = np.where(zTsn < Tmin, 0., qsn)
      vsn = np.where(zTsn > Tmax, 0., vsn)
      qsn = np.where(zTsn > Tmax, 0., qsn)

    qsnon[k,n,:,:] = qsn
    vsnon[n,:,:]   = vsn 

# sea ice enthalpy 
# qice < 0
# Convert energy CICE4 J/m2 ---> enthalpy CICE6 J/m3
# on the same vertical ice layers 
  for k in range(nilyr):
    ai_cat = aicen[n,:,:]
    iicen  = ilyr1[n]+k
    print('eicen index: iicen = {0}'.format(iicen))
    vin    = vicen[n,:,:]               # ice volume per m2 in cat=n
    vin    = np.where(vin < 1.e-30, 1.e-30, vin)
    qin    = eicen[iicen,:,:]*float(nilyr)/vin  # J/m2 --> J/m3
    qin    = np.where(ai_cat <= puny, 0.0, qin)
    qicen[k,n,:,:] = qin


# Interpolate from nilyr=4 in CICE4 to nilyr=7 ice layers in CICE6
if not cice6.nilyr == cice4.nilyr:
  qicen = mc6util.remap_enthalpy_bins(qicen, cice4.nilyr, cice6.nilyr)

# Snow enthaply interpolation has not been used
# but probably should work fine
if not cice6.nslyr == cice4.nslyr:
  print('!!! Need to check snow enthalpy interpolation \n\n!!!!')
  qsnon = mc6util.remap_enthalpy_bins(qsnon, cice4.nilyr, cice6.nilyr)

#AA = PP

# =======================================================
# 
#  scale_factor - scaling factor for shortwave radiation components
# use from CICE6 template? CICE4 scale_factor 
# is different
# scale_factor: netsw scaling factor (new netsw / old netsw)
# see: icepack_shortwave.F90
#
print(' \n\n -------------\n Creating CICE6 restart')
print('Created CICE6 restart: \n' + fl_restart6)

if not os.path.isfile(fl_restart6):
  raise Exception ('CICE6 restart missing: ' + fl_restart6)

mc6util.modify_fld_nc(fl_restart6,'uvel',uvel)
mc6util.modify_fld_nc(fl_restart6,'vvel',vvel)
mc6util.modify_fld_nc(fl_restart6,'scale_factor',scale_factor)
mc6util.modify_fld_nc(fl_restart6,'swvdr',swvdr)
mc6util.modify_fld_nc(fl_restart6,'swvdf',swvdf)
mc6util.modify_fld_nc(fl_restart6,'swidr',swidr)
mc6util.modify_fld_nc(fl_restart6,'swidf',swidf)
mc6util.modify_fld_nc(fl_restart6,'strocnxT',strocnxT)
mc6util.modify_fld_nc(fl_restart6,'strocnyT',strocnyT)
mc6util.modify_fld_nc(fl_restart6,'stressp_1',stressp_1)
mc6util.modify_fld_nc(fl_restart6,'stressp_2',stressp_2)
mc6util.modify_fld_nc(fl_restart6,'stressp_3',stressp_3)
mc6util.modify_fld_nc(fl_restart6,'stressp_4',stressp_4)
mc6util.modify_fld_nc(fl_restart6,'stressm_1',stressm_1)
mc6util.modify_fld_nc(fl_restart6,'stressm_2',stressm_2)
mc6util.modify_fld_nc(fl_restart6,'stressm_3',stressm_3)
mc6util.modify_fld_nc(fl_restart6,'stressm_4',stressm_4)
mc6util.modify_fld_nc(fl_restart6,'stress12_1',stress12_1)
mc6util.modify_fld_nc(fl_restart6,'stress12_2',stress12_2)
mc6util.modify_fld_nc(fl_restart6,'stress12_3',stress12_3)
mc6util.modify_fld_nc(fl_restart6,'stress12_4',stress12_4)
mc6util.modify_fld_nc(fl_restart6,'iceumask',iceumask)
mc6util.modify_fld_nc(fl_restart6,'fsnow',fsnow)
mc6util.modify_fld_nc(fl_restart6,'aicen',aicen)
mc6util.modify_fld_nc(fl_restart6,'vicen',vicen)
mc6util.modify_fld_nc(fl_restart6,'vsnon',vsnon)
mc6util.modify_fld_nc(fl_restart6,'Tsfcn',trcrn)
mc6util.modify_fld_nc(fl_restart6,'coszen',coszen_new)
mc6util.modify_fld_nc(fl_restart6,'iage',iage)
mc6util.modify_fld_nc(fl_restart6,'alvl',alvl)
mc6util.modify_fld_nc(fl_restart6,'vlvl',vlvl)
mc6util.modify_fld_nc(fl_restart6,'apnd',apnd)
mc6util.modify_fld_nc(fl_restart6,'hpnd',hpnd)
mc6util.modify_fld_nc(fl_restart6,'ipnd',ipnd)
mc6util.modify_fld_nc(fl_restart6,'dhs',dhs)
mc6util.modify_fld_nc(fl_restart6,'ffrac',ffrac)

#  4D fields:
# sice - ice bulk salinity
# sice - 4D field, written by layers as 3D (ncat,nj,ni)
# ufs-weather-model/CICE-interface/CICE/cicecore/cicedynB/infrastructure/io/io_netcdf
# ice_restart.F90
# ufs-weather-model/CICE-interface/CICE/cicecore/shared/ice_restart_column.F90
#aice = np.sum(aicen, axis=0)
#
# Ice salinity by layers - compute S profile using BZ99 formulation:
for ik in range(1,cice6.nilyr+1): 
  sice_lr = mc6util.sice_lr_cice4(ik, cice6.nilyr, aicen)
  fldout = 'sice{0:03d}'.format(ik)
  print('Updating ' + fldout)
  mc6util.modify_fld_nc(fl_restart6,fldout,sice_lr)

# Ice enthalpy by layers:
for ik in range(1,cice6.nilyr+1): 
  qice_lr = qicen[ik-1,:,:,:]
  fldout = 'qice{0:03d}'.format(ik)
  print('Updating ' + fldout)
  mc6util.modify_fld_nc(fl_restart6,fldout,qice_lr)

# Snow enthalpy by layers
for ik in range(1,cice6.nslyr+1):
  qsnon_lr = qsnon[ik-1,:,:,:]
  fldout = 'qsno{0:03d}'.format(ik)
  print('Updating ' + fldout)
  mc6util.modify_fld_nc(fl_restart6, fldout, qsnon_lr)

#
# Change restart date:
print('Changing global attributes: restart time to {0}/{1:02d}/{2:02d} {3} sec'.\
       format(YRnew, MMnew, MDnew, int(HRnew*3600.)))
mc6util.modify_glattr_nc(fl_restart6, 'myear',  int(YRnew))
mc6util.modify_glattr_nc(fl_restart6, 'mmonth', int(MMnew))
mc6util.modify_glattr_nc(fl_restart6, 'mday',   int(MDnew))
mc6util.modify_glattr_nc(fl_restart6, 'msec',   int(HRnew*3600.))

# Add info:
# Do not use for now - the code freezes up
inf1 = 'Restart created from GOFS3.2-93.0 CICE4: ' + cicerst4
inf2 = 'code: cice4_cice6_restart.py, NOAA EMC'
#mc6util.addnew_glattr_nc(fl_restart6, 'Info1', inf1)




