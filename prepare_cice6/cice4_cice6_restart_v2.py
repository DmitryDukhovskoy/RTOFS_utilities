# Convert binary cice4 to binary netcdf cice6 restart
#
# Dmitry Dukhovskoy NOAA NWS EMC
# March 2023
#
# December 2025
# Updated: 
#   input key arguments
#   netCDF processing based on xarray
#   velocities at E points for C grid
#
import os
import numpy as np
import sys
import importlib
import datetime
import xarray
import time
from yaml import safe_load
import argparse

PPTHN = None  # directory with custom modules, full path
if 'PPTHN' not in locals() or PPTHN is None:
  # Try default location of directory with python modules MyPython: ../
  cwd = os.getcwd()    
  parts = cwd.split(os.sep)
  if 'python' in parts:
    idx = parts.index('python')
    PPTHN = os.sep + os.path.join(*parts[:idx + 1])
  else:
    raise RuntimeError("Directory 'python' not found in current working directory path.")
  
sys.path.extend([
    os.path.join(PPTHN, 'MyPython', 'hycom_utils'),
    os.path.join(PPTHN, 'MyPython', 'draw_map'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])


import mod_time as mtime
import mod_cice6_utils as mc6util
importlib.reload(mc6util)

# Default values:
fyaml = 'cice4_cice6.yaml'
rdateT = 2025050800        # restart date and time in template cice6 restart

parser = argparse.ArgumentParser()
parser.add_argument("--fyaml", help=f"yaml file with paths, filenames, params, default={fyaml}", type=str)
parser.add_argument("--rdate6", help="Restart date in CICE6: YYYYMMDDhh", required=True, type=int)
parser.add_argument("--rdateT", help=f"Restart date in template: YYYYMMDDhh, default={rdateT}", type=int)
args = parser.parse_args()

fyaml  = args.fyaml if args.fyaml else fyaml
rdate6 = args.rdate6 if args.rdate6 else None
rdateT = args.rdateT if args.rdateT else rdateT

syst_info = os.uname() 
machine = syst_info.nodename
  
if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
elif 'an' in machine:
  print("Running on PPAN node:", machine)
  node_nm = "ppan"
elif 'ufe' in machine:
  print("Running on Ursa node:", machine)
  node_nm = "ursa"
else:
  print("Unknown machine:", machine)

# Restart CICE6 template:
dnmbT = mtime.dateint2datenum(rdateT)
YRtmp, MMtmp, MDtmp, HRtmp, _ = mtime.datevec(dnmbT, round_hrs=True)

# Output restart for CICE6:
dnmb6 = mtime.dateint2datenum(rdate6)
YRc6, MMc6, MDc6, HRc6, _ = mtime.datevec(dnmb6, round_hrs=True)

# Get input/output paths for CICE restart files:
with open(fyaml) as ff:
  PATHS = safe_load(ff)

cicerst4 = PATHS["rest_names"]["cice4"]["flnm"]
cicerstT = PATHS["rest_names"]["tmplt"]["flnm"].format(yr=YRtmp, mm=MMtmp, dd=MDtmp, hr=HRtmp)
cicerst6 = PATHS["rest_names"]["cice6"]["flnm"].format(yr=YRc6, mm=MMc6, dd=MDc6, hr=HRc6)
pthrst4  = PATHS["cice_paths"][node_nm]["cice4"]["pth"]
pthrstT  = PATHS["cice_paths"][node_nm]["tmplt"]["pth"]
pthrst6  = PATHS["cice_paths"][node_nm]["cice6"]["pth"]

# Restart files with dirs:
fl_restart4 = os.path.join(pthrst4, cicerst4)
fl_restartT = os.path.join(pthrstT, cicerstT)
fl_restart6 = os.path.join(pthrst6, cicerst6)

print(f'Input CICE4 restart:     {fl_restart4}')
print(f'Template CICE6 restart:  {fl_restartT}')
print(f'Output CICE6 restart:    {fl_restart6}')

# Read CICE params:
ice_grid4 = PATHS["cice_params"]["cice4"]["grid"]
ice_grid6 = PATHS["cice_params"]["cice6"]["grid"]

#import mod_datm_utils as mdatm
#importlib.reload(mdatm)

# Create new restart from template for writing CICE4 fields
#mdatm.cice6_newfile(fl_restartT, fl_restart6)

# Grid CICE4 - unformatted binary file
pthgrd4 = PATHS["grid_topo"][node_nm]["cice4"]["pthgrid"]
grdfl4  = PATHS["grid_topo"][node_nm]["cice4"]["filegrid"]
fgrdin4 = os.path.join(pthgrd4, grdfl4)

# CICE6 grid
pthgrd  = PATHS["grid_topo"][node_nm]["cice6"]["pthgrid"]
grdfl   = PATHS["grid_topo"][node_nm]["cice6"]["filegrid"]
fgrdin  = os.path.join(pthgrd, grdfl)

# depth:
pthdpth = PATHS["grid_topo"][node_nm]["cice6"]["pthtopo"]
dpthfl  = PATHS["grid_topo"][node_nm]["cice6"]["filedepth"]

# Check if this is .a, .b or .nc depth file:
fldptha  = fldpthb = None
topo_nc = False
topo_ab = False

if dpthfl.endswith('.nc'):
  fldpthnc = pthdpth
  fdpthin = os.path.join(pthdpth, dpthfl)
  topo_nc = True
elif dpthfl.endswith('.a'):
  fldptha = dpthfl
  fldpthb = fldptha.replace('.a', '.b')
  ftopo   = fldptha.removesuffix('.a')
  #fdpthin_a = os.path.join(pthdpth, fldptha)
  #fdpthin_b = os.path.join(pthdpth, fldpthb)
  topo_ab = True
else:
  raise ValueError(f"topo file {dpthfl} not recognized, expected *.a or *.nc")

def read_ncfield(dirflnc, varnc):
  with xarray.open_dataset(dirflnc) as dset:
    AA = dset[varnc].data.squeeze()

  return AA  

def read_rest_cice4(fid, nx, ny):
  """
    CICE4 restart fields
    Read 2D record from open file
    unformatted sequential binary file 
    bing endian
  """
  recS = np.fromfile(fid, dtype='>i4', count=1)[0]
  A     = np.fromfile(fid, dtype='>f8', count=nx*ny)
  A     = np.reshape(A,(ny,nx), order='C')
  recE = np.fromfile(fid, dtype='>i4', count=1)[0]
  if recS != recE:
    raise ValueError(f"Record length mismatch: {recS} != {recE}")

  return A

def print_minmax(sfld,A):
  print(f'   {sfld} min/max:  {np.nanmin(A)}/{np.nanmax(A)}')
  return


print(f'Creating CICE6 restart for {YRc6}/{MMc6:02d}/{MDc6:02d} {HRc6:02d}hr UTC')
print(f'CICE4 restart:     fl_restart4')
print(f'CICE6 template:    fl_restartT')
print(f'New CICE6 restart: fl_restart6')
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

#  Input parameters - check with ice_in
#  Edit mod_cice6_utils param_cice4 if needed
# Note in GOFS3.1 CICE has 1 row less than HYCOM
# ny = 3297 (and it is 3298 in HYCOM)
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

spval = 1.e30

if os.path.exists(fl_restart4):
  fid = open(fl_restart4, 'rb')
else:
  raise FileNotFoundError(f"Does not exist: {fl_restart4}")

print(f'Reading restart: {fl_restart4}')

fid.seek(0)
# Read Fortran binary
recS    = np.fromfile(fid, dtype='>i4', count=1)[0]   # record marker Start
istep   = np.fromfile(fid, dtype='>i4', count=1)[0]
runtime = np.fromfile(fid, dtype='>f8', count=1)[0]  # total elapsed time, sec
frtime  = np.fromfile(fid, dtype='>f8', count=1)[0]  # forcing time, sec
recE    = np.fromfile(fid, dtype='>i4', count=1)[0]  # record marker End

if recS != recE:
  raise ValueError(f"Record length mismatch: {recS} != {recE}")

print('Restart: step={0}, total time(yrs)={1}, forcing last update(hrs)={2}'\
       .format(istep,runtime/(3600*24*365.25),frtime/3600.))

# Read state variables:
# Tsfc is the only tracer read in this file
nx     = cice4.nx
ny     = cice4.ny
ncat   = cice4.ncat
ntilyr = cice4.ntilyr  # total # of icelrs * cat 
ntslyr = cice4.ntslyr

if topo_nc:
  Lmsk  = mc6util.read_ncfile(fdpthin,'wet')
else:
  _, Lmsk = mc6util.read_topo_ab(pthdpth, ftopo, nx, ny, lmask=True)


aicen = np.zeros((ncat,ny,nx), dtype='float64')
vicen = np.zeros((ncat,ny,nx), dtype='float64')
vsnon = np.zeros((ncat,ny,nx), dtype='float64')
trcrn = np.zeros((ncat,ny,nx), dtype='float64')

for n in range(ncat):
  print(' Category {0}'.format(n+1))
# Read ice area for category n
  A = read_rest_cice4(fid, nx, ny)
  aicen[n,:,:] = A
  print_minmax('ice area',A)

# Read ice volume/per m2/ for category n
  A = read_rest_cice4(fid, nx, ny)
  vicen[n,:,:] = A
  print_minmax('ice vol',A)

# Read snow volume/m2 for category n
  A = read_rest_cice4(fid, nx, ny)
  vsnon[n,:,:] = A
  print_minmax('snow vol',A)

# Read tracer 1 (surf T = Tsfcn by categories) - only 1 tracer in CICE4
  A = read_rest_cice4(fid, nx, ny)
  trcrn[n,:,:] = A
  print_minmax('surf T',A)

# Ice energy (J/m2):
eicen = np.zeros((ntilyr,ny,nx), dtype='float64')
print('\n Ice energy:')

for k in range(ntilyr):
  A = read_rest_cice4(fid,nx,ny)
  eicen[k,:,:] = A
  print_minmax(f'{k+1} eicen',A)

# Snow energy:
esnon = np.zeros((ntslyr,ny,nx), dtype='float64')
print('\n Snow energy:')

for k in range(ntslyr):
  A = read_rest_cice4(fid,nx,ny)
  esnon[k,:,:] = A
  print_minmax('{0} esnon'.format(k+1),A)

# Velocities:
print('\n Velocity components:')
A = read_rest_cice4(fid,nx,ny)
uvel = A.copy()
print_minmax('U vel',A)

A = read_rest_cice4(fid,nx,ny)
vvel = A.copy()
print_minmax('V vel',A)

# For C-grid need vel components on N/E points of the grid cell
# For now, simple linear interpolation
# TODO: apply higher-order interpolation polynomial for 
# deriving uvelE, vvelN
uvelE = None
vvelN = None
if ice_grid6 == 'C':
  uvelE, vvelN = mc6util.interp_uvelE_vvelN(uvel, vvel, aicen) 

# Radiation fields
# 4 radiative categories
# for calculating albedo for visible and IR wavelengths
# and penetrating sh/wave

# Scale factor to change MKS units
# for shortwave components
# default = 1
print('\n Radiation fields, W/m2: ')
A = read_rest_cice4(fid,nx,ny)
scale_factor = A.copy()
print_minmax('Scale Factor',A)

A = read_rest_cice4(fid,nx,ny)
swvdr = A.copy()
print_minmax('Sh/wave down vis. direct',A)

A = read_rest_cice4(fid,nx,ny)
swvdf = A.copy()
print_minmax('Sh/wave down vis. diff',A)

A = read_rest_cice4(fid,nx,ny)
swidr = A.copy()
print_minmax('Sh/wave down near IR dir',A)

A = read_rest_cice4(fid,nx,ny)
swidf = A.copy()
print_minmax('Sh/wave down near IR diff',A)

# Ocean stress, N/m2
print('\n Ocean stress components, N/m2:')

A = read_rest_cice4(fid,nx,ny)
strocnxT = A.copy()
print_minmax('ocean stress x-comp',A)

A = read_rest_cice4(fid,nx,ny)
strocnyT = A.copy()
print_minmax('ocean stress y-comp',A)

# Internal stress, stress tensor kg/s2
# (1) northeast, (2) northwest, (3) southwest, (4) southeast
print('\n Internal stress, kg/s2: ')

A = read_rest_cice4(fid,nx,ny)
stressp_1 = A.copy()
print_minmax('stressp_1',A)

A = read_rest_cice4(fid,nx,ny)
stressp_3 = A.copy()
print_minmax('stressp_3',A)

A = read_rest_cice4(fid,nx,ny)
stressp_2 = A.copy()
print_minmax('stressp_2',A)

A = read_rest_cice4(fid,nx,ny)
stressp_4 = A.copy()
print_minmax('stressp_4',A)

A = read_rest_cice4(fid,nx,ny)
stressm_1 = A.copy()
print_minmax('stressm_1',A)

A = read_rest_cice4(fid,nx,ny)
stressm_3 = A.copy()
print_minmax('stressm_3',A)

A = read_rest_cice4(fid,nx,ny)
stressm_2 = A.copy()
print_minmax('stressm_2',A)

A = read_rest_cice4(fid,nx,ny)
stressm_4 = A.copy()
print_minmax('stressm_4',A)

A = read_rest_cice4(fid,nx,ny)
stress12_1 = A.copy()
print_minmax('stress12_1',A)

A = read_rest_cice4(fid,nx,ny)
stress12_3 = A.copy()
print_minmax('stress12_3',A)

A = read_rest_cice4(fid,nx,ny)
stress12_2 = A.copy()
print_minmax('stress12_2',A)

A = read_rest_cice4(fid,nx,ny)
stress12_4 = A.copy()
print_minmax('stress12_4',A)

# Ice mask for dynamics
print('\n Ice Mask for Dynamics: ')

A = read_rest_cice4(fid,nx,ny)
iceumask = A.copy()
print_minmax('iceumask',A)

# For trully coupled HYCOM-CICE these fields
# are not needed
# if defined ocean mixed layer in CICE
# This is for ocean mixed layer defined (for GOFS )
print('\n Ocean mixed layer: \n');

A = read_rest_cice4(fid,nx,ny)
sst = A.copy()
print_minmax('sst',A)

A = read_rest_cice4(fid,nx,ny)
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

cice6 = mc6util.cice6()

btx = 'cice4_cice6_restart_v2.py'
  
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
ulati6 = read_ncfield(fgrdin, 'ulat')
uloni6 = read_ncfield(fgrdin, 'ulon')

# In CICE4 regional grid:
# last column (in the top ~100 rows) in lat is repeated (end-1) column
# in CICE6, these are different columns
mc6util.check_cice_grids(ulati4,uloni4,ulati6,uloni6)

# Compute coszen for CICE4 restart data
# adjusting for new restart date
# UTC --> time_zone=0
# ==> TO DO: Add land mask = 0 to match CICE6 restart 
DVnew      = [int(YRc6),int(MMc6),int(MDc6),int(HRc6)]
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
apnd = aicen*0.0

# hpnd - melt pond depth, m
# cannot be recovered from CICE4
# hpnd =0
hpnd = aicen*0.0

# ipnd - melt pond refrozen lid thickness
# ipnd cannot be rcovered from CICE4
ipnd = aicen*0.0

# dhs - change in snow thickness
# icepack_therm_vertical.F90: 
# dhs = hsn - works - hsn_new, where hsn_new - thickness of new snow (snowfall)
# Cannot be recovered for CICE4
dhs = aicen*0.0

# ffrac - fraction of fsurfn (atm-ice surface heat flux  (W/m2)) 
#         over pond used to melt ipond (pond ice)
# Cannot be recovered from CICE4
ffrac = aicen*0.0

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

# =======================================================
# 
#  scale_factor - scaling factor for shortwave radiation components
# use from CICE6 template? CICE4 scale_factor 
# is different
# scale_factor: netsw scaling factor (new netsw / old netsw)
# see: icepack_shortwave.F90
#
# Collect updated fields:
updated_vars = {
  'uvel':         uvel,
  'vvel':         vvel,
  'uvelE':        uvelE,
  'vvelN':        vvelN,
  'scale_factor': scale_factor,
  'swvdr':        swvdr,
  'swvdf':        swvdf,
  'swidr':        swidr,
  'swidf':        swidf,
  'strocnxT':     strocnxT,
  'strocnyT':     strocnyT,
  'stressp_1':    stressp_1, 
  'stressp_2':    stressp_2, 
  'stressp_3':    stressp_3, 
  'stressp_4':    stressp_4, 
  'stressm_1':    stressm_1, 
  'stressm_2':    stressm_2, 
  'stressm_3':    stressm_3, 
  'stressm_4':    stressm_4, 
  'stress12_1':   stress12_1, 
  'stress12_2':   stress12_2, 
  'stress12_3':   stress12_3, 
  'stress12_4':   stress12_4, 
  'iceumask':     iceumask,
  'fsnow':        fsnow,
  'aicen':        aicen,
  'vicen':        vicen,
  'vsnon':        vsnon,
  'Tsfcn':        trcrn,
  'coszen':       coszen_new,
  'iage':         iage,
  'alvl':         alvl,
  'vlvl':         vlvl,
  'apnd':         apnd,
  'hpnd':         hpnd,
  'ipnd':         ipnd,
  'dhs':          dhs,
  'ffrac':        ffrac
}

print(' \n\n -------------\n Creating CICE6 restart')
dst = xarray.open_dataset(fl_restartT)

for varname, new_data in updated_vars.items():
  if varname in dst:
    dst[varname] = xarray.DataArray(
      new_data,
      dims=dst[varname].dims,
      coords=dst[varname].coords
    )
  else:
    print(f"{varname} is not in {fl_restartT}")

# Add a new variable to restart file:
def add_newvar(dst, varname, A3d):
  new_fld = xarray.DataArray(A3d, 
                        dims=dst[varname].dims, 
                        coords=dst[varname].coords)
  dst[varname] = new_fld

  return dst

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
  varname = f'sice{ik:03d}'
  print(f'Updating {varname}')
  dst = add_newvar(dst, varname, sice_lr)

# Ice enthalpy by layers:
for ik in range(1,cice6.nilyr+1): 
  qice_lr = qicen[ik-1,:,:,:]
  varname = f'qice{ik:03d}'
  print(f'Updating {varname}')
  dst = add_newvar(dst, varname, qice_lr)
  #mc6util.modify_fld_nc(fl_restart6,fldout,qice_lr)

# Snow enthalpy by layers
for ik in range(1,cice6.nslyr+1):
  qsnon_lr = qsnon[ik-1,:,:,:]
  varname = f'qsno{ik:03d}'
  print(f'Updating {varname}')
  dst = add_newvar(dst, varname, qsnon_lr)
  #mc6util.modify_fld_nc(fl_restart6, fldout, qsnon_lr)

#
# Change restart date:
print(f'Changing global attributes: restart time to {YRc6}/{MMc6:02d}/{MDc6:02d} {HRc6*3600} sec')
dst.attrs['myear']  = np.int32(YRc6)
dst.attrs['mmonth'] = np.int32(MMc6)
dst.attrs['mday']   = np.int32(MDc6)
dst.attrs['msec']   = np.int32(HRc6*3600)
dst.attrs['info1']  = f"Restart created from CICE4: {cicerst4}"
dst.attrs['info2']  = f"code: {btx}"

print(f"Saving cice restart ---> {fl_restart6}")
dst.to_netcdf(fl_restart6, encoding={var: {'_FillValue': None} for var in dst.data_vars}, \
              format='NETCDF3_64BIT')

dst.close()
#ds6.close()

print(f'Created CICE6 restart: {fl_restart6}\n')

if not os.path.isfile(fl_restart6):
  raise Exception (f'ERR: CICE6 restart was NOT CREATED: {fl_restart6}')


