"""
  Layer thicknesses - SSH should be added to z*
"""
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
fldint     = "thknss"  # 
grid_shape = 'symmetr'   # MOM grid: symmetr/nonsymmetr


pthout  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_restart/'
pthtmp  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_tmp/'
pthgofs = '/work/Dmitry.Dukhovskoy/data/GOFS3.0/expt_19.0/'

grid_var = 'hgrid'
fldiout  = 'lrthknss'

print(f" INTERPOLATING on Z lev {nrun}-{expt} --> MOM {fldint}")

with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthrun  = dct[nrun][expt]["pthrun"]
pthgrid = dct[nrun][expt]["pthgrid"]
ftopo   = dct[nrun][expt]["ftopo"]
fgrid   = dct[nrun][expt]["fgrid"]

# Get lon/lat for correct variables:
if grid_var == 'hgrid':
  grid_hycom = 'ppnt'
elif grid_var == 'ugrid':
  grid_hycom = 'upnt'
elif grid_var == 'vgrid':
  grid_hycom = 'vpnt'
elif grid_var == 'qgrid':
  grid_hycom = 'qpnt'

#LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid, grdpnt=grid_hycom)
#HH = np.where(HH>=0, np.nan, HH)

huge   = 1.e25
rg     = 9806.
dnmb   = mtime.jday2dnmb(YR, jday)
dv     = mtime.datevec(dnmb)
ds     = mtime.datestr(dnmb)
MM     = dv[1]
DM     = dv[2]
rdate  = f"{YR}{MM:02d}{DM:02d}"

# Upload GOFS horiz interpolated fields:
def load_interp(pthtmp, fldiout, grid_shape, rdate):
  flintrp    = f"gofs2mom_nep_hrzi-{fldiout}_{grid_shape}_{rdate}.pkl"
  dflintrp   = os.path.join(pthtmp,flintrp)
  print(f"Loading {dflintrp}")
  with open(dflintrp, 'rb') as fid:
    Ai = pickle.load(fid)

  return Ai

#dHlr = load_interp(pthtmp, 'lrthknss', grid_shape, rdate)
SSH  = load_interp(pthtmp, 'srfhgt', grid_shape, rdate)

#j0 = 100
#i0 = 120
#dh = dHlr[:,j0,i0]

# Output interpolated fields on hycom layers:
flintrp  = f"gofs2mom_nep_vrtz-{fldiout}_{grid_shape}_{rdate}.pkl"
dflintrp = os.path.join(pthtmp,flintrp)

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
kdm         = 75
# Convert to -180/180:
hlon   = np.where(hlon > 180.0, hlon-360., hlon)

# Read restart from MOM restart template
pthrst = '/work/Dmitry.Dukhovskoy/NEP_input/restart/'
momrst = 'MOM.res.NEP_example.nc'
drst   = os.path.join(pthrst,momrst)
nc     = ncFile(drst, 'r')
dHlrM  = nc.variables['h'][:].data.squeeze()
SSHm   = nc.variables['sfc'][:].data.squeeze()
zI     = nc.variables['Interface'][:].data.squeeze()

Iall = np.where(HHM.flatten()<=0.)[0]
nall = Iall.shape[0]

Irsb = []
A3di = np.zeros((kdm,jdm,idm))*np.nan
kcc  = 0
tic  = timeit.default_timer()
ticR = timeit.default_timer()
print(f'Interpolating HYCOM GOFS --> MOM6 {fldint} {grid_shape} {grid_var}')

import mod_interp1D as mintrp
importlib.reload(mintrp)

for ikk in range(nall):
  I1 = Iall[ikk]
  jj, ii = np.unravel_index(I1,HHM.shape)  # MOM indices
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

# HYCOM GOFS indices:
#  II = np.squeeze(INDX[jj,ii,:])
#  JJ = np.squeeze(JNDX[jj,ii,:])
#  II, JJ = mblnr.sort_gridcell_indx(II,JJ)

  x0   = hlon[jj,ii]
  y0   = hlat[jj,ii]
# From GOFS interpolated onto MOM grid:
  ssh0 = SSH[jj,ii]

#  dh0  = dHlr[:,jj,ii]
#  zzH  = np.cumsum(dh0)
#  zzH  = np.insert(zzH,0,0) # interface depths
#  nhf  = len(zzH)
#  zmH  = 0.5*(zzH[:nhf-1] + zzH[1:])

# Check depth from HYCOM vs MOM:
#  hbH = -np.sum(dh0)   

#  if abs(1.-abs(hb0/hbH)) > 0.7:
#    raise Exception(f"GOFS/MOM Bottom depth error: j={jj} i={ii} {hb0} vs {hbH}")
#    print(f"GOFS/MOM Bottom depth error: j={jj} i={ii} {hbH} vs {hb0}")

# SSH correction to thickness:
# See HYCOM-tools: archv2mom6res.f
#  qq  = (ssh0 + abs(hbH))/abs(hbH)
#  dh0 = dh0*qq 

# MOM layers:
# z* levels can also be obatined from vgrid_75_2m.nc
# truncating the last layer to the local bottom depth
  hb0  = HHM[jj,ii]
  dhm  = dHlrM[:,jj,ii]
  sshm = SSHm[jj,ii]
  qqm  = (sshm + abs(hb0))/abs(hb0)
  dhm0 = dhm/qqm  # remove old ssh correction from MOM layers

# Add new ssh correction to MOM layers:
  qqm_new = (ssh0 + abs(hb0))/abs(hb0)
  dhm     = dhm0*qqm_new

  A3di[:,jj,ii] = dhm

print('END: Saving to '+dflintrp)
with open(dflintrp,'wb') as fid:
  pickle.dump(A3di,fid)


f_chck=False
if f_chck:
  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])
  ax1.plot(H1d, zmH, '.-')
  ax1.plot(M1d, zmM, '.-')







