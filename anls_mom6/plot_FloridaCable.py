"""
  Plot Flrodia Cable transports from the model and obs

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pdb
import importlib
import time
#import struct
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap

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
from mod_utils_fig import bottom_text

import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)

YRM   = 2021
YRR   = 2023
nrun  = 'MOM6'  # MOM6, RTOFS, GOFS3.1
dday  = 5        # freq of transport estimates, days - matters only for RTOFS production
sctnm = 'FlorCabl'
fld2d = 'Unrm'
f_zlm = False    # plot Zulema's transport for 2023 RTOFS production
sfx = ' '

mS=1
dS=1
if nrun == 'MOM6':
  expt = '003'
  YR   = YRM
  dday = 5
elif nrun == 'RTOFS':
  expt = 'product' # 003 or product
  sfx   = 'n-24'   # RTOFS output fields used for transp. calculation
  YR   = YRR
  if dday == 1:
    mS=9
    dS=15
elif nrun == 'GOFS3.1':
  expt = '93.0'
  YR   = YRM
  dday = 7

dnmb1 = mtime.datenum([YR,mS,dS])
dnmb2 = mtime.datenum([YR,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)

print(f"\n Plotting {nrun}-{expt} {sctnm} {fld2d} \n")

import mod_misc1 as mmisc
import mod_mom6 as mom6util
import mod_colormaps as mcmp
import mod_read_hycom as mhycom
importlib.reload(mom6util)
importlib.reload(mmisc)
importlib.reload(mcmp)

import yaml

with open('paths_expts.yaml') as ff:
  dct = yaml.safe_load(ff)

pthrun  = dct[nrun][expt]["pthrun"]
pthoutp = dct[nrun][expt]["pthoutp"]
pthgrid = dct[nrun][expt]["pthgrid"]
ftopo   = dct[nrun][expt]["ftopo"]
fgrid   = dct[nrun][expt]["fgrid"]


if nrun == 'MOM6':
  floutp  = f"mom6-{expt}_{fld2d}VFlx_{dv1[0]}" + \
            f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
  fgrd_mom  = os.path.join(pthgrid, fgrid)
  ftopo_mom = os.path.join(pthgrid, ftopo)
  HH  = mom6util.read_mom6depth(ftopo_mom)
elif nrun == 'RTOFS':
  if fld2d == 'salt':
    fld = 'salin'
  elif fld2d == 'potT':
    fld = 'temp'
  floutp  = f"rtofs-{expt}_{fld}xsct_{dv1[0]}" + \
            f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
  _, _, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
elif nrun == 'GOFS3.1':
  if fld2d == 'salt':
    fld = 'saln'
  elif fld2d == 'potT':
    fld = 'temp'
  floutp  = f"gofs31-930_{fld}xsct_{dv1[0]}" + \
            f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
  _, _, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)


ffout = pthoutp + floutp
print('Loading ' + ffout)

with open(ffout, 'rb') as fid:
  F2D, UFLX = pickle.load(fid)

# Time
TM   = F2D.TM
if sfx == 'n-24':
  TM = TM-1

# Depth-integrated: full grid
VFlx  = UFLX.trnsp*1e-6  # 1D depth-integrated flow, Sv
VFtot = np.nansum(VFlx, axis=1)

# Cable transports:
pthobs = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/OBS/FloridaCable/'
flobs  = f'FloridaCable_{YR}.txt'
dflobs = os.path.join(pthobs,flobs)
# Zulema's estimate from Sept 14
ftrnsp_zlm = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/rtofs.prod/data/' + \
             'transpFlorCur_prod_n-24.txt'

VFobs, TMobs = mom6vld.read_flcable(YR,dflobs)
VFobs        = np.where(VFobs>1.e20, np.nan, VFobs)
if f_zlm:
  VFzlm        = mom6vld.read_flcableZ(ftrnsp_zlm)
  dayz1        = mtime.datenum([2023,9,14])
  TMz          = np.arange(dayz1, dayz1+len(VFzlm))


btx = 'plot_FloridaCable.py'
plt.ion()

import mod_utils as mutil
import plot_sect as psct

dv1 = mtime.datevec(TM[0])
dv2 = mtime.datevec(TM[-1])
YR1 = dv1[0]
MM1 = dv1[1]
DD1 = dv1[2]
YR2 = dv2[0]
MM2 = dv2[1]
DD2 = dv2[2]
stl = f"0.08 {nrun}-{expt} {sctnm}, Vol. Flux, Sv   " + \
      f"{YR1}/{MM1:02d}/{DD1:02d}-{YR2}/{MM2:02d}/{DD2:02d}"

# Plot depth-integrated transport
Tday    = TM-TM[0]
if f_zlm:
  Tdayz   = TMz-TM[0]
Tdayobs = TMobs-TM[0]

# Average over the same time interval
it0 = 0
if f_zlm:
  it1   = np.where(Tdayobs == 0)[0][0]
  if Tdayz[-1] <= Tday[-1]:
    it2   = np.where(Tdayobs == Tdayz[-1])[0][0]
    it3   = np.where(Tday == Tdayz[-1])[0][0]
    it4   = len(Tdayz)-1
  else:
    it2   = np.where(Tdayobs == Tday[-1])[0][0]
    it3   = len(Tday)-1
    it4   = np.where(Tdayz == Tday[-1])[0][0]
else:
  if Tdayobs[0] <= 0:
    it1 = np.where(Tdayobs == 0)[0][0]
  else:
    it0 = np.where(Tday == Tdayobs[0])
    it1 = 0

  if Tday[-1] > Tdayobs[-1]:
    it3 = np.where(Tday == Tdayobs[-1])[0][0] 
    it2 = len(Tdayobs)
  else:
    it3 = len(Tday)
    it2 = np.where(Tdayobs == Tday[-1])[0][0]

mnVF  = np.mean(VFtot[it0:it3])
stdVF = np.std(VFtot)
if f_zlm:
  mnVFz = np.mean(VFzlm[it0:it4])
mnVFo = np.nanmean(VFobs[it1:it2])



# Depth-integrated transport
plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])
ln1, = ax1.plot(Tday, VFtot, color=[0.,0.4,0.8], label="RTOFS DD")
if f_zlm:
  ln2, = ax1.plot(Tdayz, VFzlm, color=[0.8,0.3,0], label="RTOFS ZG")
ln3, = ax1.plot(Tdayobs, VFobs, color=[0.3,0.3,0.3], label="Cable")

ax1.set_xlim([0, len(Tday)])
ax1.set_xticks(np.arange(0, len(Tday), 5))
ax1.grid(True)
ax1.set_xlabel(f'Days since {dv1[0]}/{dv1[1]}/{dv1[2]}')
ax1.set_ylabel(f'Transport, Sv')
ax1.set_title(stl)

sinfo = f'mean Flx {nrun}-{expt} (DD) = {mnVF:6.1f} Sv\n'
if f_zlm:
  sinfo = sinfo + 'mean Flx RTOFS-prod (ZG) = {0:6.1f} Sv\n'.format(mnVFz)
sinfo = sinfo + 'mean Flx Cable           = {0:6.1f} Sv'.format(mnVFo)


ax2 = plt.axes([0.6, 0.1, 0.3, 0.2])
if f_zlm:
  lgd = plt.legend(handles=[ln1,ln2,ln3], loc='upper left')
else:
  lgd = plt.legend(handles=[ln1,ln3], loc='upper left')

ax2.text(0., 0., sinfo)
ax2.axis('off')

bottom_text(btx,pos=[0.02, 0.03])

