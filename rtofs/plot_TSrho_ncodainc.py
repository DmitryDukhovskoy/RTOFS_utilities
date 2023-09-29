"""
 Analyze and plot T/S profiles at a given location
 Compute HYCOM rho compare with target densities
 From ncoda increment update files
 or any other 2 profiles from different archv files
 just change input directories accordingly

"""
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

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text

import mod_utils as mutil
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_read_ncoda as rncoda
importlib.reload(rncoda)
#import mod_extrargo as exargo
#importlib.reload(exargo)
from mod_utils import tsarray

np.set_printoptions(precision=3)

rdate0 = '20230418'  # date of incrementally updated field
expt   = 'paraD1'
sfx    = 'n-24'
rdate  = rdate0

# location to plot:
iPrf = 601
jPrf = 1608
# Region for T/S :
# Sulu Sea:
IIr = [550, 580, 604, 610, 573, 542]
JJr = [1618, 1639, 1639, 1608, 1576, 1610]

import mod_time as mtime
#jday = rncoda.rdate2julian(rdate)
#yr, mo, mday, hr = rncoda.parse_rdate(rdate)
jday             = mtime.rdate2jday(rdate)
yr, mo, mday, hr = mtime.parse_rdate(rdate)


plt.ion()

pthbs   = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/'
pthbin  = pthbs+'hycom/ncoda_archv_inc/'
pthbgr  = pthbs+'hycom/ncoda_archv_inc/'  # background field dir

#fhcm   = 'rtofs_glo.t00z.' + sfx + '.archv'
flupdt = 'archv_1_inc.{0}_{1:03d}_00'.format(yr,jday)
flbgr  = 'archv.{0}_{1:03d}_00'.format(yr,jday)

fina  = pthbin + flupdt + '.a'
finb  = pthbin + flupdt + '.b'
fbgra = pthbgr + flbgr + '.a'
fbgrb = pthbgr + flbgr + '.b'

#importlib.reload(mmisc)
import mod_time as mtime
YR     = int(rdate0[0:4])
MM     = int(rdate0[4:6])
DD     = int(rdate0[6:8])
#  yrday  = mmisc.date_yearday(YR,MM,DD)
yrday  = mtime.rdate2jday(rdate0)
dnmb0  = mtime.rdate2datenum(rdate0)

#
# Date of plotted fields:
if sfx == 'n-24':
  dnmbP = dnmb0-1
elif sfx[0] == 'f':
  hr = int(sfx[1:])
  dnmbP = dnmb0+float(hr)/24.

dvP = mtime.datevec(dnmbP)
YRp = dvP[0]
MMp = dvP[1]
DDp = dvP[2]
HRp = dvP[3]

import mod_read_hycom as mhycom
importlib.reload(mhycom)
from mod_utils_fig import bottom_text

print('Processing '+fina)
IDM, JDM, KDM = mhycom.hycom_dim(fina,finb)


pthgrid= '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
lonPrf = LON[jPrf,iPrf]
latPrf = LAT[jPrf,iPrf]


huge = 1.e20
rg   = 9806.

# Read layer pressures:
#ZZ, ZM = mhycom.zz_zm_fromDP(dH)
ZZ, ZM = mhycom.get_zz_zm(fina,finb)
ZZprf = ZZ[:,jPrf,iPrf]
ZMprf = ZM[:,jPrf,iPrf] 

# Read S or T:
Tprf  = np.zeros((KDM))*np.nan
Sprf  = np.zeros((KDM))*np.nan
for lvl in range (1,KDM+1):
  F,n1,m1,l1 = mhycom.read_hycom(fina,finb,'temp',rLayer=lvl)
  F[np.where(F>huge)] = np.nan
  Tprf[lvl-1] = F[jPrf,iPrf]

for lvl in range (1,KDM+1):
  F,n1,m1,l1 = mhycom.read_hycom(fina,finb,'salin',rLayer=lvl)
  F[np.where(F>huge)] = np.nan
  Sprf[lvl-1] = F[jPrf,iPrf]

#
# Background field:
ZZ, ZM = mhycom.get_zz_zm(fbgra,fbgrb)
ZZBprf = ZZ[:,jPrf,iPrf]
ZMBprf = ZM[:,jPrf,iPrf]

# Read S or T:
TBprf  = np.zeros((KDM))*np.nan
SBprf  = np.zeros((KDM))*np.nan
for lvl in range (1,KDM+1):
  F,n1,m1,l1 = mhycom.read_hycom(fbgra,fbgrb,'temp',rLayer=lvl)
  F[np.where(F>huge)] = np.nan
  TBprf[lvl-1] = F[jPrf,iPrf]

for lvl in range (1,KDM+1):
  F,n1,m1,l1 = mhycom.read_hycom(fbgra,fbgrb,'salin',rLayer=lvl)
  F[np.where(F>huge)] = np.nan
  SBprf[lvl-1] = F[jPrf,iPrf]

#
# Compute sigma using HYCOM eq. of state:
import mod_eqState as meqst
sig2  = meqst.sig2_17t_db(Tprf, Sprf)
sigb2 = meqst.sig2_17t_db(TBprf, SBprf)

# Read target densities:
TDENS = mhycom.read_targ_dens(finb, 'salin', '=', 4)

# ===================
#
#   Plot profiles
# 
# ====================
# Define T, S limits:
zz0  = -3000.
nInt = 5
th1, th2, dTh = mutil.prof_limits(Tprf, ZMprf, nT=nInt, z0=zz0)
ta1, ta2, dTa = mutil.prof_limits(TBprf, ZMBprf, nT=nInt, z0=zz0)
tlim1 = min([th1,ta1])
tlim2 = max([th2,ta2])
dltT  = 0.1*(int((tlim2-tlim1)*10/nInt))

sh1, sh2, dSh = mutil.prof_limits(Sprf,ZMprf,z0=zz0,ndec=1)
sa1, sa2, dSa = mutil.prof_limits(SBprf,ZMBprf,z0=zz0,ndec=1)
slim1 = min([sh1,sa1])
slim2 = max([sh2,sa2])
dltS  = 0.1*(int((slim2-slim1)*10/nInt))

# Plot T/S 
clra  = [0.,0.4,0.8]  # RTOFS parallel run
clrh  = [0.8,0.2,0]   # RTOFS production
clra0 = [0.,0.9,0.4]  # GDEM climatology
clr2  = [0.8,0.,1]  # 

ctl1 = 'Temp '+rdate0
ctl2 = 'Salin '+rdate0

xpos1  = [0.07, 0.15, 0.34, 0.8]
xpos2  = [0.48, 0.15, 0.34, 0.8]
xpos23 = [0.09, 0.165, 0.18, 0.1]
xpos3  = [0.83, 0.6, 0.15, 0.3]
xpos4  = [0.6, 0.02, 0.4, 0.12] 

sinfo = expt + ' ' + rdate0 + '\n'
sinfo = sinfo + flupdt + '\n'
sinfo = sinfo + 'iF={0} jF={1}\n'.format(iPrf+1,jPrf+1)
sinfo = sinfo + 'lon={0:5.2f}  lat={1:5.2f}'.format(lonPrf,latPrf)

#####
# Plot T/S Profiles

fig1 = plt.figure(1,figsize=(9,8), constrained_layout=False)
fig1.clf()
#
# T profile
#ax1 = plt.axes(xpos1)
ax1 = fig1.add_axes(xpos1)
ax1.plot(Tprf,ZMprf, '.-', color=clra, linewidth=2, label='incrupd')
ax1.plot(Tprf,ZMprf, '.-', color=clrh, linewidth=2, label='bkground')

ax1.set_yticks(np.arange(-4000,0,200))
ax1.set_xticks(np.arange(tlim1,tlim2,dltT))
ax1.set_ylim(zz0,0)
ax1.set_xlim(tlim1,tlim2)
ax1.grid(True)
ax1.set_title(ctl1)

# S profile
#ax2 = plt.axes(xpos2)
ax2 = fig1.add_axes(xpos2)
#plt.plot(Sar, ZSar, '.-', color=clra, linewidth=2, label="Argo")
ln1, = ax2.plot(Sprf, ZMprf, '.-', color=clra, linewidth=2, label='incrupd')
ln2, = ax2.plot(SBprf, ZMBprf, '.-', color=clrh, linewidth=2, label='bkground')

ax2.set_yticks(np.arange(-4000,0,200))
ax2.set_xticks(np.arange(slim1,slim2,dltS))
ax2.set_ylim(zz0,0)
ax2.set_xlim(slim1,slim2)
ax2.grid(True)
ax2.set_title(ctl2)

#ax3 = plt.axes(xpos3)
ax3 = fig1.add_axes(xpos3)
lgd = plt.legend(handles=[ln1,ln2], loc='upper left')
ax3.axis('off')


ax24 = fig1.add_axes(xpos4)
ax24.text(0,0,sinfo)
ax24.axis('off')

btx = 'plot_TSrho_ncodainc.py'
bottom_text(btx, pos=[0.07, 0.03])


#####
# Plot target densities and densitites from T/S profiles:
# Allowed density deviation from the target density:
# see blkdat.input: sigjump 
# 0.02   'sigjmp' = minimum density jump across interfaces  (kg/m**3)
sigjump = 0.02

Hb = np.floor(HH[jPrf,iPrf])
zlim   = Hb-10.
slgm1  = (np.floor(min(sig2)*10))/10.
slgm2  = (np.ceil(np.max(sig2)*10))/10.


ctl21 = 'sigma2 & targ dens ' + rdate0
fig2 = plt.figure(2, figsize=(9,8))
fig2.clf()

nsigma = 14   # sigma/z-layers
ax21 = plt.axes(xpos1)
ln21, = ax21.plot(sig2,ZMprf,'.-', color=clra, linewidth=3, label='incrupd')
ln22, = ax21.plot(TDENS, ZMprf, '.-', color=clra0, linewidth=2, label="TDENS")
ln23, = ax21.plot(sigb2,ZMBprf,'.-',color=clrh, linewidth=2, label="bkground")
ax21.set_yticks(np.arange(-4000,0,500))
ax21.set_ylim(zlim,0)
ax21.set_xlim(slgm1,slgm2)
ax21.set_ylim(-2900,0)
ax21.grid(True)
ax21.set_title(ctl21)

ctl22='sigma2-TDENS'
ax22 = plt.axes(xpos2)
ax22.plot(sig2-TDENS, ZMprf,'.-', color=clra, linewidth=2)
ax22.plot(sigb2-TDENS, ZMBprf,'.-', color=clrh, linewidth=2)
ax22.set_xlim(-10.*sigjump, 10.*sigjump)
ax22.plot([-sigjump, -sigjump],[zlim, 0],'r--')
ax22.plot([sigjump, sigjump],[zlim, 0],'r--')
ax22.set_ylim(-2900,0)
ax22.set_yticks(np.arange(-4000,0,500))
ax22.set_ylim(zlim,0)
ax22.grid(True)

for ikl in range(KDM):
  zm0 = ZMprf[ikl]
  sm0 = sig2[ikl]-TDENS[ikl]
  if zm0 > -300 or np.isnan(zm0):
    continue

  stt = '  Lr {0}'.format(ikl+1)
  ax22.text(sm0, zm0-10., stt, color=clra)


for ikl in range(KDM):
  zm0 = ZMBprf[ikl]
  sm0 = sigb2[ikl]-TDENS[ikl]
  if zm0 > -300 or np.isnan(zm0):
    continue

  stt = '  Lr {0}'.format(ikl+1)
  ax22.text(sm0, zm0-10, stt, color=clrh)

ax22.set_title(ctl22)

ax23 = plt.axes(xpos23)
lgd2 = plt.legend(handles=[ln21,ln22,ln23], loc='lower left')
ax23.axis('off')

ax24 = plt.axes(xpos4)
ax24.text(0,0,sinfo)
ax24.axis('off')

bottom_text(btx, pos=[0.07, 0.03])

