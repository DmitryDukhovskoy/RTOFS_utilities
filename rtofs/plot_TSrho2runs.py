"""
 Analyze and plot T/S profiles at a given location
 Compute HYCOM rho compare with target densities
 from any other 2 profiles from different archv files
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

# Run 1:
rdate1  = '20230326'  # date of incrementally updated field
expt1   = 'paraD1'
sfx1    = 'n-24'
#pthbs1  = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/wcoss2.paraD1/'
pthbs1  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/paraD/'
pthbin1 = pthbs1 + 'rtofs.' + rdate1 + '/'
flnm1   = 'rtofs_glo.t00z.{0}.archv'.format(sfx1)
fin1a   = pthbin1 + flnm1 + '.a'
fin1b   = pthbin1 + flnm1 + '.b'
# Run 2:
rdate2  = '20230418'
expt2   = 'paraD2'
sfx2    = 'n-24'
#pthbs2  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/paraD1/'
pthbs2  = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/wcoss2.paraD1/'
pthbin2 = pthbs2 + 'rtofs.' + rdate2 + '/'
flnm2   = 'rtofs_glo.t00z.{0}.archv'.format(sfx2)
fin2a   = pthbin2 + flnm2 + '.a'
fin2b   = pthbin2 + flnm2 + '.b'

# location to plot:
iPrf = 601
jPrf = 1608
# Region for T/S :
# Sulu Sea:
IIr = [550, 580, 604, 610, 573, 542]
JJr = [1618, 1639, 1639, 1608, 1576, 1610]

import mod_time as mtime
jday1                = mtime.rdate2jday(rdate1)
yr1, mo1, mday1, hr1 = mtime.parse_rdate(rdate1)
dnmb1                = mtime.rdate2datenum(rdate1)
jday2                = mtime.rdate2jday(rdate2)
yr2, mo2, mday2, hr2 = mtime.parse_rdate(rdate2)
dnmb2                = mtime.rdate2datenum(rdate2)


plt.ion()


#
# Date of plotted fields:
if sfx1[0:2] == 'n-':
  hr = int(sfx1[2:])
  dnmb1P = dnmb1-float(hr)/24.
elif sfx1[0] == 'f':
  hr = int(sfx1[1:])
  dnmb1P = dnmb1+float(hr)/24.

if sfx2[0:2] == 'n-':
  hr = int(sfx2[2:])
  dnmb2P = dnmb2-float(hr)/24.
elif sfx2[0] == 'f':
  hr = int(sfx2[1:])
  dnmb2P = dnmb2+float(hr)/24.

dv1P = mtime.datevec(dnmb1P)
YR1p = dv1P[0]
MM1p = dv1P[1]
DD1p = dv1P[2]
HR1p = dv1P[3]
dv2P = mtime.datevec(dnmb2P)
YR2p = dv2P[0]
MM2p = dv2P[1]
DD2p = dv2P[2]
HR2p = dv2P[3]

import mod_read_hycom as mhycom
importlib.reload(mhycom)
from mod_utils_fig import bottom_text

print('Get dimensions from '+fin1a)
IDM, JDM, KDM = mhycom.hycom_dim(fin1a,fin1b)

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
ZZ, ZM = mhycom.get_zz_zm(fin1a,fin1b)
ZZ1prf = ZZ[:,jPrf,iPrf]
ZM1prf = ZM[:,jPrf,iPrf] 

ZZ, ZM = mhycom.get_zz_zm(fin2a,fin2b)
ZZ2prf = ZZ[:,jPrf,iPrf]
ZM2prf = ZM[:,jPrf,iPrf] 

# Read S & T for output 1 and output 2:
T1prf  = np.zeros((KDM))*np.nan
S1prf  = np.zeros((KDM))*np.nan
T2prf  = np.zeros((KDM))*np.nan
S2prf  = np.zeros((KDM))*np.nan
for lvl in range (1,KDM+1):
  F,n1,m1,l1 = mhycom.read_hycom(fin1a,fin1b,'temp',rLayer=lvl)
  F[np.where(F>huge)] = np.nan
  T1prf[lvl-1] = F[jPrf,iPrf]

for lvl in range (1,KDM+1):
  F,n1,m1,l1 = mhycom.read_hycom(fin1a,fin1b,'salin',rLayer=lvl)
  F[np.where(F>huge)] = np.nan
  S1prf[lvl-1] = F[jPrf,iPrf]

for lvl in range (1,KDM+1):
  F,n1,m1,l1 = mhycom.read_hycom(fin2a,fin2b,'temp',rLayer=lvl)
  F[np.where(F>huge)] = np.nan
  T2prf[lvl-1] = F[jPrf,iPrf]

for lvl in range (1,KDM+1):
  F,n1,m1,l1 = mhycom.read_hycom(fin2a,fin2b,'salin',rLayer=lvl)
  F[np.where(F>huge)] = np.nan
  S2prf[lvl-1] = F[jPrf,iPrf]
#
#
# Compute sigma using HYCOM eq. of state:
import mod_eqState as meqst
sig2_1  = meqst.sig2_17t_db(T1prf, S1prf)
sig2_2 = meqst.sig2_17t_db(T2prf, S2prf)

# Read target densities:
TDENS = mhycom.read_targ_dens(fin1b, 'salin', '=', 4)


def put_lrnumb(axx,ZMprf,xfld,zzS,clr0):
  for ikl in range(KDM):
    zm0 = ZMprf[ikl]
    sm0 = xfld[ikl]
    if zm0 > zzS or np.isnan(zm0):
      continue

    stt = '  Lr {0}'.format(ikl+1)
    axx.text(sm0, zm0-10., stt, color=clr0)


# ===================
#
#   Plot profiles
# 
# ====================
# Define T, S limits:
zz0  = -3000.
nInt = 5
th1, th2, dTh = mutil.prof_limits(T1prf, ZM1prf, nT=nInt, z0=zz0)
ta1, ta2, dTa = mutil.prof_limits(T2prf, ZM2prf, nT=nInt, z0=zz0)
tlim1 = min([th1,ta1])
tlim2 = max([th2,ta2])
dltT  = 0.1*(int((tlim2-tlim1)*10/nInt))

sh1, sh2, dSh = mutil.prof_limits(S1prf,ZM1prf,z0=zz0,ndec=1)
sa1, sa2, dSa = mutil.prof_limits(S2prf,ZM2prf,z0=zz0,ndec=1)
slim1 = min([sh1,sa1])
slim2 = max([sh2,sa2])
dltS  = 0.1*(int((slim2-slim1)*10/nInt))

# Plot T/S 
clra  = [0.,0.4,0.8]  # RTOFS run 1
clrh  = [0.8,0.2,0]   # RTOFS run 2
clra0 = [0.,0.9,0.4]  # 
clr2  = [0.8,0.,1]  # 

ctlT = 'Temp ' + rdate1 + sfx1 + ' '+ rdate2 + sfx2 
ctlS = 'Salin ' + rdate1 + sfx1 + ' '+ rdate2 + sfx2 

xpos1  = [0.07, 0.15, 0.34, 0.8]
xpos2  = [0.48, 0.15, 0.34, 0.8]
xpos23 = [0.09, 0.165, 0.18, 0.1]
xpos3  = [0.83, 0.6, 0.15, 0.3]
xpos4  = [0.5, 0.02, 0.45, 0.12] 

sinfo = 'Run1: ' + expt1 + ' ' + rdate1 + ' '
sinfo = sinfo + flnm1 + '\n'
sinfo = sinfo + 'Run2: ' + expt2 + ' ' + rdate2 + ' '
sinfo = sinfo + flnm2 + '\n'
sinfo = sinfo + 'iF={0} jF={1}\n'.format(iPrf+1,jPrf+1)
sinfo = sinfo + 'lon={0:5.2f}  lat={1:5.2f}'.format(lonPrf,latPrf)

Hb   = HH[jPrf,iPrf]
zlim = np.floor(Hb-10.)

#####
# Plot T/S Profiles

fig1 = plt.figure(1,figsize=(9,8), constrained_layout=False)
fig1.clf()
#
# T profile
#ax1 = plt.axes(xpos1)
ax1 = fig1.add_axes(xpos1)
ax1.plot(T1prf,ZM1prf, '.-', color=clra, linewidth=2, label='Run1')
ax1.plot(T2prf,ZM2prf, '.-', color=clrh, linewidth=2, label='Run2')
ax1.plot([tlim1,tlim2],[Hb,Hb],'--',color=[0,0,0])
# Put layer numbers:
zzS = -300.
put_lrnumb(ax1,ZM1prf,T1prf,zzS,clra)
put_lrnumb(ax1,ZM2prf,T2prf,zzS,clrh)

ax1.set_yticks(np.arange(-4000,0,200))
ax1.set_xticks(np.arange(tlim1,tlim2,dltT))
ax1.set_ylim(zlim,0)
ax1.set_xlim(tlim1,tlim2)
ax1.grid(True)
ax1.set_title(ctlT)

# S profile
#ax2 = plt.axes(xpos2)
ax2 = fig1.add_axes(xpos2)
#plt.plot(Sar, ZSar, '.-', color=clra, linewidth=2, label="Argo")
ln1, = ax2.plot(S1prf, ZM1prf, '.-', color=clra, linewidth=2, label='Run1')
ln2, = ax2.plot(S2prf, ZM2prf, '.-', color=clrh, linewidth=2, label='Run2')
ax2.plot([slim1,slim2],[Hb,Hb],'--',color=[0,0,0])
put_lrnumb(ax2,ZM1prf,S1prf,zzS,clra)
put_lrnumb(ax2,ZM2prf,S2prf,zzS,clrh)

ax2.set_yticks(np.arange(-4000,0,200))
ax2.set_xticks(np.arange(slim1,slim2,dltS))
ax2.set_ylim(zlim,0)
ax2.set_xlim(slim1,slim2)
ax2.grid(True)
ax2.set_title(ctlS)

#ax3 = plt.axes(xpos3)
ax3 = fig1.add_axes(xpos3)
lgd = plt.legend(handles=[ln1,ln2], loc='upper left')
ax3.axis('off')

ax24 = fig1.add_axes(xpos4)
ax24.text(0,0,sinfo)
ax24.axis('off')

btx = 'plot_TSrho2runs.py'
bottom_text(btx, pos=[0.02, 0.02])


#####
# Plot target densities and densitites from T/S profiles:
# Allowed density deviation from the target density:
# see blkdat.input: sigjump 
# 0.02   'sigjmp' = minimum density jump across interfaces  (kg/m**3)
sigjump = 0.02
ctl21 = 'sigma2: ' + rdate1 + sfx1 + ' '+ rdate2 + sfx2
ctl22 = 'sigma2-TDENS: ' + rdate1 + sfx1 + ' '+ rdate2 + sfx2

slgm1  = (np.floor(min(sig2_1)*10))/10.
slgm2  = (np.ceil(np.max(sig2_1)*10))/10.+0.25

fig2 = plt.figure(2, figsize=(9,8))
fig2.clf()

nsigma = 14   # sigma/z-layers
ax21 = plt.axes(xpos1)
ln21, = ax21.plot(sig2_1,ZM1prf,'.-', color=clra, linewidth=2, label='Run1')
ln22, = ax21.plot(TDENS, ZM1prf, '.-', color=clra0, linewidth=2, label="TDENS")
ln23, = ax21.plot(sig2_2,ZM2prf,'.-',color=clrh, linewidth=2, label='Run2')
ax21.plot([slgm1,slgm2],[Hb,Hb],'--',color=[0,0,0])
put_lrnumb(ax21,ZM1prf,sig2_1,zzS,clra)
put_lrnumb(ax21,ZM2prf,sig2_2,zzS,clrh)

ax21.set_yticks(np.arange(-4000,0,200))
ax21.set_ylim(zlim,0)
ax21.set_xlim(slgm1,slgm2)
ax21.grid(True)
ax21.set_title(ctl21)

ax22 = plt.axes(xpos2)
ax22.plot(sig2_1-TDENS, ZM1prf,'.-', color=clra, linewidth=2)
ax22.plot(sig2_2-TDENS, ZM2prf,'.-', color=clrh, linewidth=2)
ax22.set_xlim(-10.*sigjump, 10.*sigjump)
ax22.plot([-sigjump, -sigjump],[zlim, 0],'r--')
ax22.plot([sigjump, sigjump],[zlim, 0],'r--')
ax22.plot([-10.*sigjump, 10.*sigjump],[Hb,Hb],'--',color=[0,0,0])
put_lrnumb(ax22,ZM1prf,sig2_1-TDENS,zzS,clra)
put_lrnumb(ax22,ZM2prf,sig2_2-TDENS,zzS,clrh)

ax22.set_yticks(np.arange(-4000,0,200))
ax22.set_ylim(zlim,0)
ax22.grid(True)
ax22.set_title(ctl22)

ax23 = plt.axes(xpos23)
lgd2 = plt.legend(handles=[ln21,ln22,ln23], loc='lower left')
ax23.axis('off')

ax24 = plt.axes(xpos4)
ax24.text(0,0,sinfo)
ax24.axis('off')

bottom_text(btx, pos=[0.01, 0.02])

