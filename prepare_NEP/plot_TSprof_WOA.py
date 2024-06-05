"""
  Plot MOM-NEP output T/S profiles
  with WOA23 profiles
  for specified regions
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
#import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mom6util
import mod_mom6_valid as mvalid 
#import mod_valid_utils as mvutil
#importlib.reload(mcmp)
importlib.reload(mvalid)

nrun   = "MOM6_NEP"
exptGF = "NEP_BGCphys_GOFS"  # GPFS restart
exptGL = "NEP_BGCphys"       # GLORYS-based restart
dnmb0  = mtime.datenum([1993,5,1])
#regn   = 'BeringS'
#regn   = 'GAlaska'
regn   = 'SouthCA'

with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

# Read MOM6 NEP domain nx=342, ny=816
pthgrid_mom = dct["MOM6_NEP"][exptGF]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"][exptGF]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"][exptGF]["fgrid"]
pthoutp     = dct["MOM6_NEP"][exptGF]["pthoutp"]
pthoutpGF   = pthoutp
pthoutpGL   = dct["MOM6_NEP"][exptGL]["pthoutp"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom)
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom)
LONM, LATM  = mom6util.read_mom6grid(dfgrid_mom)
HHM         = mom6util.read_mom6depth(dftopo_mom)
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]

# Regions:
if regn == 'BeringS':
  IJ = [[155, 710],
        [191, 710],
        [190, 625],
        [155, 625]]
  hmin = -40.
elif regn == 'GAlaska':
  IJ = [[175, 535],
        [235, 535],
        [235, 443],
        [175, 443]]
  hmin = -500.
elif regn == 'SouthCA':
  IJ = [[213, 216],
        [262, 216],
        [262, 41],
        [213, 41]]
  hmin = -500.

IJ = np.array(IJ)

def get_prof(dflmom, varplt, JP, IP):
  nc   = ncFile(dflmom,'r')
  A3d  = nc.variables[varplt][:].data.squeeze()
  A3d  = np.where(A3d > 1.e10, np.nan, A3d)
  dH   = nc.variables['h'][:].data.squeeze()
  #ZZ, ZM = mom6util.get_zz_zm(dflmom, sshvar='ssh')

  # Remove ssh from lr thicknesses:
  ssh = nc.variables['ssh'][:].squeeze().data
  dH  = mom6util.remove_ssh_lrthk(ssh, HHM, dH)
  _, ZM = mom6util.zz_zm_fromDP(dH, ssh, f_ssh = False)

  ZM0  = ZM[:,jdeep,ideep].squeeze()
  PRF  = A3d[:,JP,IP]
  ZPRF = ZM[:,JP,IP] 

  return ZM0, PRF, ZPRF

import mod_misc1 as mmisc
# Find regional mean and demean:
IIH, JJH  = np.meshgrid(np.arange(idm), np.arange(jdm))
_, IP, JP = mmisc.inpolygon_v2(IIH, JJH, IJ[:,0],IJ[:,1])

hmm   = HHM[JP,IP]
hmm   = np.where(hmm > hmin, np.nan, hmm)
isub  = np.where(~np.isnan(hmm))[0] 
IP    = IP[isub]
JP    = JP[isub]

jdeep, ideep = np.where(HHM == np.min(HHM))
hbtm = np.floor(np.min(HHM[JP,IP]))
hbtm = max([hbtm, -500])

DV      = mtime.datevec(dnmb0)
YR      = DV[0]
MM      = DV[1]
DD      = DV[2]
_, jday = mtime.dnmb2jday(dnmb0)
jday    = int(jday)

# Derive profiles:
flmom   = f'ocean_{YR}_{jday:03d}_12.nc'
dflmom  = os.path.join(pthoutpGF, f"{YR}",f"{MM:02d}", flmom)
ZMGF0, TGF, _ = get_prof(dflmom, "potT", JP, IP)
ZMGF0, SGF, _ = get_prof(dflmom, "salt", JP, IP)
dflmom  = os.path.join(pthoutpGL, f"{YR}",f"{MM:02d}", flmom)
ZMGL0, TGL, _ = get_prof(dflmom, "potT", JP, IP)
ZMGL0, SGL, _ = get_prof(dflmom, "salt", JP, IP)

# Profiles WOA:
# use seasonal profiles decadal means
if MM >=1 and MM <=3:
  seas = 13
elif MM >=4 and MM <=6:
  seas = 14
elif MM >=7 and MM <= 9:
  seas = 15
else:
  seas = 16

import mod_regmom as regmom
Xp = LONM[JP,IP]
Yp = LATM[JP,IP] 
TWOA, SWOA, ZWOA = regmom.derive_TSprof_WOA23(seas, YR, Xp, Yp)

alf = 5.
tGF_md  = np.nanmedian(TGF, axis=1)
tGF_lp  = np.nanpercentile(TGF, alf, axis=1)
tGF_up  = np.nanpercentile(TGF, (100.-alf), axis=1)
tGL_md  = np.nanmedian(TGL, axis=1)
tGL_lp  = np.nanpercentile(TGL, alf, axis=1)
tGL_up  = np.nanpercentile(TGL, (100.-alf), axis=1)
tWOA_md = np.nanmedian(TWOA, axis=1)
tWOA_lp = np.nanpercentile(TWOA, alf, axis=1)
tWOA_up = np.nanpercentile(TWOA, (100.-alf), axis=1)

sGF_md  = np.nanmedian(SGF, axis=1)
sGF_lp  = np.nanpercentile(SGF, alf, axis=1)
sGF_up  = np.nanpercentile(SGF, (100.-alf), axis=1)
sGL_md  = np.nanmedian(SGL, axis=1)
sGL_lp  = np.nanpercentile(SGL, alf, axis=1)
sGL_up  = np.nanpercentile(SGL, (100.-alf), axis=1)
sWOA_md = np.nanmedian(SWOA, axis=1)
sWOA_lp = np.nanpercentile(SWOA, alf, axis=1)
sWOA_up = np.nanpercentile(SWOA, (100.-alf), axis=1)

clrgf   = [0., 0.4, 0.8] 
clrpgf  = [0.9, 0.95, 1]
clrgl   = [0.9, 0.3, 0]
clrpgl  = [1, 0.95, 0.9]
clrwoa  = [0., 0.9, 0.7]
clrpwoa = [0.9, 1, 0.92]

# Approximate limits:
iz = min(np.where(ZWOA <= hbtm)[0])
tlim1 = np.round(np.nanmin(tWOA_lp[:iz]), decimals=1)-0.1
tlim2 = np.round(np.nanmax(tWOA_up[:iz]), decimals=1)
slim1 = np.round(np.nanmin(sWOA_lp[:iz]), decimals=1)-0.1
slim2 = np.round(np.nanmax(sWOA_up[:iz]), decimals=1)

# Polygons with percentiles:
innan   = max(np.where(~np.isnan(tWOA_lp))[0])
tWOA_lp = np.where(np.isnan(tWOA_lp), tWOA_lp[innan], tWOA_lp)
innan   = max(np.where(~np.isnan(tWOA_up))[0])
tWOA_up = np.where(np.isnan(tWOA_up), tWOA_up[innan], tWOA_up)
innan   = max(np.where(~np.isnan(sWOA_lp))[0])
sWOA_lp = np.where(np.isnan(sWOA_lp), sWOA_lp[innan], sWOA_lp)
innan   = max(np.where(~np.isnan(sWOA_up))[0])
sWOA_up = np.where(np.isnan(sWOA_up), sWOA_up[innan], sWOA_up)
tprct_y = np.concatenate((ZWOA, np.flip(ZWOA)))
tprct_x = np.concatenate((tWOA_lp, np.flip(tWOA_up)))
sprct_y = np.concatenate((ZWOA, np.flip(ZWOA)))
sprct_x = np.concatenate((sWOA_lp, np.flip(sWOA_up)))


#==================
plt.ion()

fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
plt.clf()

# T profile
ax1 = plt.axes([0.08, 0.1, 0.25, 0.8])
ln1, = ax1.plot(tGF_md,  ZMGF0, color=clrgf,  linewidth=2, label="GOFSrst")
ln2, = ax1.plot(tGL_md,  ZMGL0, color=clrgl,  linewidth=2, label="GLORYSrst")
ln3, = ax1.plot(tWOA_md, ZWOA,  color=clrwoa, linewidth=2, label="WOA23")
ln4, = ax1.plot(tWOA_lp, ZWOA, color=clrpwoa, linewidth=2, label="5-95 IPR WOA")
ax1.plot(tWOA_up, ZWOA, color=clrpwoa, linewidth=2)
ax1.fill(tprct_x, tprct_y, facecolor=clrpwoa)

plt_prct = False
if plt_prct:
  ax1.plot(tGF_lp, ZMGF0, color=clrpgf, linewidth=2)
  ax1.plot(tGF_up, ZMGF0, color=clrpgf, linewidth=2)
  ax1.plot(tGL_lp, ZMGL0, color=clrpgl, linewidth=2)
  ax1.plot(tGL_up, ZMGL0, color=clrpgl, linewidth=2)


ax1.set_ylim([hbtm, 0])
ax1.set_xlim([tlim1, tlim2])
ax1.grid(True)
ctl1 = f"Median T {regn} {YR}/{MM}/{DD}"
ax1.set_title(ctl1)

# S profile
ax2 = plt.axes([0.4, 0.1, 0.25, 0.8])
ax2.plot(sGF_md,  ZMGF0, color=clrgf,  linewidth=2)
ax2.plot(sGL_md,  ZMGF0, color=clrgl,  linewidth=2)
ax2.plot(sWOA_md, ZWOA,  color=clrwoa, linewidth=2)
ax2.plot(sWOA_lp, ZWOA, color=clrpwoa, linewidth=2)
ax2.plot(sWOA_up, ZWOA, color=clrpwoa, linewidth=2)
ax2.fill(sprct_x, sprct_y, facecolor=clrpwoa)


ax2.set_ylim([hbtm, 0])
ax2.set_xlim([slim1, slim2])
ax2.grid(True)
ctl2 = f"Median S {regn}"
ax2.set_title(ctl2)

# Plot map:
LMSK = HHM.copy()
LMSK = np.where(LMSK<0.,0.,1.)
lcmp = mutil.clrmp_lmask(2, clr_land=[0.2,0.2,0.2])
ax3 = plt.axes([0.65, 0.4, 0.3, 0.3])
ax3.pcolormesh(LMSK, shading='flat',\
                cmap=lcmp, \
                vmin=0, \
                vmax=1)
ax3.axis('scaled')
ax3.contour(HHM, [-5000,-4000,-3000,-2000,-1000], \
            colors=[(0.9,0.9,0.9)], linestyles='solid',linewidths=1)
ax3.plot(IP, JP, marker='.', color=[0.,0.6,1])

# Set y, x axis not visible
#ahd = plt.gca()
xax = ax3.axes.get_xaxis()
xax = xax.set_visible(False)
yax = ax3.axes.get_yaxis()
yax = yax.set_visible(False)
ax3.set_title('Selected profiles')

ax5 = plt.axes([0.7, 0.1, 0.25, 0.25])
lgd = plt.legend(handles=[ln1,ln2,ln3,ln4], loc='upper left')
#ax5.text(0,0.02,ss2)
ax5.axis('off')


#fig1.colorbar(im,ax=ax1,orientation='horizontal')
btx = 'plot_TSprof_WOA.py'
bottom_text(btx, pos=[0.02, 0.02])





