"""
  Plot SSH and Gulfstream front using SSH from RTOFS 
  or GOFS3.1
  0.08 grid

"""
import os
import numpy as np
import sys
import importlib
import matplotlib
import matplotlib.pyplot as plt
import datetime

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hausdorff')

from mod_utils_fig import bottom_text
import mod_read_hycom as mhycom
import mod_misc1 as mmisc
import mod_hausdorff_distance as mmhd
import mod_utils_fig as mufig
import mod_utils as mutil
import mod_rtofs as mrtofs
import mod_gulfstream as mgulf
import mod_time as mtime
import mod_read_ncoda as mncoda
importlib.reload(mncoda)
importlib.reload(mmisc)
importlib.reload(mgulf)
importlib.reload(mtime)

expt    = 'gofs'  # product or gofs
rdate0  = '20240101'
sfx     = 'n-24'

# Input directories with RTOFS and 2DVAR ssh:
pthbase  = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/'
if expt == 'product':
  pthrtofs = f'{pthbase}rtofs/hycom/rtofs.{rdate0}/'
elif expt == 'gofs':
  pthrtofs = f'{pthbase}gofs3.1/gofs.{rdate0}/'
pthgrid  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

fhcm   = f'archv.{rdate0}00'
fina   = pthrtofs + fhcm + '.a'
finb   = pthrtofs + fhcm + '.b'

yrday  = mtime.rdate2jday(rdate0)
dnmb0  = mtime.rdate2datenum(rdate0)
# Date of plotted fields:
if sfx == 'n-24':
  dnmbP = dnmb0-1
  hr = 0
elif sfx[0] == 'f':
  hr = int(sfx[1:])
  dnmbP = dnmb0+float(hr)/24.
rdateP = mtime.dnumb2rdate(dnmbP, ihours=False)

dvP = mtime.datevec(dnmbP)
YRp = dvP[0]
MMp = dvP[1]
DDp = dvP[2]
HRp = dvP[3]


huge = 1.e20
rg   = 9806.
IDM, JDM, KDM = mhycom.hycom_dim(fina,finb)


ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

# Define region of interest:
# For the HYCOM  0.08 Global grid !!! 
II = [2565, 2565, 2950, 2950]
JJ = [1809, 2190, 2190, 1809]

# Region to exclude shallow regions
ISH = [2562, 2569, 2578, 2587, 2590, 2609, 2721, 3012, 3012, 2562]
JSH = [1808, 1810, 1810, 1823, 1866, 1911, 1991, 2012, 1772, 1770]


SSH, _, _, _ = mhycom.read_hycom(fina,finb,'srfhgt')
SSH[SSH>huge] = np.nan
SSH = SSH/9.806

ii1 = min(II)
ii2 = max(II)
jj1 = min(JJ)
jj2 = max(JJ)
SSHmn = np.nanmean(SSH[jj1:jj2,ii1:ii2])

f_dmean = True
if f_dmean:
  SSH = SSH-SSHmn

# Gulfstream region:
xlim1 = 2520
xlim2 = 3020
ylim1 = 1760
ylim2 = 2250

import mod_colormaps as mclrmp
cmpr = mclrmp.colormap_ssh()
cmpr.set_bad([0,0,0])

dssh       = 0.2
ssh_cntrs  = np.arange(0,1.5,dssh)
ssh_ncntrs = np.arange(-1.2,-0.01,dssh)
sshg       = 0.05    # Gulf Stream contour
clr_cntrs  = [(0.6,0.6,0.6)]
clr_ncntrs = [(0,0.4,0.5)]
clr_gs     = [0, 0, 0]    # GS contour

# Find Gulfstream
X, Y = np.meshgrid(np.arange(IDM), np.arange(JDM))
MS, IM, JM    = mmisc.inpolygon_v2(X,Y,II,JJ)  # Gulfstream region
Msh, Ish, Jsh = mmisc.inpolygon_v2(X,Y,ISH,JSH)  # shelves to exclude

z0 = -100.
SSHsub = SSH.copy()
#SSHsub[np.where( (HH<0) & (HH>=z0) & (LAT>24.2) )] = -2.
SSHsub[np.where(MS==0)] = np.nan
SSHsub[np.where((Msh==1) & (HH>z0))] = np.nan
# Patches near Bahamas
SSHsub[1808:1815, 2571:2586] = np.nan

GSCNT = mgulf.derive_contour(SSHsub, tz0=sshg)
IGS_rtofs = GSCNT[:,0]
JGS_rtofs = GSCNT[:,1]


# Plot fields:
plt.ion()
fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
fig1.clf()

ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
im1 = ax1.pcolormesh(SSH, cmap=cmpr, vmin=-1., vmax=1.)
ax1.set_title(f'{expt} RTOFS ssh {rdate0} {fhcm}')

ax1.axis('scaled')
ax1.set_xlim([xlim1,xlim2])
ax1.set_ylim([ylim1,ylim2])
f_ticksoff = True
if f_ticksoff:
  ax1.set_xticklabels([])
  ax1.set_yticklabels([])
  ax1.set_xticks([])
  ax1.set_yticks([])

ax1.contour(SSH, ssh_cntrs, linestyles='solid',
            colors=clr_cntrs, linewidths=1)
ax1.contour(SSH, ssh_ncntrs, linestyles='solid',
            colors=clr_ncntrs, linewidths=1)
ax1.plot(IGS_rtofs, JGS_rtofs, color=clr_gs)

lons = np.linspace(-180,180,73)
lats = np.linspace(-90,90,37)
XX   = np.arange(IDM)
YY   = np.arange(JDM)

ax1.contour(LON, lons, linestyles='dotted', colors=[(0.8,0.8,0.8)])
ax1.contour(LAT, lats, linestyles='dotted', colors=[(0.8,0.8,0.8)])


cax2 = fig1.add_axes([0.91, 0.2, 0.015, 0.6])
clb2 = fig1.colorbar(im1, cax=cax2, extend='both')
cax2.set_yticklabels(cax2.get_yticks())
ticklabs = clb2.ax.get_yticklabels()
#clb2.ax.set_yticklabels(ticklabs,fontsize=10)
clb2.ax.set_yticklabels(["{:.2f}".format(i) for i in clb2.get_ticks()], fontsize=10)
clb2.ax.tick_params(direction='in', length=6)

btx = 'plot_ssh_gulfstream.py'
bottom_text(btx, pos=[0.08, 0.02])






