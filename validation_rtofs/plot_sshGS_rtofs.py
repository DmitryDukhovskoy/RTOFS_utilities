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

expt      = 'product'
rdate0    = '20240101'
sfx       = 'n-24'
f_gofs    = True
f_sshanls = True

sshg      = 0.05    # Gulf Stream contour

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

# Input directories with RTOFS and 2DVAR ssh:
pthbase  = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/'
pthrtofs = f'{pthbase}rtofs.prod/rtofs.{rdate0}/'
#pthanls  = f'{pthbase}rtofs.prod/rtofs.{rdate0}/glbl_var/restart/'
pthanls  = f'{pthbase}rtofs.prod/rtofs.{rdate0}/glbl_var/restart/'
pthgofs  = f'{pthbase}gofs3.1/gofs.{rdateP}/'
pthgrid  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

huge = 1.e20
rg   = 9806.

# Define region of interest:
# For the HYCOM  0.08 Global grid !!! 
II = [2565, 2565, 2950, 2950]
JJ = [1809, 2190, 2190, 1809]

# Region to exclude shallow regions
ISH = [2562, 2569, 2578, 2587, 2590, 2609, 2721, 3012, 3012, 2562]
JSH = [1808, 1810, 1810, 1823, 1866, 1911, 1991, 2012, 1772, 1770]

frtofsa = 'rtofs_glo.t00z.' + sfx + '.archv.a'
frtofsb = 'rtofs_glo.t00z.' + sfx + '.archv.b'
dfrtofsa = os.path.join(pthrtofs,frtofsa)
dfrtofsb = os.path.join(pthrtofs,frtofsb)

SSH, _, _, _ = mhycom.read_hycom(dfrtofsa, dfrtofsb, 'srfhgt')
SSH[SSH>huge] = np.nan
SSH = SSH/9.806

IDM = SSH.shape[1]
JDM = SSH.shape[0]

ii1 = min(II)
ii2 = max(II)
jj1 = min(JJ)
jj2 = max(JJ)
SSHmn = np.nanmean(SSH[jj1:jj2,ii1:ii2])

f_dmean = True
if f_dmean:
  SSH = SSH-SSHmn


# Find Gulfstream
IGS_rtofs, JGS_rtofs = mgulf.derive_GSssh(dfrtofsa, dfrtofsb, HH, \
                                       II, JJ, ISH, JSH, SSH=SSH, sshg=0.05)

if f_gofs:
  fhcm   = f'archv.{rdateP}00'
  fina   = pthgofs + fhcm + '.a'
  finb   = pthgofs + fhcm + '.b'
  IGS_gofs, JGS_gofs = mgulf.derive_GSssh(fina,finb, HH,\
                                          II, JJ, ISH, JSH, sshg=0.05) 

if f_sshanls:
  fanls  = f'seahgt_sfc_1o3241x2441_{rdateP}00_0000_analfld'
  IGS_anls, JGS_anls = mgulf.derive_GSssh_anls(pthanls, fanls, \
                                   II, JJ, ISH, JSH, LON, LAT, HH, sshg=0.05)


# Plot fields:
print(' ===== Plotting ====')
import mod_colormaps as mclrmp
cmpr = mclrmp.colormap_ssh()
cmpr.set_bad([0,0,0])

dssh       = 0.2
ssh_cntrs  = np.arange(0,1.5,dssh)
ssh_ncntrs = np.arange(-1.2,-0.01,dssh)
clr_cntrs  = [(0.6,0.6,0.6)]
clr_ncntrs = [(0,0.4,0.5)]
clr_gs     = [0, 0, 0]    # GS contour
clr_gofs   = [0.2, 0.9, 0]
clr_anls   = [0.9, 0., 0.8]

# Gulfstream region:
xlim1 = 2520
xlim2 = 3020
ylim1 = 1760
ylim2 = 2250


plt.ion()
fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
fig1.clf()

ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
im1 = ax1.pcolormesh(SSH, cmap=cmpr, vmin=-1., vmax=1.)
ax1.set_title(f'RTOFS-product ssh & GS front {rdateP}')

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
ln1, = ax1.plot(IGS_rtofs, JGS_rtofs, color=clr_gs, label="RTOFS")
if f_gofs:
  ln2, = ax1.plot(IGS_gofs, JGS_gofs, color=clr_gofs, label="GOFS")

if f_sshanls:
  ln3, = ax1.plot(IGS_anls, JGS_anls, color=clr_anls, label="Anls")

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


ax3 = plt.axes([0.73, 0.82, 0.25, 0.18])
if f_gofs and f_sshanls:
  lgd = plt.legend(handles=[ln1,ln2,ln3], loc='upper right')
elif f_gofs:
  lgd = plt.legend(handles=[ln1,ln2], loc='upper right')
elif f_sshanls:
  lgd = plt.legend(handles=[ln1,ln3], loc='upper right')
else:
  lgd = plt.legend(handles=[ln1], loc='upper right')
ax3.axis('off')


btx = 'plot_sshGS_rtofs.py'
bottom_text(btx, pos=[0.08, 0.02])






