"""
  Plot extracted 2D field along a section
  Note variables are on staggered grid in the output files

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

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/mom6_utils')

import mod_time as mtime
from mod_utils_fig import bottom_text

import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)

expt  = '003'
hg    = 1.e15
#sctnm = 'Fram79'
sctnm = 'DavisS2'
#sctnm = 'Yucatan2'  # slented section
fld2d = 'Unrm'

#fld2d = 'salt'
#fld2d = 'potT'
dnmb1 = mtime.datenum([2021,1,1])
dnmb2 = mtime.datenum([2021,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)

pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'
pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/' + \
          'MOM6_CICE6/expt{0}/'.format(expt)
#floutp = 'mom6-{4}_u2dsect_{0}{1:02d}-{2}{3:02d}_{5}.pkl'.\
#         format(dv1[0], dv1[1], dv2[0], dv2[1], expt, sctnm)
floutp = f"mom6-{expt}_{fld2d}VFlx_{dv1[0]}" + \
         f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
ffout = pthoutp + floutp

print('Loading ' + ffout)

with open(ffout, 'rb') as fid:
  F2D, UFLX = pickle.load(fid)

# 2D fields are at half-grid points
TM   = F2D.TM
Unrm = F2D.Fld2D  # 2D Flow: Time x depth x Width
IIhf   = F2D.Iindx
JJhf   = F2D.Jindx
XX   = F2D.LON
YY   = F2D.LAT
Lsgm = F2D.Lsgm
Hbhf = F2D.Hbtm
ZZi  = F2D.ZZi
# Depth-integrated: full grid
VFlx = UFLX.trnsp*1e-6  # 1D depth-integrated flow, Sv
Hbtm = UFLX.Hbtm
Xdst = np.cumsum(Lsgm)
#Xdst = np.insert(Xdst,0,0)
#mtime.datestr(TM)

KN,JN,IN = np.where(np.isnan(Unrm))
f_nans = False
if len(KN) >  0:
  f_nans = True
  Unrm = np.where(np.isnan(Unrm), 1.e30, Unrm)

pthgrid   = pthrun + 'INPUT/'
fgrd_mom  = pthgrid + 'regional.mom6.nc'
ftopo_mom = pthgrid + 'ocean_topog.nc'

import mod_misc1 as mmisc
import mod_mom6 as mom6util
STR = mom6vld.ocean_straits()
nlegs = STR[sctnm]["nlegs"]
I1    = STR[sctnm]["xl1"]
I2    = STR[sctnm]["xl2"]
J1    = STR[sctnm]["yl1"]
J2    = STR[sctnm]["yl2"]
Ni    = np.zeros((nlegs))
Nj    = np.zeros((nlegs))
IJ    = np.zeros((nlegs+1,2))
for kk in range(nlegs):
  Ni[kk]     = I2[kk]+1
  Nj[kk]     = J2[kk]+1
  IJ[kk,0]   = I1[kk]
  IJ[kk,1]   = J1[kk]
  IJ[kk+1,0] = I2[kk]
  IJ[kk+1,1] = J2[kk]

HH        = mom6util.read_mom6depth(ftopo_mom)
LON, LAT  = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
DX, DY    = mom6util.dx_dy(LON,LAT)
SGMT      = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive')
II        = SGMT.I_indx
JJ        = SGMT.J_indx

# Time-mean section:
Fav = np.nanmean(Unrm, axis=0).squeeze()
Fav[JN,IN] = 0.0

# For plotting - project slanted sections on
# X or Y axis
# Interpolate over the gaps for smooth picture
f_proj = 'X' 
if f_proj == 'X' or f_proj == 'x':
  print("Projecting Unrm on X-axis for plotting")
  Favi = mom6vld.project2X(IIhf,JJhf,Fav)
elif f_proj == 'Y' or f_proj == 'y':
  Favi = mom6vld.project2Y(IIhf,JJhf,Fav)
else:
  Favi = Fav.copy()

# Truncate 2D field to bottom depth
# to mask out filled land values 
#Favi = mom6util.fill_bottom(Favi, ZZi, Hbhf)

STR = mom6vld.ocean_straits()
nlegs = STR[sctnm]["nlegs"]
I1    = STR[sctnm]["xl1"]
I2    = STR[sctnm]["xl2"]
J1    = STR[sctnm]["yl1"]
J2    = STR[sctnm]["yl2"]
Ni    = np.zeros((nlegs))
Nj    = np.zeros((nlegs))
IJ    = np.zeros((nlegs+1,2))
for kk in range(nlegs):
  Ni[kk]     = I2[kk]+1
  Nj[kk]     = J2[kk]+1
  IJ[kk,0]   = I1[kk]
  IJ[kk,1]   = J1[kk]
  IJ[kk+1,0] = I2[kk]
  IJ[kk+1,1] = J2[kk]

btx = 'plot_UVsection.py'
plt.ion()

import mod_utils as mutil
import plot_sect as psct

cmpr = mutil.colormap_ssh(nclrs=100)
if fld2d == 'salt':
  cmpr = mutil.colormap_salin(clr_ramp=[0.,0.4,0.6])
  rmin = STR[sctnm]["smin"]
  rmax = STR[sctnm]["smax"] 
elif fld2d == 'potT':
  cmpr = mutil.colormap_temp()
  rmin = STR[sctnm]["tmin"]
  rmax = STR[sctnm]["tmax"] 
elif fld2d == 'Unrm':
  cmpr = mutil.colormap_ssh(nclrs=100)
  rmin = STR[sctnm]["umin"]
  rmax = STR[sctnm]["umax"] 

dv1 = mtime.datevec(TM[0])
dv2 = mtime.datevec(TM[-1])
YR1 = dv1[0]
MM1 = dv1[1]
DD1 = dv1[2]
YR2 = dv2[0]
MM2 = dv2[1]
DD2 = dv2[2]
stl = f"0.08 MOM6-CICE6-{expt} {sctnm}, {fld2d}  Mean: " + \
      f"{YR1}/{MM1:02d}/{DD1:02d}-{YR2}/{MM2:02d}/{DD2:02d}"

# For plotting - smooth 2D field:
Favi = mom6vld.box_fltr(Favi, npnts=3)

# Contour:
cntr1=[-0.2,-0.15,-0.1,-0.05]
cntr2=[0,0.25,0.5,0.75,1.,1.25,1.5,1.75]

ni = len(IIhf)
#XI = np.arange(0, ni, 1, dtype=int)
XI = (np.cumsum(Lsgm) - Lsgm[0])*1.e-3# distance along section, km
mom6vld.plot_xsect(XI, Hbhf, ZZi, Favi, HH, stl=stl,\
                   rmin = rmin, rmax = rmax, clrmp=cmpr,\
                   IJs=IJ, btx=btx, cntr1=cntr1, cntr2=cntr2)

# Plot depth-integrated transport
Tday = TM-TM[0]

# Mean flux for plotting:
# Interpolate over zigzagging segments
# and smooth for plotting
VFlxi = mom6vld.project2X(II,JJ,VFlx) 
VFlxi = mom6vld.runmn1D(VFlxi, axis=0)

alf  = 10.
mVF  = np.nanmean(VFlxi, axis=0)
pL   = np.percentile(VFlxi, alf, axis=0)
pU   = np.percentile(VFlxi, (100.-alf), axis=0)
#
# Total transport:
VFtot = np.nansum(VFlx, axis=1)
mnVF  = np.mean(VFtot)
stdVF = np.std(VFtot)

# segment lengths for whole grid cells
# and Longitudes for whole grid cells
nnI   = len(mVF)
Lsgm1 = np.zeros((nnI))
Lon1  = np.zeros((nnI))
Lat1  = np.zeros((nnI))
for ii in range(nnI):
  ix1       = ii*2
  ix2       = ix1+1
  Lsgm1[ii] = Lsgm[ix1] + Lsgm[ix2] 
  Lon1[ii]  = 0.5*(XX[ix1] + XX[ix2])
  Lat1[ii]  = 0.5*(YY[ix1] + YY[ix2]) 
#

#XXI = np.arange(0,nnI,1)
XXI = (np.cumsum(Lsgm1)-Lsgm1[0])*1.e-3 # distance along section, km
#XXI = Lon1

plt.ion()
fig2 = plt.figure(2,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
ax1.plot(XXI, mVF)
ax1.plot(XXI, pL, color=[0.8, 0.85, 1.])
ax1.plot(XXI, pU, color=[0.8, 0.85, 1.])

xl1 = np.min(XXI)
xl2 = np.max(XXI)

ax1.set_xlim([xl1, xl2])
ax1.set_xticks(np.linspace(np.floor(xl1),np.ceil(xl2),10))
ax1.grid(True)

ctl = ('0.08 MOM6-CICE6-{0} {1}, depth.intgr. VFlux (Sv),  \n'.format(expt, sctnm)\
        + 'mean & IDR {0}/{1}/{2} - {3}/{4}/{5}'.\
       format(dv1[0], dv1[1], dv1[2], dv2[0], dv2[1], dv2[2]))
ax1.set_title(ctl)

# Plot bottom:
verts = [(xl1-1.,-8000),(xl1-1., Hbtm[0]),*zip(XXI,Hbtm),\
         (xl2+1.,Hbtm[-1]),(xl2+1.,-8000)]
poly = Polygon(verts, facecolor='0.', edgecolor='none')

ax2 = plt.axes([0.1, 0.36, 0.8, 0.1])
ax2.cla()
#ax2.plot(XX, Hbtm)
ax2.add_patch(poly)
ax2.set_xlim([xl1, xl2])
ax2.set_yticks(np.arange(-5000.,0.,500.))
ax2.set_ylim([np.floor(np.min(Hbtm)), 0])
#ax2.set_xticks(np.arange(np.floor(xl1),np.ceil(xl2),2.))
ax2.set_xticks(np.linspace(np.floor(xl1),np.ceil(xl2),10))
ax2.grid(True)

ax3 = plt.axes([0.1, 0.1, 0.6, 0.2])
ax3.plot(Tday, VFtot)
ax3.set_xlim([Tday[0], Tday[-1]])
ax3.set_xticks(np.arange(0,Tday[-1],30))
ax3.grid(True)
ax3.set_title('VolTransport, Sv')
ax3.set_xlabel('Days, {0}'.format(dv1[0]))

ax4 = plt.axes([0.72, 0.15, 0.15, 0.2])
sinfo = 'mean Flux = {0:7.2f} Sv\n'.format(mnVF)
sinfo = sinfo + 'StDev = {0:6.3f} Sv'.format(stdVF)
ax4.text(0., 0., sinfo)
ax4.axis('off')

bottom_text(btx,pos=[0.02, 0.03])






