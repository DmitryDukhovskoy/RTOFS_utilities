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

YRM   = 2021
YRR   = 2023
nrun  = 'MOM6'  # MOM6, RTOFS, GOFS3.1
#sctnm = 'Fram79'
#sctnm = 'Fram79s2'
#sctnm = 'DavisStr'
#sctnm  = 'DavisStr2' # straight line
#sctnm = 'DavisS2'   # slanted section
#sctnm = 'Yucatan2'  # slanted section
sctnm = 'FlorCabl'
#sctnm = 'BarentsS'
#sctnm = 'BeringS'
#sctnm = 'DenmarkS'
#sctnm = 'IclShtl'
#sctnm = 'ShtlScot'
#sctnm = 'LaManch'
#sctnm = 'NAtl39'
fld2d = 'Unrm'

mS=1
dS=1
if nrun == 'MOM6':
  expt = '003'
  YR   = YRM
elif nrun == 'RTOFS':
  expt = 'product' # 003 or product
  YR   = YRR
#mS=9
#dS=15
elif nrun == 'GOFS3.1':
  expt = '93.0'
  YR   = YRM

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

hg    = 1.e15

if nrun == 'MOM6':
  pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
            '008mom6cice6_' + expt + '/'
  pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/' + \
             'MOM6_CICE6/expt{0}/'.format(expt)
  floutp = f"mom6-{expt}_{fld2d}VFlx_{dv1[0]}" + \
            f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
  pthgrid   = pthrun + 'INPUT/'
  fgrd_mom  = pthgrid + 'regional.mom6.nc'
  ftopo_mom = pthgrid + 'ocean_topog.nc'
  HH  = mom6util.read_mom6depth(ftopo_mom)
elif nrun == 'RTOFS':   
  pthrun  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/wcoss2.prod/'
  pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/RTOFS_production/'
  floutp  = f"rtofs-{expt}_{fld2d}xsct_{dv1[0]}" + \
            f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
  pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
  ftopo   = 'regional.depth'
  fgrid   = 'regional.grid'
  _, _, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
elif nrun == 'GOFS3.1':   
  pthrun  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/GOFS3.1/restart/'
  pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/GOFS3.1/'
  floutp  = f"gofs31-930_{fld2d}xsct_{dv1[0]}" + \
            f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
  pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
  ftopo   = 'regional.depth'
  fgrid   = 'regional.grid'
  _, _, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

ffout = pthoutp + floutp
print('Loading ' + ffout)

with open(ffout, 'rb') as fid:
  F2D, UFLX = pickle.load(fid)

# 2D fields are at half-grid points
TM   = F2D.TM
Unrm = F2D.Fld2D  # 2D Flow: Time x depth x Width
II   = F2D.Iindx
JJ   = F2D.Jindx
XX   = F2D.LON
YY   = F2D.LAT
Lsgm = F2D.Lsgm
Hbtm = F2D.Hbtm
ZZi  = F2D.ZZi
# Depth-integrated: full grid
VFlx = UFLX.trnsp*1e-6  # 1D depth-integrated flow, Sv
#Hbtm = UFLX.Hbtm
Xdst = np.cumsum(Lsgm)
#Xdst = np.insert(Xdst,0,0)
#mtime.datestr(TM)

KN,JN,IN = np.where(np.isnan(Unrm))
f_nans = False
if len(KN) >  0:
  f_nans = True
  Unrm = np.where(np.isnan(Unrm), 1.e30, Unrm)

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


# Check segments etc:
# Draw section line with all segments, norm vectors and
# half-segments
f_chcksgm = False
if f_chcksgm:
  LON, LAT  = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
  DX, DY    = mom6util.dx_dy(LON,LAT)
  SGMT      = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive')
  Vnrm1     = SGMT.half_norm1
  Vnrm2     = SGMT.half_norm2
  fgnmb=2
  ax2 = mom6vld.plot_section_map(II, JJ, IJ, Vnrm1, Vnrm2, IIhf, JJhf, \
                     HH, fgnmb=fgnmb, btx='plot_UVsection.py')
  A = stop

# Time-mean section:
Fav = np.nanmean(Unrm, axis=0).squeeze()
Fav[JN,IN] = 0.0

# For plotting - U is already projected onto unit normal 
# of the main section(s)
Favi = Fav.copy()

# Truncate 2D field to bottom depth
# to mask out filled land values 
#Favi = mom6util.fill_bottom(Favi, ZZi, Hbhf)

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
stl = f"0.08 {nrun}-{expt} {sctnm}, {fld2d}  Mean: " + \
      f"{YR1}/{MM1:02d}/{DD1:02d}-{YR2}/{MM2:02d}/{DD2:02d}"

# For plotting - smooth 2D field:
#Favi = mom6vld.box_fltr(Favi, npnts=3)

# Contour:
cntr1 = STR[sctnm]["ucntr1"] 
cntr2 = STR[sctnm]["ucntr2"]
#cntr1=[-0.2,-0.15,-0.1,-0.05]
#cntr2=[0,0.25,0.5,0.75,1.,1.25,1.5,1.75]

ni = len(II)
#XI = np.arange(0, ni, 1, dtype=int)
XI = (np.cumsum(Lsgm) - Lsgm[0])*1.e-3# distance along section, km
mom6vld.plot_xsect(XI, Hbtm, ZZi, Favi, HH, stl=stl,\
                   rmin = rmin, rmax = rmax, clrmp=cmpr,\
                   IJs=np.vstack((II,JJ)).T, btx=btx, cntr1=cntr1, cntr2=cntr2)

# Plot depth-integrated transport
Tday = TM-TM[0]

#VFlxi = VFlx.copy()
# Mean flux for plotting:
# Interpolate over zigzagging segments
# and smooth for plotting
VFlxi = mom6vld.project2X(II,JJ,VFlx) 
for ii in range(3):
  VFlxi = mom6vld.runmn1D(VFlxi, npnts=3, axis=0)

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

#XXI = np.arange(0,nnI,1)
XXI = (np.cumsum(Lsgm)-Lsgm[0])*1.e-3 # distance along section, km
#XXI = Lon1


# Depth-integrated transport

plt.ion()
fig2 = plt.figure(2,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
ax1.plot(XXI, mVF)
ax1.plot(XXI, pL, color=[0.8, 0.85, 1.])
ax1.plot(XXI, pU, color=[0.8, 0.85, 1.])

f_setrgn = False
if f_setrgn:
# Bind the button_press_event with the onclick() method
  fig2.canvas.mpl_connect('button_press_event', onclick)

xl1 = np.min(XXI)
xl2 = np.max(XXI)

ax1.set_xlim([xl1, xl2])
ax1.set_xticks(np.linspace(np.floor(xl1),np.ceil(xl2),10))
ax1.grid(True)

ctl = f'{nrun}-{expt} {sctnm} depth.intgr. VFlux (Sv), \n' + \
      f'mean & IDR {dv1[0]}/{dv1[1]}/{dv1[2]} - {dv2[0]}/{dv2[1]}/{dv2[2]}'
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



f_topo = False
if f_topo:
  fig1 = plt.figure(1, figsize(9,8))
  plt.clf()
  ax11 = plt.axes([0.1, 0.1, 0.8, 0.8])
  ax11.contour(HH,[0.],colors=[(0,0.,0.)], linestyles='solid')
  ax11.contour(HH,[-500.],colors=[(0,0.6,0.9)], linestyles='solid')





