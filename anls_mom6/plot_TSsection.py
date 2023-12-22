"""
  Plot extracted 2D field along a section
  zig-zagging sections are allowed

  extracted in extrTSxsect_polysegm.py
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
#sctnm = 'DavisStr'
sctnm = 'Yucatan2'  # slented section

fld2d = 'salt'
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
  F2D = pickle.load(fid)

TM   = F2D.TM
A2dT = F2D.Fld2D  # 2D fields: Time x depth x Width
II   = F2D.Iindx
JJ   = F2D.Jindx
XX   = F2D.LON
YY   = F2D.LAT
Lsgm = F2D.Lsgm
Hbtm = F2D.Hbtm
ZZi  = F2D.ZZi

Xdst = np.cumsum(Lsgm)
Xdst = np.insert(Xdst,0,0)
#mtime.datestr(TM)

KN,JN,IN = np.where(np.isnan(A2dT))
f_nans = False
if len(KN) >  0:
  f_nans = True
  A2dT = np.where(np.isnan(A2dT), 1.e30, A2dT)

pthgrid   = pthrun + 'INPUT/'
fgrd_mom  = pthgrid + 'regional.mom6.nc'
ftopo_mom = pthgrid + 'ocean_topog.nc'

import mod_misc1 as mmisc
import mod_mom6 as mom6util
HH  = mom6util.read_mom6depth(ftopo_mom)

# Time-mean section:
A2d = np.nanmean(A2dT, axis=0).squeeze()
A2d[JN,IN] = 0.0

# For plotting - project slanted sections on
# X or Y axis
# Interpolate over the gaps for smooth picture
#f_proj = 'X' 
#if f_proj == 'X' or f_proj == 'x':
#  A2di = mom6vld.project2X(II,JJ,A2d)
#elif f_proj == 'Y' or f_proj == 'y':
#  A2di = mom6vld.project2Y(II,JJ,A2d)
A2di = A2d.copy()

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

ni = len(II)
#XI = np.arange(0, ni, 1, dtype=int)
XI = np.cumsum(Lsgm)*1.e-3 # distance along section, km
mom6vld.plot_xsect(XI, Hbtm, ZZi, A2di, HH, stl=stl,\
                   rmin = rmin, rmax = rmax, clrmp=cmpr,\
                   IJs=IJ, btx=btx)


