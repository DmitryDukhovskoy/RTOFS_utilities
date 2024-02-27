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
sys.path.append('/home/Dmitry.Dukhovskoy/python/anls_mom6')

import mod_read_hycom as mhycom
import mod_time as mtime
from mod_utils_fig import bottom_text

import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)

expt  = 'product'
hg    = 1.e15
#sctnm = 'Fram79'
sctnm = 'DavisStr'
fld2d = 'Unrm'
#fld2d = 'salin'
#fld2d = 'temp'

dnmb1 = mtime.datenum([2021,1,1])
dnmb2 = mtime.datenum([2021,12,31])
dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)

pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/wcoss2.prod/'
pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/RTOFS_production/'
floutp = f"rtofs-{expt}_{fld2d}xsct_{dv1[0]}" + \
         f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
ffout = pthoutp + floutp

STR = mom6vld.ocean_straits()
i1  = STR[sctnm]["xl1"] 
i2  = STR[sctnm]["xl2"]
j1  = STR[sctnm]["yl1"]
j2  = STR[sctnm]["yl2"]

pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)




print('Loading ' + ffout)

with open(ffout, 'rb') as fid:
  F2D = pickle.load(fid)

FLD  = F2D.Fld2D
TM   = F2D.TM
Hbtm = F2D.Hbtm
Xsct = F2D.LON
Ysct = F2D.LAT
ZZ   = F2D.ZZi

KN,JN,IN = np.where(np.isnan(FLD))
f_nans = False
if len(KN) >  0:
  f_nans = True
  FLD = np.where(np.isnan(FLD), 1.e30, FLD)

import mod_mom6 as mom6util
import mod_utils as mutil
import plot_sect as psct

# Time-mean section:
Fav = np.nanmean(FLD, axis=0).squeeze()
Fav[JN,IN] = np.nan

btx = 'plot_2Dfld_section_rtofs.py'
plt.ion()


cmpr = mutil.colormap_ssh(nclrs=100)
if fld2d == 'salin':
  cmpr = mutil.colormap_salin(clr_ramp=[0.,0.4,0.6])
  rmin = STR[sctnm]["smin"]
  rmax = STR[sctnm]["smax"] 
elif fld2d == 'temp':
  cmpr = mutil.colormap_temp()
  rmin = STR[sctnm]["tmin"]
  rmax = STR[sctnm]["tmax"] 
elif fld2d == 'Unrm':
  cmpr = mutil.colormap_ssh(nclrs=100)
  rmin = STR[sctnm]["umin"]
  rmax = STR[sctnm]["umax"] 

YR1 = dv1[0]
MM1 = dv1[1]
DD1 = dv1[2]
YR2 = dv2[0]
MM2 = dv2[1]
DD2 = dv2[2]
stl = f"RTOFS-{expt} {sctnm}, {fld2d}  Mean: " + \
      f"{YR1}/{MM1:02d}/{DD1:02d}-{YR2}/{MM2:02d}/{DD2:02d}"

mom6vld.plot_xsect(Xsct, Hbtm, ZZ, Fav, HH, stl=stl,\
                   rmin = rmin, rmax = rmax, clrmp=cmpr,\
                   ijsct=[i1,i2,j1,j2], btx=btx)


