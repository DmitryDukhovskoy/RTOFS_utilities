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

YRM   = 2021
YRR   = 2023

fld2d = 'salt'
#fld2d = 'potT'
nrun  = 'GOFS3.1'  # MOM6, RTOFS, GOFS3.1

#sctnm = 'Fram79s2'
#sctnm = 'DavisS2'
#sctnm = 'Yucatan2'  # slanted section
#sctnm = 'BarentsS'
#sctnm = 'BeringS'
#sctnm = 'DenmarkS'
#sctnm = 'IclShtl'
#sctnm = 'ShtlScot'
#sctnm = 'LaManch'
#sctnm = 'NAtl39'
#======= Ocean Sections =====
#sctnm = 'BaffNAFram'
sctnm = 'AlaskaIcld' 


if nrun == 'MOM6':
  expt = '003'
  YR   = YRM
elif nrun == 'RTOFS':
  expt = 'product' # 003 or product
  YR   = YRR 
elif nrun == 'GOFS3.1':
  expt = '93.0'
  YR   = YRM

# Years for MOM6 - 2021, 2020 not extracted
# RTOFS - 2023, 2022 
dnmb1 = mtime.datenum([YR,1,1])
dnmb2 = mtime.datenum([YR,12,31])

dv1   = mtime.datevec(dnmb1)
dv2   = mtime.datevec(dnmb2)
hg    = 1.e15


print(f"\n Plotting {nrun}-{expt} {sctnm} {fld2d} \n")

import mod_misc1 as mmisc
import mod_mom6 as mom6util
import mod_colormaps as mcmp
import mod_read_hycom as mhycom
importlib.reload(mom6util)
importlib.reload(mmisc)
importlib.reload(mcmp)

if nrun == 'MOM6':
  pthrun  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
            '008mom6cice6_' + expt + '/'
  pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/' + \
            'MOM6_CICE6/expt{0}/'.format(expt)
  floutp  = f"mom6-{expt}_{fld2d}VFlx_{dv1[0]}" + \
            f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
  pthgrid   = pthrun + 'INPUT/'
  fgrd_mom  = pthgrid + 'regional.mom6.nc'
  ftopo_mom = pthgrid + 'ocean_topog.nc'
  HH  = mom6util.read_mom6depth(ftopo_mom)
elif nrun == 'RTOFS':
  if fld2d == 'salt':
    fld = 'salin'
  elif fld2d == 'potT':
    fld = 'temp'

  pthrun  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/wcoss2.prod/'
  pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/RTOFS_production/'
  floutp  = f"rtofs-{expt}_{fld}xsct_{dv1[0]}" + \
            f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
  pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
  ftopo   = 'regional.depth'
  fgrid   = 'regional.grid'
  _, _, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
elif nrun == 'GOFS3.1':
  if fld2d == 'salt':
    fld = 'saln'
  elif fld2d == 'potT':
    fld = 'temp'

  pthrun  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/GOFS3.1/restart/'
  pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/GOFS3.1/'
  floutp  = f"gofs31-930_{fld}xsct_{dv1[0]}" + \
            f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
  pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
  ftopo   = 'regional.depth'
  fgrid   = 'regional.grid'
  _, _, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)


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

# Time-mean section:
A2d = np.nanmean(A2dT, axis=0).squeeze()
A2d[JN,IN] = np.nan
A2di = A2d.copy()

# For plotting - fill land/bottom and smooth 2D field:
#A2di = mom6util.fill_bottom(A2di, ZZi, Hbtm)
A2di = mom6vld.box_fltr(A2di, npnts=3)
#A2di = mom6util.bottom2nan(A2di, ZZi, Hbtm)


STR = mom6vld.ocean_straits()
if sctnm in STR:
  oc_strait = True
else:
  STR = mom6vld.ocean_sections()
  if sctnm in STR:
    oc_strait = False
  else:
    raise Exception(f"Name {sctnm} is not defined as a strait or section")

IJ = np.zeros((len(II),2))
IJ[:,0] = II
IJ[:,1] = JJ

btx = 'plot_TSsection.py'
plt.ion()

import mod_utils as mutil
import plot_sect as psct
#importlib.reload(mcmp)

if fld2d == 'salt':
#  cmpr = mcmp.colormap_salin(clr_ramp=[0.94,0.95,1])
  cmpr = mcmp.colormap_haline()
  rmin = STR[sctnm]["smin"]
  rmax = STR[sctnm]["smax"] 
  cntr1 = STR[sctnm]["scntr"]
elif fld2d == 'potT':
#  cmpr = mcmp.colormap_temp()
  cmpr = mcmp.colormap_temp2()
  rmin = STR[sctnm]["tmin"]
  rmax = STR[sctnm]["tmax"] 
  cntr1 = STR[sctnm]["tcntr"]


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

ni = len(II)
XI = np.arange(0, ni, 1, dtype=int)
XI = np.cumsum(Lsgm)*1.e-3 # distance along section, km
mom6vld.plot_xsect(XI, Hbtm, ZZi, A2di, HH, stl=stl,\
                   rmin = rmin, rmax = rmax, clrmp=cmpr,\
                   IJs=IJ, btx=btx, btm_midpnt=True, cntr2=cntr1)


