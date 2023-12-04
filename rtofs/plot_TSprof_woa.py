"""
  Plot T/S profiles at a given locations
  WOA18 climatologies

  Cara Manning
  Baffin Bay T/S anomaly in central Baffing 2017-2019

The station locations are:
BB2: 72.75 N, 67 W
224: 70.43 N, 62.96 W

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
from copy import copy

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
import mod_read_ncoda as rncoda
from mod_utils_fig import bottom_text

MM = 8
latPrf = 72.75
lonPrf = -67.

# -sea* files: land filled with S/T values - for interpolation
pthwoa  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/WOA18/'
flS     = f"woa18_025grd_saln{MM:02d}-sea.dat"
flSstd  = f"woa18_025grd_sstd{MM:02d}-sea.dat"
flT     = f"woa18_025grd_temp{MM:02d}-sea.dat"
flTstd  = f"woa18_025grd_tstd{MM:02d}-sea.dat"
flGrid  = 'woa18_depth_gridFort.dat' 
fld_saln = os.path.join(pthwoa, flS)
fld_Sstd = os.path.join(pthwoa, flSstd)
fld_temp = os.path.join(pthwoa, flT)
fld_Tstd = os.path.join(pthwoa, flTstd)
fld_Grid = os.path.join(pthwoa, flGrid)

import mod_WODdata as mwod
importlib.reload(mwod)
kdm, idm, jdm, ZZ, LON, LAT = mwod.read_WOA18_grid(fld_Grid)
#S3d = mwod.read_WOA18_3D(fld_saln, idm, jdm, kdm)
Swoa = mwod.woa_profile(fld_saln, LON, LAT, idm, jdm, kdm, lonPrf, latPrf)
Sstd = mwod.woa_profile(fld_Sstd, LON, LAT, idm, jdm, kdm, lonPrf, latPrf)
Twoa = mwod.woa_profile(fld_temp, LON, LAT, idm, jdm, kdm, lonPrf, latPrf)
Tstd = mwod.woa_profile(fld_Tstd, LON, LAT, idm, jdm, kdm, lonPrf, latPrf)

import plot_sect as psct
importlib.reload(psct)

fgnmb=1
plt.ion()
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.35, 0.8])
ax2 = plt.axes([0.55, 0.1, 0.35, 0.8])

stl = f"WOA18 T MM={MM} Lon={lonPrf:.2f} Lat={latPrf:.2f}"
xl1=-2
xl2=5.5
psct.plot_profile(ax1, Twoa, ZZ, 10., stl=stl, zlim=-2100, \
                  xl1=xl1, xl2=xl2)
psct.plot_profile(ax1, Twoa-Tstd, ZZ, 10., stl=stl, zlim=-2100, \
                  xl1=xl1, xl2=xl2, lclr=(0.7, 0.7, 0.7), f_line=True)
psct.plot_profile(ax1, Twoa+Tstd, ZZ, 10., stl=stl, zlim=-2100, \
                  xl1=xl1, xl2=xl2, lclr=(0.7, 0.7, 0.7), f_line=True)


stl = f"WOA18 S MM={MM} Lon={lonPrf:.2f} Lat={latPrf:.2f}"
sxl1=29.5
sxl2=34.95
psct.plot_profile(ax2, Swoa, ZZ, 10., stl=stl, zlim=-2100, \
                  xl1=sxl1, xl2=sxl2)
psct.plot_profile(ax2, Swoa-Sstd, ZZ, 10., stl=stl, zlim=-2100, \
                  xl1=sxl1, xl2=sxl2, lclr=(0.7, 0.7, 0.7), f_line=True)
psct.plot_profile(ax2, Swoa+Sstd, ZZ, 10., stl=stl, zlim=-2100, \
                  xl1=sxl1, xl2=sxl2, lclr=(0.7, 0.7, 0.7), f_line=True)

btx = 'plot_TSprof_woa.py'
bottom_text(btx)



