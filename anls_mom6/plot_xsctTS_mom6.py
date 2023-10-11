"""
  Plot sections with vertical distribution of T or S
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
from copy import copy
import matplotlib.colors as colors
from matplotlib.patches import Polygon

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_cice6_utils as mc6util
#import mod_valid_utils as mvutil


expt    = '003'
YR      = 2020
jday    = 1
HR      = 12
hg      = 1.e15
isct    = [0]  # specify sections to plot
fld     = 'potT'  # salt or potT

XSCT = mutil.rtofs_sections()
SCTnames = list(XSCT.keys())

pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'

dnmb = mtime.jday2dnmb(YR,jday)
DV   = mtime.datevec(dnmb)
MM   = DV[1]
DD   = DV[2]

pthbin = pthrun + 'tarmom_{0}{1:02d}/'.format(YR,MM)
flin   = pthbin + 'ocnm_{0}_{1:03d}_{2}.nc'.format(YR,jday,HR)


import mod_mom6 as mom6util
importlib.reload(mom6util)
pthgrid  = pthrun + 'INPUT/'
fgrd_mom = pthgrid + 'regional.mom6.nc'
ftopo_mom= pthgrid + 'ocean_topog.nc'
LON, LAT = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
HH       = mom6util.read_mom6depth(ftopo_mom)
Lmsk     = mom6util.read_mom6lmask(ftopo_mom)
idma     = LON.shape[0]
jdma     = LAT.shape[0]
ijdma    = idma*jdma

A2D = mom6util.read_mom6(flin, fld, rLayer=5)

