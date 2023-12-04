"""
  Check MOM6 increments
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

pthnc = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/GDAS/ICSDIR/C384O008/gdas.20210702/00/ocean/'
flnc  = 'MOM.inc.TSzh.nc'
foutp = os.path.join(pthnc, flnc)

import mod_misc1 as mmisc
import mod_mom6 as mom6util
importlib.reload(mom6util)

nc   = ncFile(foutp, 'r')
zh   = nc.variables['zh'][:].data
zlvl = nc.variables['Level'][:].data
 
pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_003/'
pthgrid   = pthrun + 'INPUT/'
fgrd_mom  = pthgrid + 'regional.mom6.nc'
ftopo_mom = pthgrid + 'ocean_topog.nc'
LON, LAT  = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
HH        = mom6util.read_mom6depth(ftopo_mom)
Lmsk      = mom6util.read_mom6lmask(ftopo_mom)

jj = 1300
ii = 470
dH = zh[:,jj,ii]
dzz = np.diff(zlvl)

mmisc.print_1col(dH)
mmisc.print_1col(zlvl)
mmisc.print_2col(dH,dzz)
mmisc.print_2col(dH,zlvl)


