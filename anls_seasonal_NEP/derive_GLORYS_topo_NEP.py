"""
  Derive approximate Bottom/land masks for GLORYS reanalysis 
  Subsampled for NEP region by Liz

  Use unfilled fields, save inpython binary pickle file

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import time
import timeit
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap

#PPTHN = '/home/Dmitry.Dukhovskoy/python'
PPTHN = []
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mmom6
import xarray
from mod_misc1 import print_1col

# Region domain
lat1 = 10.
lat2 = 81.
lon1 = 156.5
lon2 = 255.5

YR1  = 2004
YR2  = 2015
HR   = 12
MM   = 1
MD   = 1
dday = 5


pthoutp   = '/work/Dmitry.Dukhovskoy/data/glorys_topo_NEP/'
pthglorys = '/archive/e1n/datasets/GLORYS/'

def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  dmm = np.copy(dmm0)
  return dmm

def lookup_ncvar(nc):
  ii=0
  for var in nc.variables.values():
    ii+=1
    print('--------\n')
    print('Var # {0}'.format(ii))
    print(var)


# Get lon/lat region boundaries:
pthinput = os.path.join(pthglorys,'2015','nep_10')
flnm = f'GLORYS_REANALYSIS_NEP_2015-12-01.nc'
dflnm = os.path.join(pthinput,flnm)
ds_glorys = xarray.open_dataset(dflnm)
S3d = ds_glorys['so'].data.squeeze()

lon = ds_glorys['longitude'].data
lat = ds_glorys['latitude'].data
Glon, Glat = np.meshgrid(lon ,lat)
jdm, idm = Glon.shape[:2]

ZM  = -ds_glorys['depth'].data
ZZ  = mmom6.zm2zz(ZM)
kdm = len(ZM)

a2d = S3d[0,:].squeeze()
HH = np.where(np.isnan(a2d),100.,-1)  # land mask at 0 m
for kk in range(1,kdm):
  am1 = S3d[kk-1,:].squeeze()
  a2d = S3d[kk,:].squeeze()  
  J,I = np.where( (np.isnan(a2d)) & (~np.isnan(am1)) )
  HH[J,I] = ZZ[kk]

jdm, idm = HH.shape
 
ftopo_regn = f"GLORYS12_topoNEP_{jdm}x{idm}.pkl"
dftopo_regn = os.path.join(pthoutp,ftopo_regn)
if not os.path.isfile(dftopo_regn): 
  print(f"Saving GLORYS region topo, grid --> {dftopo_regn}")
  with open(dftopo_regn,'wb') as fid:
    pickle.dump(HH, fid)

