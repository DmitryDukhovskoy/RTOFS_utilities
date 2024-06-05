"""
  Quick plot of saved fields in ocean_daily.nc
"""
import os
import numpy as np
import matplotlib.pyplot as plt
#import torch
import sys
import pdb
import netCDF4
import importlib
from netCDF4 import Dataset as ncFile
import timeit
import pickle
import yaml


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
import mod_mom6 as mom6util
import mod_misc1 as mmsc1

f_cont = True  # continue from last saved record

grd=0.25
if grd==0.25:
  cgrd=4
woa='woa23'
seas=15    # season: 1-12 monthly, 13-winter (Jan-Mar), 14-spring (Apr-Jun), ...

pthout = '/work/Dmitry.Dukhovskoy/run_output/NEP_nudged/'

btx = 'read_oceandaily.nc'

tfnm=f'{woa}_decav_t{seas:02d}_{cgrd:02d}.nc'
sfnm=f'{woa}_decav_s{seas:02d}_{cgrd:02d}.nc'

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

flnm = '20040101.20040101.ocean_daily.nc'
fnc  = os.path.join(pthout, flnm)
nc   = ncFile(fnc)

vrnm = 'ssh'
ssh  = read_field(fnc, vrnm)
AA  = ssh.squeeze()
rmin = -0.2
rmax = 0.2
cmpr = mutil.colormap_ssh(nclrs=200)


plt.ion()
fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
plt.clf()

ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
im1 = ax1.pcolormesh(AA, \
                 cmap=cmpr,\
                 vmin=rmin, \
                 vmax=rmax)

ax1.axis('scaled')
ax1.set_xlim([0, idm-1])
ax1.set_ylim([0, jdm-1])

ctlt = f'flnm {vrnm}'
ax1.set_title(ctlt)




