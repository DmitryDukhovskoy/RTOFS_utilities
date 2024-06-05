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


pthout = '/work/Dmitry.Dukhovskoy/data/mom6_nep_restart/'
#pthout = '/work/Dmitry.Dukhovskoy/NEP_input/restart/'

btx = 'check_restart.nc'


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


# Read MOM6 NEP domain nx=342, ny=816
with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom)
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom)
hlon, hlat  = mom6util.read_mom6grid(dfgrid_mom, grid='symmetr', grdpnt='hgrid')
HHM         = mom6util.read_mom6depth(dftopo_mom)
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]

#flnm = 'MOM.res.NEP_20090401.nc'
flnm = 'MOM.res.GOFS3.0_1993010100.nc'
fnc  = os.path.join(pthout, flnm)
nc   = ncFile(fnc)

vrnm = 'Salt'
A3d   = read_field(fnc, vrnm)
A3d   = A3d.squeeze()
rmin = 5.
rmax = 35.
#cmpr = mutil.colormap_ssh(nclrs=200)
cmpr = mutil.colormap_salin2(nclrs=200)


k=0
A2d = A3d[k,:,:].squeeze()

plt.ion()
fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
plt.clf()

ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
im1 = ax1.pcolormesh(A2d, \
                 cmap=cmpr,\
                 vmin=rmin, \
                 vmax=rmax)
ax1.contour(HHM,[(0)], colors=[(0,0,0)])

ax1.axis('scaled')
ax1.set_xlim([0, idm-1])
ax1.set_ylim([0, jdm-1])

ctlt = f'flnm {vrnm} Lr={k+1}'
ax1.set_title(ctlt)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')







