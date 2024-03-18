# Calculate 1st baroclinic Rossby R and corresponding phase speed 
# Using WKB approach described in Chelton 1996
import os
import numpy as np
import matplotlib.pyplot as plt
import torch
import sys
import netCDF4
from netCDF4 import Dataset as ncFile

# sys.path.append('

grd=0.25
#ftmp='/data/ncei2/w18C/analysis/all_0/{0}/mean/M02001'.format(grd)
#fsal='/data/ncei2/w18C/analysis/all_0/{0}/mean/s013'.format(grd)

#nlat=720
#nlon=1440
#nz=100
#dt=np.dtype((np.float32,(nlat,nlon)))

urlT='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/temperature/decav/{0}/'.\
      format(grd)
urlS='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/salinity/decav/{0}/'.\
      format(grd)

tfnm='woa18_decav_t13_04.nc'
sfnm='woa18_decav_s16_04.nc'

def read_field(furl,varnm):
  print("Reading {0}".format(furl))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].squeeze()
  dmm = np.copy(dmm0)
  return dmm

def lookup_ncvar(nc):
  ii=0
  for var in nc.variables.values():
    ii+=1
    print('--------\n')
    print('Var # {0}'.format(ii))
    print(var)


iz=0 
tvar='t_an'
svar='s_an'

furl=urlT+tfnm
T=read_field(furl,tvar)
T[T>1.e10]=np.nan
kdm=T.shape[0]
idm=T.shape[1]
jdm=T.shape[2]

ZZ=read_field(furl,'depth')
ZZ=-abs(ZZ)
LAT=read_field(furl,'lat')
LON=read_field(furl,'lon')

furlS=urlS+sfnm
S=read_field(furlS,svar)
S[S>1.e10]=np.nan

aa=T[iz,:,:].squeeze()
ss=S[iz,:,:].squeeze()


plt.ion()
fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
im = ax1.pcolormesh(aa,shading='flat')
ax1.axis('scaled')
ax1.set(xlim=(0,idm),ylim=(0,jdm))
ax1.autoscale(enable=True, axis='both', tight=True)

ax1.set_title('Var {0}, z level={1}'.format(tvar,iz))
fig1.colorbar(im,ax=ax1,orientation='horizontal')


plt.show()

cmpS = plt.cm.get_cmap('turbo')
fig2 = plt.figure(2,figsize=(8,8), constrained_layout=False)
plt.clf()
ax2 = plt.axes([0.1, 0.1, 0.8, 0.8])
im2= ax2.pcolormesh(ss,cmap=cmpS,shading='flat',vmin=30.,vmax=36.)
ax2.axis('scaled')
ax2.set(xlim=(0,idm),ylim=(0,jdm))
ax2.autoscale(enable=True, axis='both', tight=True)

ax2.set_title('Var {0}, z level={1}'.format(svar,iz))
fig2.colorbar(im2,ax=ax2,orientation='horizontal')

plt.show()




