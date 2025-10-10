import numpy as np
import os
import importlib
import xarray
import matplotlib.pyplot as plt
#from yaml import safe_load
from pathlib import Path, PurePath

pthob = '/gpfs/f5/cefi/world-shared/NEP_input/obcs/glorys/'
fobc  = os.path.join(pthob,'nep_10km_glorys_obcs_1993.nc')

dset = xarray.open_dataset(fobc)

# On a supergrid
# Segment 1
it   = 10
isgm = 4
dz1  = dset[f'dz_so_segment_{isgm:03d}'].isel(time=it).data.squeeze()
lon1 = dset[f'lon_segment_{isgm:03d}'].data
lat1 = dset[f'lat_segment_{isgm:03d}'].data
xx1  = np.arange(0,len(lon1))


varnm = 'thetao'
zz1 = -np.cumsum(dz1, axis=0)  # note these are interface depths, need to add 0 at the top
A1  = dset[f'{varnm}_segment_{isgm:03d}'].isel(time=it).data.squeeze()


plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.24, 0.85, 0.7])
rmin = -2.
rmax = 20.
im1 = ax1.pcolormesh(xx1,zz1,A1, \
                 vmin=rmin, \
                 vmax=rmax)

plt.colorbar(im1)
ax1.set_title(f'{varnm} segm={isgm} min/max lat={min(lat1):4.1f}/{max(lat1):4.1f}' +
              f' lon={min(lon1):4.1f}/{max(lon1):4.1f}')

#ax1.axis('scaled')
#ax1.set_xlim([1010, 2280])
#ax1.set_ylim([1500, 3050])



