"""
  Plot Maria domain
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
#import struct
import datetime
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
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
import mod_regmom as mrgm
#import mod_valid_utils as mvutil
importlib.reload(mcmp)
importlib.reload(mrgm)

grid_shape = 'symmetr'
grid_var = 'hgrid'

with open('paths_files_obcs.yaml') as ff:
  dct = yaml.safe_load(ff)

pthgrid_hycom = dct["RTOFS"]["pthgrid"]
ftopo_hycom   = dct["RTOFS"]["ftopo"]
fgrid_hycom   = dct["RTOFS"]["fgrid"]

# Read MOM6 domain:
pthgrid_mom  = dct["MOM6"]["pthgrid"]
ftopo_mom    = dct["MOM6"]["ftopo"]
fgrid_mom    = dct["MOM6"]["fgrid"]
pthobc_mom   = dct["MOM6"]["pthobc"]
dftopo_mom   = os.path.join(pthgrid_mom, ftopo_mom)
dfgrid_mom   = os.path.join(pthgrid_mom, fgrid_mom)
hlon, hlat   = mom6util.read_mom6grid(dfgrid_mom, grid=grid_shape, grdpnt=grid_var)
HHM          = mom6util.read_mom6depth(dftopo_mom)
jdm          = np.shape(HHM)[0]
idm          = np.shape(HHM)[1]
# Convert to 0:360:
hlon = (hlon + 360) % 360
#hlon_180 = (hlon + 180) % 360 - 180


pthrtofs_input = '/gpfs/f6/drsa-hurr1/world-shared/save/Maria.Aristizabal/Scripts_to_prep_MOM6/Files_to_create_MOM6_OBC_for_Dmitry'

import mod_colormaps as mclrmp
clrmp_name = 'winter'
clr_ramp   = [1, 1, 1]   # add white color at the end of the colormap
clrmp = mclrmp.addendclr_colormap(clrmp_name, clr_ramp, nramp=0.1, ramp_start=False)
clrmp.set_bad(color=[0.0, 0., 0.])
rmin = -7000.
rmax = 0.

from mpl_toolkits.basemap import Basemap, cm

lon0 = 220.
lat0 = 50.
res  = 'l'
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
xR, yR = m(hlon,hlat)
PMsk = ( (xR > 1e20) | (yR > 1e20) )
AA = HHG.copy()
AA = np.insert(AA, 0, AA[:,-1], axis=1)
AA = np.insert(AA, -1, AA[-1,:], axis=0)
AA[PMsk] = np.nan
xR[PMsk]   = 1.e30
yR[PMsk]   = 1.e30

ny, nx = HHG.shape
AA = AA[0:ny, 0:nx]
AA = np.where(HHG >= 0., np.nan, AA)


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

f_ortho = False
if f_ortho:
  m.drawcoastlines()
  im1 = m.pcolormesh(xR, yR, AA, cmap=clrmp, vmin=rmin, vmax=rmax)

  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

HHM = np.where(HHM>=0, np.nan, HHM)
im1 = ax1.pcolormesh(HHM, cmap=clrmp, vmin=rmin, vmax=rmax)

sttl = f"NN domain"
ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='min')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.0f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_domain.py'
bottom_text(btx, fsz=8)





