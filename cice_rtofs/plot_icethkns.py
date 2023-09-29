"""
  Plot sea ice thickness
  from CICE output netcdf files

  For production forecasts, use script to download CICE output
  /home/Dmitry.Dukhovskoy/scripts/rtofs/get_production_cice.sh

"""
import os
import numpy as np
import sys
import importlib
import matplotlib.pyplot as plt
import datetime
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset as NetCDFFile
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.colors import ListedColormap
import importlib
#import mod_colormaps


sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

import mod_utils_fig as mufig
import mod_time as mtime
import mod_read_hycom as mhycom
import mod_colormaps as mclrs

#rdate = '20230307'
rdate = '20220905'
YR, MM, DM, HR = mtime.parse_rdate(rdate)

pthscr  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/'
pthgrid = pthscr+'hycom_fix/'
pthin   = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_paraCd/rtofs.'\
          + rdate + '/'
flnm    = 'cice_inst.{0}-{1:02d}-{2:02d}-00000.nc'.format(YR,MM,DM)
fpinp   = pthin+flnm

fld   = 'aice'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

print("Plotting {3}: {0}/{1}/{2}".format(YR,MM,DM,fld))

nc = NetCDFFile(fpinp)
#aice = nc.variables["aice"][:].squeeze().data
A2D = nc.variables[fld][:].squeeze().data
A2D = np.where(A2D > 1.e20, np.nan, A2D)

print(fld + " shape=",A2D.shape)

# Select Arctic region:
import mod_cice_utils as mcice
#importlib.reload(mcice)
LONA = mcice.sub_region2D(LON, region='Arctic')
LATA = mcice.sub_region2D(LAT, region='Arctic')
HHA  = mcice.sub_region2D(HH, region='Arctic')
A2D = mcice.sub_region2D(A2D, region='Arctic')

#
# Create colormaps for sea ice conc
importlib.reload(mcice)
cmpice = mcice.colormap_conc()
#cmpice = mclrs.create_colormap(CLR,Ncmp)
cmpice.set_bad(color=[0.3, 0.3, 0.3])
cmpice.set_under(color=[1, 1, 1])
rmin = 0.15
rmax = 1.

m = Basemap(projection='npstere',boundinglat=50,lon_0=0,resolution='l')


plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.2, 0.75, 0.75])

# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,360.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)


ny = A2D.shape[0]; nx = A2D.shape[1]

lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly spaced grid.
x, y = m(lons, lats) # compute map proj coordinates.
xh, yh = m(LONA,LATA)  # hycom coordinates on the projections
im1 = ax1.pcolormesh(xh, yh, A2D, shading='flat',
                     cmap=cmpice,
                     vmin=rmin, vmax=rmax)
                  

stl = 'RTOFSv2 CICE {0} {1}/{2:02d}/{3:02d}'.format(fld,YR,MM,DM) 
ax1.set_title(stl)

ax2 = fig1.add_axes([ax1.get_position().x0, ax1.get_position().y0-0.1, 
                     ax1.get_position().width,0.02])
clb = plt.colorbar(im1, cax=ax2, orientation='horizontal')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.tick_params(direction='in', length=12)



btx = 'plot_iceconc.py'
mufig.bottom_text(btx, pos=[0.1, 0.03])










