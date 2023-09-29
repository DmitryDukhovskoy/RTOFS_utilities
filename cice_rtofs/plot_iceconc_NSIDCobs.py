"""
  Plot sea ice concentration
  Observed data 
  Sea Ice Concentration NSIDC Near-Real-Time Climate Data Record V2, Arctic
  or Antarctic
  Polar Portal NSIDC
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

# initial date of f/cast
#rdate = '20230307'
fdate = '20220909'    # date to plot
fnmb  = mtime.rdate2datenum(fdate)
DVf   = mtime.datevec(fnmb)


pthscr  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/'
pthgrid = pthscr+'hycom_fix/'
ftopo   = 'regional.depth'
fgrid   = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

print("Plotting NSIDC CDRv2: {0}/{1}/{2}".format(DVf[0],DVf[1],DVf[2]))


# Select Arctic region for RTOFS
import mod_cice_utils as mcice
#importlib.reload(mcice)
LONA = mcice.sub_region2D(LON, region='Arctic')
LATA = mcice.sub_region2D(LAT, region='Arctic')
HHA  = mcice.sub_region2D(HH, region='Arctic')

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
ny = LONA.shape[0]; nx = LONA.shape[1]
xP, yP = m(0.,90.)
# Read observations
pthi  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/OBS/iceconc_nsidc_cdr/'
fice  = 'nsidc_iceCDR_09012022_09092022_NH.nc'
fpice = pthi + fice 
  
Xobs, Yobs, Cice = mcice.read_iceconc_cdr(fpice, fnmb, xP, yP, lon0=0.)


plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.2, 0.75, 0.75])

# draw parallels.
parallels = np.arange(0.,88.,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,359.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)


im1 = ax1.pcolormesh(Xobs, Yobs, Cice, shading='flat',
                     cmap=cmpice,
                     vmin=rmin, vmax=rmax)

stl = 'NOAA NSIDCS CDRv2 IceConc {0}/{1:02d}/{2:02d}'.\
       format(DVf[0],DVf[1],DVf[2]) 
ax1.set_title(stl)

ax2 = fig1.add_axes([ax1.get_position().x0, ax1.get_position().y0-0.1, 
                     ax1.get_position().width,0.02])
clb = plt.colorbar(im1, cax=ax2, orientation='horizontal')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_iceconc_NSIDCobs.py'
mufig.bottom_text(btx, pos=[0.1, 0.03])


  



