"""
  Plot sea ice concentration
  from CICE output netcdf files

  Southern Hemisphere

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

# initial date of f/cast
#rdate = '20230307'
rdate = '20220901'  
hr    = 192   # f/cast hour
f_obs = True  # plot observed ice edge, NSIDC CDRv2 25 km 


YR, MM, DM, HR = mtime.parse_rdate(rdate)

# Find forecast date:
rnmb  = mtime.rdate2datenum(rdate)
ndays = np.floor(hr/24)
fnmb  = rnmb + hr/24.
fdate = mtime.dnumb2rdate(fnmb)
DVf   = mtime.datevec(fnmb)


pthscr  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/'
pthgrid = pthscr+'hycom_fix/'
pthin   = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/wcoss2.prod/rtofs.'\
          + rdate + '/'
flnm    = 'rtofs_glo_2ds_f{0:03d}_ice.nc'.format(hr)
fpinp   = pthin+flnm

fld   = 'ice_coverage'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

print("Plotting {3}: {0}/{1}/{2}".format(YR,MM,DM,fld))

nc = NetCDFFile(fpinp)
#aice = nc.variables["aice"][:].squeeze().data
A2D = nc.variables[fld][:].squeeze().data
A2D = np.where(A2D > 1.e20, np.nan, A2D)


print(fld + " shape=",A2D.shape)

# Select Antarctic region:
import mod_cice_utils as mcice
#importlib.reload(mcice)
LONA = mcice.sub_region2D(LON, region='Antarctic')
LATA = mcice.sub_region2D(LAT, region='Antarctic')
HHA  = mcice.sub_region2D(HH, region='Antarctic')
A2D = mcice.sub_region2D(A2D, region='Antarctic')


#
# Create colormaps for sea ice conc
importlib.reload(mcice)
cmpice = mcice.colormap_conc()
#cmpice = mclrs.create_colormap(CLR,Ncmp)
cmpice.set_bad(color=[0.3, 0.3, 0.3])
cmpice.set_under(color=[1, 1, 1])
rmin = 0.15
rmax = 1.

m = Basemap(projection='spstere',boundinglat=-50,lon_0=0,resolution='l')
ny = A2D.shape[0]; nx = A2D.shape[1]
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly spaced grid.
x, y = m(lons, lats) # compute map proj coordinates.
xh, yh = m(LONA,LATA)  # hycom coordinates on the projections

mh, nh = xh.shape[0], xh.shape[1]
xP, yP = m(0.,-90.)
# Read observations
if f_obs:
  pthi  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/OBS/iceconc_nsidc_cdr/'
  fice  = 'nsidc_iceCDR_09012022_09152022_SH.nc'
  fpice = pthi + fice 
  
  Xobs, Yobs, Cice = mcice.read_iceconc_cdr(fpice, fnmb, 
                     xP, yP, regn='Antarctic',lon0=0.)


plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.2, 0.75, 0.75])

# draw parallels.
parallels = np.arange(-90,-10,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,359.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)


im1 = ax1.pcolormesh(xh, yh, A2D, shading='flat',
                     cmap=cmpice,
                     vmin=rmin, vmax=rmax)

# Plot contour of observed ice
if f_obs:
  ax1.contour(Xobs,Yobs,Cice,[0.15], linestyles='solid', 
              colors=[(0,0,0)], linewidths=1)
                  

stl = 'RTOFSv2 CICE IceConc InitDate: {1}/{2:02d}/{3:02d}, F/cast: {4} hrs'.\
       format(fld,YR,MM,DM,hr) 
ax1.set_title(stl)

ax2 = fig1.add_axes([ax1.get_position().x0, ax1.get_position().y0-0.1, 
                     ax1.get_position().width,0.02])
clb = plt.colorbar(im1, cax=ax2, orientation='horizontal')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.tick_params(direction='in', length=12)

# Info
ssinf = 'RTOFS run {0}/{1:02d}/{2:02d}\n'.format(YR,MM,DM)
ssinf = ssinf + 'F/cast date: {0}/{1:02d}/{2:02d}:{3:02d}\n'.\
           format(DVf[0],DVf[1],DVf[2],DVf[3])
ssinf = ssinf + flnm
ax3 = plt.axes([0.105, 0.205, 0.23, 0.08])
xE  = 0.2
yE  = 0.58
if f_obs:
  ax3.plot([0.05,xE],[yE,yE],'-',color=[0,0,0])
  ax3.text(xE+0.2*xE, yE, 'NOAA/NSIDC CDRv2',verticalalignment='center')
ax3.text(0.03, 0.04, ssinf)
ax3.set_xlim([0, 1.2])
ax3.set_ylim([0, 0.8])
ax3.set_xticks([])
ax3.set_yticks([])

btx = 'plot_iceconc_productSH.py'
mufig.bottom_text(btx, pos=[0.1, 0.03])



# Plot observed ice concentration:
#  im2 = ax1.pcolormesh(Xobs, Yobs, Cice, shading='flat',
#                     cmap=cmpice,
#                     vmin=rmin, vmax=rmax)

  



