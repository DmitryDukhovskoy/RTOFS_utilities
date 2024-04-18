"""
  Derive mean SSH 
  GOFS reanalysis
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
import matplotlib.animation as animation
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
import mod_mom6 as mom6util
import mod_gofs31 as mgofs

# Region domain
lat1 = 22.
lat2 = 48.
lon1 = 235.42
lon2 = 251.

YR   = 1994
HR   = 12

# Western coast for GLBv grid:
is1 = 0
is2 = 1000
js1 = 1775
js2 = 2800


pthoutp  = '/work/Dmitry.Dukhovskoy/data/gofs31_btmu_westcoast/'
pthpkl  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_tmp/'

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

varnm = 'surf_el'

# Load mean field:
floutp  = f"gofs31_53X_sshmean_westcoast_{YR}.pkl"
dfloutp = os.path.join(pthoutp,floutp)

print(f"Loading mean SSH --> {dfloutp}")
with open(dfloutp,'rb') as fid:
  [Amn, HH, LON, LAT] = pickle.load(fid)

cmpr = mutil.colormap_ssh(nclrs=200)
rmin = -0.25
rmax = 0.25
cmpr.set_bad(color=[0.1, 0.1, 0.1])


# Set up orthographic projection
from mpl_toolkits.basemap import Basemap, cm

# Add extra row/col for plotting with pcolormesh
lonw = LON.copy()
latw = LAT.copy()
lonw = np.insert(lonw, -1, lonw[:,-1]+0.01, axis=1)
lonw = np.insert(lonw, -1, lonw[-1,:]+0.01, axis=0)
latw = np.insert(latw, -1, latw[-1,:]+0.01, axis=0)
latw = np.insert(latw, -1, latw[:,-1]+0.01, axis=1)

lon0 = 220.
lat0 = 50.
res  = 'l'
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
xR, yR = m(lonw, latw)
PMsk = ( (xR > 1e20) | (yR > 1e20) )

plt.ion()
fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)

istart = -1e30
ndays = 365
iloop = 0
def animate(irc):
  global iloop
  if irc > ndays:
    print('Step {0} Stopping ...'.format(irc))
    anim.event_source.stop()
    iloop += 1
    return

  if iloop > 0:
# To stop FuncAnimation from looping after 1 time
    return

  jday   = irc+1
  dnmb   = mtime.jday2dnmb(YR, jday)
  dv     = mtime.datevec(dnmb)
  ds     = mtime.datestr(dnmb)
  MM     = dv[1]
  DM     = dv[2]

  nexpt, _  = mgofs.gofs31_expt53X(dnmb)
  expt      = int(nexpt*10)
  print(f'Plotting {irc} {YR}/{MM}/{DM} expt:{expt}')

  urlB  = 'https://tds.hycom.org/thredds/dodsC/datasets/GLBv0.08/'
  urlF  = f'{urlB}expt_53.X/data/{YR}'
  urlT  = f'{urlB}expt_53.X/topo/'
  ftopo = 'depth_GLBv0.08_11.nc'
  fnmnc = f'hycom_GLBv0.08_{expt}_{YR}{MM:02d}{DM:02d}{HR:02d}_t000.nc'
  furl  = os.path.join(urlF, fnmnc)
  fturl = os.path.join(urlT, ftopo)

# Concatenate data over the 180W line - where this grid begins
#  lon1d = read_field(furl, 'lon')
#  lon_topo = read_field(fturl, 'Longitude')
#  lon_topo = np.where(lon_topo>=180., lon_topo-360., lon_topo)
#  istart   = np.argwhere(lon_topo == lon1d[0])[0][0]

  try:
    print(f'URL: {furl}')
    ssh  = read_field(furl, varnm)
  except:
    print(f"Could not read {furl}, irc={irc} skipping ...")
    return

  ssh = np.where(ssh < -100., np.nan, ssh)
  A2d = ssh[js1:js2+1, is1:is2+1]
  A2d = A2d - Amn

  ctitle = f'GOFS3.1-{nexpt:3.1f} dltSSH {YR}/{MM:02d}/{DM:02d}'

  data = A2d.copy()
  data = np.insert(data, -1, data[:,-1], axis=1)
  data = np.insert(data, -1, data[-1,:], axis=0)
  data[PMsk] = np.nan
  xR[PMsk]   = 1.e30
  yR[PMsk]   = 1.e30

  nx   = A2d.shape[1]
  ny   = A2d.shape[0]
  data = data[0:ny, 0:nx]


  plt.clf()

  ax1 = plt.axes([0.04, 0.04, 0.8, 0.8])
  m.drawcoastlines()
  im1 = m.pcolormesh(xR, yR, data, shading='flat', cmap=cmpr,\
                     vmin=rmin, vmax=rmax)

  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  ax1.set_title(ctitle)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='both')
  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=10)

  #fig1.colorbar(im,ax=ax1,orientation='horizontal')
  btx = 'anim_dltssh_gofs.py'
  bottom_text(btx, pos=[0.02, 0.02])

  return im1

# Keep frames = # of recrods + 1 to trigger the stop option to avoid
# repeated loops, repeat=False does not work
#anim = animation.FuncAnimation(fig1, animate, frames=10,  repeat=False)
anim = animation.FuncAnimation(fig1, animate, frames=ndays+1,  repeat=False)
fig1.canvas.draw()
anim.event_source.stop() 
#plt.draw()
#plt.show()

# Define the meta data for the movie:
#FFMpegWrite = animation.writers['ffmpeg']
metadata = dict(title=f'GOFS3.1-53X dltSSH', artist='Matplotlib',
                comment='created: anim_dltssh_gofs.py')
FFwriter = animation.FFMpegWriter(fps=5, metadata=metadata)
flout = f'gofs31_53X_dltSSH_{YR}.mp4'
drflout = os.path.join(pthpkl, flout)
print('Saving ' + drflout)
anim.save(drflout, writer = FFwriter)

anim.event_source.stop()


