# Create animation of SSH from MOM6 run
# Not polar regions
#  
# matplotlib.org/stable/gallery/animation/simple_anim.html
#
# To save an mp4 animation:
# Make sure that python environment has ffmpeg installed:
#(https://stackoverflow.com/questions/40509002/
#  ffmpeg-is-not-being-detected-by-spyder): 
# "You need to install a copy of ffmpeg that can be recognized by Anaconda. 
# Please run this command in a terminal to do that":
# e.g. for me:
# eval "$($PYPATH/bin/conda shell.bash hook)"
# conda activate anls
# conda install -c conda-forge ffmpeg
# 
# Transfer: 
# (1) hera ---> Niagara trusted 
# scp 008mom6cice6_003_ssh.mp4 Dmitry.Dukhovskoy@dtn-niagara.fairmont.rdhpcs.noaa.gov:/collab1/data/Dmitry.Dukhovskoy/.
# (2) Niagara trusted ---> Niagara untrusted
# (3) from local machine:
# DSD@MBpro:~\> rsync -av -e ssh --progress Dmitry.Dukhovskoy@udtn-niagara.fairmont.rdhpcs.noaa.gov:/collab1/data_untrusted/Dmitry.Dukhovskoy/008mom6cice6_003_ssh.mp4 .
#
#


import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pdb
import importlib
#import struct
import datetime
#import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import matplotlib.animation as animation
import time
from netCDF4 import Dataset as ncFile

#sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_cice6_utils as mc6util
#import mod_valid_utils as mvutil

expt  = '003'
YR1   = 2020
mo1   = 1
mday1 = 1
YR2   = 2020
mo2   = 12
mday2 = 31
dday  = 5

HR      = 12
hg      = 1.e15
regn_nm = 'GOM' # GOM, Carib, GulfStr, SOcean, Kurosh, Agulhas
#regn_nm = 'GulfStr'
#regn_nm = 'NPac1'
#regn_nm = 'GINSea'
#regn_nm = 'Agulhas'
#regn_nm = 'Arctic'  # no subset needed 
#regn_nm = 'Antrct'  # no subset needed 

pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'
pthout = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/figs/MOM6_CICE6/expt' + \
          expt + '/animations/' 

# Create time array for animation:
import mod_mom6 as mom6util
importlib.reload(mom6util)
date1 = [YR1, mo1, mday1]
date2 = [YR2, mo2, mday2]
TPLT  = mom6util.create_time_array(date1, date2, dday)

pthgrid  = pthrun + 'INPUT/'
fgrd_mom = pthgrid + 'regional.mom6.nc' 
ftopo_mom= pthgrid + 'ocean_topog.nc'
LON, LAT = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
HH       = mom6util.read_mom6depth(ftopo_mom)
Lmsk     = mom6util.read_mom6lmask(ftopo_mom)
idma     = LON.shape[1]
jdma     = LON.shape[0]
ijdma    = idma*jdma


#importlib.reload(mvutil)
importlib.reload(mutil)
REGNS  = mutil.rtofs_reg2Dmaps()
REGNMS = list(REGNS.keys())
xl1    = REGNS[regn_nm]["xl1"]
xl2    = REGNS[regn_nm]["xl2"]
yl1    = REGNS[regn_nm]["yl1"]
yl2    = REGNS[regn_nm]["yl2"]
Ip     = REGNS[regn_nm]["Ip"]
Jp     = REGNS[regn_nm]["Jp"]

plt.ion()
# Points inside the region
import mod_misc1 as mmisc
#importlib.reload(mmisc)
X, Y = np.meshgrid(np.arange(idma), np.arange(jdma))
Rmsk, IRg, JRg = mmisc.inpolygon_v2(X,Y,Ip,Jp)  # region


def is_empty(self):
  """
  Return True if there is no visible artist in the figure
  """

  children = self.get_children()
  if len(children) < 1:
      return True
  
  for child in self.get_children():
      if child.get_visible():
          return False

  return True

from matplotlib import cm
from copy import copy

dssh       = 0.2
ssh_cntrs  = np.arange(0,1.5,dssh)
ssh_ncntrs = np.arange(-1.2,-0.01,dssh)
clr_cntrs  = [(0,0,0)]
clr_ncntrs = [(0,0.4,0.5)]
clrmp = copy(plt.cm.coolwarm)
clrmp.set_bad(color=[0.3,0.3,0.3])
rmin = -0.5
rmax = 0.5


def init():
  im1 = ax1.pcolormesh(Lmsk, vmin=rmin, vmax=rmax, cmap=clrmp)
  cr1 = ax1.contour(Lmsk, [0.9], linestyles='solid', 
                colors=[(0,0,0)], linewidths=1)
  cr2 = ax1.contour(Lmsk, [0.1], linestyles='solid', 
                colors=[(0,0,0)], linewidths=1)

  return im1, cr1, cr2

# Not polar projections:
# Plot region and select points if needed:
fig1 = plt.figure(1,figsize=(9,8), constrained_layout=False)
fig1.clf()
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
cax = fig1.add_axes([0.92, 0.3, 0.015, 0.4])
im1 = ax1.pcolormesh(Lmsk, vmin=rmin, vmax=rmax, cmap=clrmp)
cr1 = ax1.contour(Lmsk, [0.9], linestyles='solid',
              colors=[(0,0,0)], linewidths=1)
cr2 = ax1.contour(Lmsk, [0.1], linestyles='solid',
              colors=[(0,0,0)], linewidths=1)
clb = fig1.colorbar(im1, cax=cax, extend='both')
cax.set_yticklabels(cax.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=5)

btx = 'anim_ssh.py'
bottom_text(btx, pos=[0.08, 0.02])
plt.draw()
plt.show()

icc = 0
while is_empty(fig1):
  print('Waiting for figure to be created {0} ...'.format(icc))
  time.sleep(5)
  icc += 1
  if icc > 100:
    raise Exception('Error: Waiting too long')



#ims  = []
nrec = TPLT.shape[0]
#nrec  = 3
iloop = 0
def animate(irc):
  global iloop
#  global anim
  if irc == nrec:
    print('Step {0} Stopping ...'.format(irc))
    anim.event_source.stop()
    iloop += 1
    return

  if iloop > 0:
# To stop FuncAnimation from looping after 1 time
    return

  dnmb = TPLT['dnmb'][irc]
  DV   = mtime.datevec(dnmb)
  YR   = DV[0]
  MM   = DV[1]
  DD   = DV[2]
  jday = int(TPLT['yrday'][irc])
  print('Plotting {0} {1}/{2}/{3}'.format(irc,YR,MM,DD))

  pthbin = pthrun + 'tarmom_{0}{1:02d}/'.format(YR,MM)
  flin   = pthbin + 'ocnm_{0}_{1:03d}_{2}.nc'.format(YR,jday,HR)
  nc     = ncFile(flin,'r')

  ssh   = nc.variables['SSH'][:].squeeze().data
  ssh   = np.where(ssh > hg, np.nan, ssh)

  # Demean ssh using regional mean:
  ssh_sub = ssh[JRg,IRg]
  ssh_mn  = np.nanmean(ssh_sub)
  dSSH    = ssh - ssh_mn

  stl  = 'SSH, {5}, 0.08MOM6-CICE6-{3}, {0}/{1}/{2}, dssh={4:.3f}'.\
       format(DV[0],DV[1],DV[2],expt,dssh, regn_nm) 

  ax1.cla()
  im1  = ax1.pcolormesh(dSSH, vmin=rmin, vmax=rmax, cmap=clrmp)
  cr1  = ax1.contour(dSSH, ssh_cntrs, linestyles='solid', 
              colors=[(0,0,0)], linewidths=1)
  cr2  = ax1.contour(dSSH, ssh_ncntrs, linestyles='solid', 
              colors=[(0,0.4,0.5)], linewidths=1)

  ax1.clabel(cr1, cr1.levels, inline=False, fmt='%3.1f', fontsize=10)
  ax1.clabel(cr2, cr2.levels, inline=False, fmt='%3.1f', fontsize=10)
  ax1.axis('scaled')
  ax1.set_xlim([xl1,xl2])
  ax1.set_ylim([yl1,yl2])
  ax1.set_xticks([])
  ax1.set_yticks([])

  ax1.set_title(stl)


#  print('ok')

  return im1, cr1, cr2,

# Keep frames = # of recrods + 1 to trigger the stop option to avoid
# repeated loops, repeat=False does not work
#anim = animation.FuncAnimation(fig1, animate, frames=len(TPLT),  repeat=False)
anim = animation.FuncAnimation(fig1, animate, frames=nrec+1,  repeat=False)
# anim.event_source.stop() 
plt.draw()
#plt.show()

# Define the meta data for the movie:
#FFMpegWrite = animation.writers['ffmpeg']
metadata = dict(title='SSH, MOM6-CICE6 ' + expt, artist='Matplotlib',
                comment='created: anim_ssh.py')
FFwriter = animation.FFMpegWriter(fps=5, metadata=metadata) 
flout = '008mom6cice6_{0}_ssh.mp4'.format(expt)
drflout = pthout + flout
print('Saving ' + drflout)
anim.save(drflout, writer = FFwriter)

anim.event_source.stop()



