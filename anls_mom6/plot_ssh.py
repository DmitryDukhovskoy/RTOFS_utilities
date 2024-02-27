# Plot SSH from MOM6 run
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
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


expt    = '003'
YR      = 2020
jday    = 45
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

dnmb = mtime.jday2dnmb(YR,jday)
DV   = mtime.datevec(dnmb)
MM   = DV[1]
DD   = DV[2]

pthbin = pthrun + 'tarmom_{0}{1:02d}/'.format(YR,MM)
flin   = pthbin + 'ocnm_{0}_{1:03d}_{2}.nc'.format(YR,jday,HR)
nc     = ncFile(flin,'r')

ssh   = nc.variables['SSH'][:].squeeze().data
LON   = nc.variables['xh'][:].squeeze().data
LAT   = nc.variables['yh'][:].squeeze().data
idma  = LON.shape[0]
jdma  = LAT.shape[0]
ijdma = idma*jdma

import mod_mom6 as mom6util
importlib.reload(mom6util)
pthgrid  = pthrun + 'INPUT/'
fgrd_mom = pthgrid + 'regional.mom6.nc' 
ftopo_mom= pthgrid + 'ocean_topog.nc'
LON, LAT = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
HH       = mom6util.read_mom6depth(ftopo_mom)
Lmsk     = mom6util.read_mom6lmask(ftopo_mom)

# Get lon/lat from hycom
# Note HYCOM has 1 extra row 
#import mod_read_hycom as mhycom
#pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/GLBb0.08_expt93.0/'
#ftopo   = 'depth_GLBb0.08_09m11'
#fgrid   = 'regional.grid'
#LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
#Lmsk    = np.where(HH>0, 0, 1)

ssh = np.where(ssh > hg, np.nan, ssh)

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

# Points inside the region
import mod_misc1 as mmisc
#importlib.reload(mmisc)
X, Y = np.meshgrid(np.arange(idma), np.arange(jdma))
Rmsk, IRg, JRg = mmisc.inpolygon_v2(X,Y,Ip,Jp)  # region

# Demean ssh using regional mean:
ssh_sub = ssh[JRg,IRg]
ssh_mn  = np.nanmean(ssh_sub)
dSSH    = ssh - ssh_mn



# Function to print mouse click event coordinates
def onclick(event):
   print([event.xdata, event.ydata])


plt.ion()

from matplotlib import cm
from copy import copy

dssh       = 0.2
ssh_cntrs  = np.arange(0,1.5,dssh)
ssh_ncntrs = np.arange(-1.2,-0.01,dssh)
clr_cntrs  = [(0,0,0)]
clr_ncntrs = [(0,0.4,0.5)]
stl = 'SSH, {5}, 0.08MOM6-CICE6-{3}, {0}/{1}/{2}, dssh={4:.3f}'.\
       format(DV[0],DV[1],DV[2],expt,dssh, regn_nm) 
clrmp = copy(plt.cm.coolwarm)
clrmp.set_bad(color=[0.3,0.3,0.3])
rmin = -0.5
rmax = 0.5


if regn_nm == 'Arctic' or regn_nm == 'Antrct':
  Nclrs = 200
  ixm=np.linspace(0, 1, num=Nclrs, endpoint=True) 
  CLR = clrmp(ixm)
  CLR = CLR[:,0:3]

  import mod_colormaps as mclrs
  CMP = mclrs.create_colormap(CLR, Nclrs)

  mc6util.plot_polar_2D(LON, LAT, dSSH, region=regn_nm, nfg=1, \
           rmin=rmin, rmax=rmax, cmpice=CMP, stl=stl, \
           cntr1=ssh_cntrs, clr1=clr_cntrs, \
           cntr2=ssh_ncntrs, clr2=clr_ncntrs)

else:
# Not polar projections:
# Plot region and select points if needed:
  fig1 = plt.figure(1,figsize=(9,8), constrained_layout=False)
  fig1.clf()


  ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8],)
  im1 = ax1.pcolormesh(dSSH, vmin=rmin, vmax=rmax, cmap=clrmp)
  ax1.contour(dSSH, ssh_cntrs, linestyles='solid', 
              colors=[(0,0,0)], linewidths=1)
  ax1.contour(dSSH, ssh_ncntrs, linestyles='solid', 
              colors=[(0,0.4,0.5)], linewidths=1)
  ax1.axis('scaled')
  ax1.set_xlim([xl1,xl2])
  ax1.set_ylim([yl1,yl2])

  ax1.set_title(stl)

  cax = fig1.add_axes([0.92, 0.3, 0.015, 0.4])
  clb = fig1.colorbar(im1, cax=cax, extend='both')
  cax.set_yticklabels(cax.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=5)

btx = 'plot_ssh.py'
bottom_text(btx, pos=[0.08, 0.02])

f_setrmu = False
if f_setrmu:
# Bind the button_press_event with the onclick() method
  fig1.canvas.mpl_connect('button_press_event', onclick)


