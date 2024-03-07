# Plot SSH from MOM6 run
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
from netCDF4 import Dataset as ncFile

PPTHN = '/home/Dmitry.Dukhovskoy/python'
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
import mod_cice6_utils as mc6util
#import mod_valid_utils as mvutil


expt    = '003'
YR      = 2021
jday    = 45
HR      = 12
hg      = 1.e15

pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'

dnmb = mtime.jday2dnmb(YR,jday)
DV   = mtime.datevec(dnmb)
MM   = DV[1]
DD   = DV[2]

pthbin = pthrun + 'tarmom_{0}{1:02d}/'.format(YR,MM)
flname = 'ocnm_{0}_{1:03d}_{2}.nc'.format(YR,jday,HR)
flin   = pthbin + flname
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

ssh = np.where(ssh > hg, np.nan, ssh)


from matplotlib import cm
#clrs   = cm.get_cmap('viridis',200)
#cmpr   = cm.get_cmap('PuOr',200)
cmpr = mutil.colormap_ssh(nclrs=100)
cmpr.set_bad(color=[0.2,0.2,0.2])
rmin = -1.2
rmax = 1.2

plt.ion()

btx  = 'plot_ssh_orthproj.py'
sttl = 'MOM6-CICE6 {0} SSH {1}/{2}/{3}, {4}'.format(expt, YR, MM, DD, flname)

importlib.reload(mutil)
mutil.plot_orthographic_proj(LON, LAT, ssh, cmpr=cmpr, sttl=sttl, \
        lon0=-50, lat0=40, rmin=rmin, rmax=rmax, btx=btx, fgnmb=1)


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


