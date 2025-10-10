# Plot HYCOM SSH from GOFS3.1 
# used as IC for MOM6
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
#import struct
import datetime
#import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
from netCDF4 import Dataset as ncFile

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
import mod_cice6_utils as mc6util
#import mod_valid_utils as mvutil


expt    = '93.0'
YR      = 2020
jday    = 1
HR      = 12
rg      = 9806.
hg      = 1.e15
#regn_nm = 'GOM' # GOM, Carib, GulfStr, SOcean, Kurosh, Agulhas
#regn_nm = 'GulfStr'
#regn_nm = 'NPac1'
#regn_nm = 'GINSea'
regn_nm = 'Agulhas'
#regn_nm = 'Arctic'  # no subset needed 
#regn_nm = 'Antrct'  # no subset needed 

dnmb = mtime.jday2dnmb(YR,jday)
DV   = mtime.datevec(dnmb)
MM   = DV[1]
DD   = DV[2]

pthhcm = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/GLBb0.08_expt93.0/'
flhcm  = '930_archv.2020_001_00'
fina   = pthhcm+flhcm+'.a'
finb   = pthhcm+flhcm+'.b'

# Read ssh:
ssh,nn,mm,ll  = mhycom.read_hycom(fina,finb,'srfhgt')
ssh[ssh>hg] = np.nan
ssh = ssh/9.806


idm, jdm, kdm = mhycom.hycom_dim(fina,finb)
ijdm = idm*kdm


pthgrid = pthhcm
ftopo   = 'depth_GLBb0.08_09m11'
fgrid   = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

# Subset region
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
X, Y = np.meshgrid(np.arange(idm), np.arange(jdm))
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
stl  = 'SSH, {5}, GOFS3.1-{3}, {0}/{1}/{2}, dssh={4:.3f}'.\
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


btx = 'plot_ssh_gofs.py'
bottom_text(btx, pos=[0.08, 0.02])

f_setrmu = False
if f_setrmu:
# Bind the button_press_event with the onclick() method
  fig1.canvas.mpl_connect('button_press_event', onclick)


