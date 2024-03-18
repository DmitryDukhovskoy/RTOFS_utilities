# 1st baroclinic Rossby radius 
# 2 approaches: WKB-theory and numerically solve eig/value problem
#
# Use bottom topography for depths deeper than WOA last depth level (-5500m)
#  Topography interpolated from ETOPO2
#
#
# Calculate N2 following Chelton 1996
# N2 is calculated at the middle of the depth interval
# density is adiabatically adjusted to the mid-grid depth
#
# !!!!!!!!!!!!!!!
#
# Problems: (1) max depth in climat -5500 m - need to bottom
# (2) a few cases when min(eig/value)>0 --- check out
#
# !!!!!!!!!!!!!!!!
#
##reset -f 
from IPython import get_ipython
get_ipython().magic('reset -sf')

import os
import numpy as np
import matplotlib.pyplot as plt
#import torch
import sys
import pdb
import importlib
import matplotlib.mlab as mlab
from matplotlib.colors import ListedColormap
#from mpl_toolkits.basemap import Basemap, shiftgrid
#from scipy.interpolate import griddata

sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
import mod_utils_fig
# importlib.reload(mod_utils_fig)
from mod_utils_fig import bottom_text
from mod_utils_fig import plot_fld2D

idm = 1440
jdm = 720

grd=0.25
if grd==0.25:
  cgrd=4
woa='woa18'
seas=13    # season: 1-12 - monthly, 13-16 - spring, summer, ...

pthout = '/data/ncei2/work/dmitry.dukhovskoy/data_output/'

btx = 'plot_rrossby_topo.py'

#


# ------------
# Read TOPO
# -------------
pthout = '/data/ncei2/work/dmitry.dukhovskoy/topo/woa/'
fout = pthout+'etopo2woa025_bilinear.dat'
ijdm = idm*jdm
fid = open(fout,'rb')
try:
  fid.seek(0)
  HW = np.fromfile(fid,dtype='<f4', count=ijdm)
finally:
  fid.close()
HW = HW.reshape(jdm,idm)


# Loading saved Rossby R.
fout1 = pthout+'Rrossby_num_v2.dat'
fid = open(fout1,'rb')
data = np.fromfile(fid,dtype='f', count=idm*jdm)
fid.close()
RsbNum = data.reshape((jdm,idm))

ctitle='1st Barocl Rossby Radius, WOA18'
cl1=5.
cl2=245.

plt.ion()
cmpS = plt.cm.get_cmap('turbo')
fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
im = ax1.pcolormesh(RsbNum,cmap=cmpS,shading='flat')
#ax1.axis('scaled')
ax1.set(xlim=(0,idm),ylim=(550,jdm))
#ax1.autoscale(enable=True, axis='both', tight=True)

ax1.set_title(ctitle)
fig1.colorbar(im,ax=ax1,orientation='horizontal')

im.set_clim(cl1,cl2)


plt.show()
figout='plot_field.jpg'
plt.savefig(figout)
plt.close(fig1)

wwwpth='/net/www-dev/www/ncei_wod/'
wwwfig='pthn_fig.jpg'
cmd = ('mv {0} {1}{2}'.format(figout,wwwpth,wwwfig)) 
os.system(cmd)  


from mod_plot_anls import zonalavrg
ctl2='1st Barocl Rossby R zonal avrg, {0}, {1}'.format(tfnm,sfnm)
zonalavrg(RsbNum,ctl2,lat,LMsk,btx=btx,ifg=1)

# Test: plot eigenfunctions
f_plt=0
if f_plt>0:
  plt.ion()
  fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
  plt.clf()

  im=6

  Rrsb = RsbNum[jj,ii]
  x0 = lon[ii]
  y0 = lat[jj]
  ww = W[im]
  vvk = V[:,im]
  zzk = Z_phi[0:kbtm+1]  
  nzk = zzk.shape[0]
# 
# Add surface and bottom to eig/functions
  zzV=np.zeros(nzk+2)
  zzV[0] = 0.  
  zzV[1:nzk+1]=zzk
  zzV[nzk+1]=zbtm

# Add 0 at the ends for eig/functions:
  vvV = np.zeros(zzV.shape[0])
  vvV[1:nzk+1] = vvk

  plt.plot(vvV,zzV,'.-')
  ctl = 'Eig/vector ii={2}, jj={3}, {4:5.2f}E, {5:5.2f}N, im={0}, Rr={1:6.0f} km'.\
         format(im,Rrsb,ii,jj,x0,y0)
  plt.title(ctl)





