"""
  Create HYCOM relaxation file for climatology fields
  2D array with relaxation weights s^-1 (e-folding scale)
  Sulu Sea
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
import mod_read_hycom as mrhycom
#import mod_write_hycom as mwhycom
import mod_utils_fig as mufig
import mod_misc1 as mmisc
importlib.reload(mmisc)

# Read topo
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo   = 'regional.depth'
fgrid   = 'regional.grid'

LON, LAT, HH = mrhycom.read_grid_topo(pthgrid,ftopo,fgrid)

IDM  = LON.shape[1]
JDM  = LON.shape[0]
IJDM = IDM*JDM
npad = 4096-IJDM%4096
toto = np.zeros((npad)).astype('float32')

# Relaxaion - stron in W. Caribbean
# deeper than zz0 
# Gaussian shape for relaxation 
# decreasing rls weights (s^-1) from
# the region of max rlx to the bndry
# West - bounded by land (isobath zz0)
zz0 = -250.
regn = 'SuluSea'
rlx_days = 3.    # days strongest relaxation
rlx_dayw = 30.   # days weekest rlx 
rlx_min = 1./(rlx_dayw*86400.)  # weekest rlx, sec^-1
rlx_max = 1./(rlx_days*86400.)  # strongest rlx, sec^-1

plt.ion()
# TO plot/select rmu regions:
f_plt = False
if f_plt:
  fgnmb = 1
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
  lvl = np.arange(-5000,0,500)
  im1 = ax1.contour(HH,levels=lvl)
  ax1.contour(HH,[0],colors=[(0,0,0)])

  ax1.axis('scaled')
  ax1.set_xlim([450,750])
  ax1.set_ylim([1500,1700])

  # Function to print mouse click event coordinates
  def onclick(event):
     print([event.xdata, event.ydata])

  # Set new rmu region
  # by selection coordinate of relax zones
  f_setrmu = True
  if f_setrmu:
  # Bind the button_press_event with the onclick() method
    fig1.canvas.mpl_connect('button_press_event', onclick)
  #
  # Can click on the map and read coordinates from the screen


# Region of max rlx:
# Sulu Sea:
II = np.array([539, 567, 586, 608, 613, 616, 611, 602, 584, 576, \
               569, 561, 544, 536])
JJ = np.array([1610, 1638, 1641, 1642, 1632, 1612, 1605, 1601, 1584, 1578, \
               1576, 1576, 1588, 1605])

# Region of non-zero relaxation
# Can be same as max relax or bigger to create some transitional zone 
# with weaker relax around max relax zone, if needed
IR = np.array([539, 567, 586, 608, 613, 616, 611, 602, 598, 586, \
               573, 546, 533, 535])
JR = np.array([1610, 1638, 1641, 1642, 1632, 1612, 1605, 1601, 1586, 1579, \
               1569, 1575, 1588, 1603])

# Find points inside the polygon
# similar to matlab inpolygon:
X, Y = np.meshgrid(np.arange(IDM), np.arange(JDM))
MSKx, IMx, JMx = mmisc.inpolygon_v2(X,Y,II,JJ)  # strong relaxation region
MSKr, IIr, JJr = mmisc.inpolygon_v2(X,Y,IR,JR)  # non-zero rlx region


# Derive contour coordinates of max rlx
cs = plt.contour(MSKx,[0.95])
for item in cs.collections:
   for i in item.get_paths():
      v  = i.vertices
      Xc = v[:, 0]
      Yc = v[:, 1]

# Find min distance
nI   = Xc.shape[0]
Dmin = np.zeros((JDM,IDM))+1.e20
IX0  = Dmin*0-1
IY0  = Dmin*0-1
print('Creating rmu ...')
for kk in range(nI):
  if kk%50 == 0:
    print(' processed {0:5.2f}% ...'.format(kk/nI*100.))
  x0, y0 = Xc[kk], Yc[kk]
  D = np.sqrt((X-x0)**2+(Y-y0)**2)
  Dmin = np.where(Dmin>D,D,Dmin)
  jp, ip = np.where(Dmin == D)
  IX0[jp,ip] = x0    # loc of the closes pnt on rlx bdnry
  IY0[jp,ip] = y0

cc   = 15.  # larger cc - wider transition zone from high-low relaxation
RMU  = rlx_max*np.exp(-((X-IX0)**2+(Y-IY0)**2)/(2.*cc**2))
RMU[JMx,IMx] = rlx_max
RMU[HH>zz0]  = rlx_min
RMU[RMU<rlx_min] = rlx_min
RMU[np.where(MSKr<1)] = 0.
RMU[np.where(HH>=0)] = np.nan

#
# Smooth
import mod_misc1 as mutil
nfltr = 5
if nfltr > 0:
  Aflt = RMU.copy()
  Aflt[np.where(np.isnan(Aflt))]=0.

  for iflt in range(nfltr):
    Aflt = mutil.box_fltr(Aflt,i1=520,i2=640,j1=1560,j2=1650, nbx=25)

  RMU = Aflt.copy()
  RMU[np.where(MSKr<1)] = 0.
  RMU[np.where(HH>=0)] = np.nan




f_save = True
if f_save:
# Relaxation files:
  pthrmu = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_paraD/relax/rmu/'
  flout  = 'rmu_{0}_{1:02d}day'.format(regn,int(rlx_days))
#  flout = 'rmu_{0:02d}days_Carib'.format(int(np.floor(rlx_days)))
  fina  = pthrmu + flout + '.a'
  finb  = pthrmu + flout + '.b'

# Write rmu.b:
  sln2 = 'Sulu Sea {0:d} > {1}m days e-folding time\n \n \n'.\
           format(int(np.floor(rlx_days)),int(np.floor(np.abs(zz0))))
  sln5 = 'i/jdm = {0} {1}\n'.format(IDM,JDM)
  sln6 = '     rmu: range =    {0:10.7E}   {1:10.7E}\n'.\
         format(np.nanmin(RMU), np.nanmax(RMU))

  print('Saving '+finb)
  with open(finb,'w') as fidb:
    fidb.write('Climatology relaxation mask\n')
    fidb.write(sln2)
    fidb.write(sln5)
    fidb.write(sln6)  

# Write rmu.a:
  print('Saving '+fina)
  with open(fina,'wb') as fida:
    A2D = RMU.squeeze().astype('>f4')
    A2D[np.where(np.isnan(A2D))] = 0.
    A2D.tofile(fida, format=">f4")
    toto.tofile(fida, format=">f4")


clrmp = copy(plt.cm.gist_ncar_r)
clrmp.set_bad(color=[0.7,0.7,0.7])

f_pltrmu = True
if f_pltrmu:
  fgnmb = 1
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
  im1 = ax1.pcolormesh(RMU, vmin=0, vmax=3.86e-6, cmap=clrmp)
  ax1.axis('scaled')
  ax1.set_xlim([450, 700])
  ax1.set_ylim([1525,1700])
  lvl = np.arange(-5000,0,1000)
  ax1.contour(HH,levels=lvl,linestyles='solid',colors=[(0.5,0.5,0.5)])

  stl='rmu, max rlx={0:6.3f} min={1:6.3f} day'.\
       format(rlx_days,rlx_dayw)
  ax1.set_title(stl)

# Colorbar:
  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='max')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  tval = clb.get_ticks()
# Convert to days-1
  dval = tval.copy()
  clb_lbl = []
  for ll in range(tval.shape[0]):
    vv = tval[ll]
    if vv==0:
      vdd = 0.
    else:
      vdd = 1./vv*(1./86400.)

    dval[ll] = vdd
    dvals = '{0:4.2f}'.format(vdd)
    clb_lbl.append(dvals)

  clb.set_ticks(tval)
  clb.set_ticklabels(clb_lbl)

  btx = 'create_rmu_sulu.py'
  mufig.bottom_text(btx)






