"""
  Utility subroutines for validation mom6
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pdb
import importlib
#import struct
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
from mod_utils_fig import bottom_text


def ocean_straits():
  STR = {
    "Fram80"  : {
      "xl1"   : 3379,
      "xl2"   : 3519,
      "yl1"   : 3007,
      "yl2"   : 3007
    },
    "Fram79"  : {
      "xl1"   : 3353,
      "xl2"   : 3555,
      "yl1"   : 2965,
      "yl2"   : 2965
    },
  }

  return STR

def minmax_clrmap(dmm,pmin=10,pmax=90,cpnt=0.01,fsym=False):
  """
  Find min/max limits for colormap 
  discarding pmin and 1-pmax min/max values
  cpnt - decimals to leave
  """
  dmm = dmm[~np.isnan(dmm)]
  a1  = np.percentile(dmm,pmin)
  a2  = np.percentile(dmm,pmax)
  cff = 1./cpnt
  rmin = cpnt*(int(a1*cff))
  rmax = cpnt*(int(a2*cff))

  if fsym and (rmin<0. and rmax>0.) :
    dmm = max([abs(rmin),abs(rmax)])
    rmin = -dmm
    rmax = dmm

  return rmin,rmax


def plot_xsect(XX, Hb, ZZ, A2d, HH, fgnmb=1, stl='Vert Section', rmin=[], rmax=[], \
               clrmp=[], btx=[], dcntr=2, lstart=25, zstart=-1000.,\
               ijsct=[]):
  """
  Plot vertical section of a scalar field
  """

  if not clrmp:
#  clrmp = copy(plt.cm.rainbow)
    clrmp = copy(plt.cm.Spectral_r)

  clrmp.set_bad(color=[0.7,0.7,0.7])

  nx = A2d.shape[1]
  ll = A2d.shape[0]

  if not rmin:
    rmin, rmax = minmax_clrmap(A2d)

  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
  #ax1.plot(XX,Hb)
  im1 = ax1.pcolormesh(XX, ZZ, A2d, shading='flat', \
                 cmap=clrmp,\
                 vmin=rmin, \
                 vmax=rmax)

  # Patch bottom:
  Bmin = np.floor(np.min(Hb))-100.
  verts = [(np.min(XX),-8000),*zip(XX,Hb),(np.max(XX),-8000)]
  poly = Polygon(verts, facecolor='0.6', edgecolor='0.6')
  ax1.add_patch(poly)

  # Plot interfaces
  dZZ = abs(np.diff(ZZ, axis=0))
  for kk in range(lstart, ll+1):
    dmm = ZZ[kk,:]
    if kk > 0:
      dzz = dZZ[kk-1,:]
#      if all (dzz == 0):  # all layers at the bottom
#        break
    ax1.plot(XX, dmm, color=[0,0,0], linewidth=0.5)

  iB = np.where(Hb == np.min(Hb))[0][0]
  if dcntr>0:
    cc=-1
    for kk in range(ll+1):
      zlr = ZZ[kk,iB]
      if kk > 0:
        dzz = dZZ[kk-1,iB]
      else:
        dzz = 0.

      if zlr <= zstart:
        cc += 1
        if cc%dcntr == 0 and dzz > 1.e-18:
          zm = zlr+0.5*dzz
          txt = '{0}'.format(kk)
          ax1.text(XX[iB],zm,txt,fontsize=10)

  ax1.set_xlim([np.min(XX), np.max(XX)])
  ax1.set_ylim([Bmin, 0])

# indices in Fortran convention:
#  stl = ('{5}, {4}, Hybrd layers, lat={0:5.2f}, jF={1}, iF={2}:{3}'.\
#         format(y0,j1+1,i1+1,i2+1,drfnm,xsct))
  ax1.set_title(stl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                       0.02, ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

# Draw a section on the map:
  xdm = HH.shape[1]
  ydm = HH.shape[0]
  if ijsct:
    i1 = ijsct[0]
    i2 = ijsct[1]
    j1 = ijsct[2]
    j2 = ijsct[3]
    ax2 = plt.axes([0.75, 0.025, 0.18, 0.18])
    ax2.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
    ax2.plot((i1,i2),(j1,j1),'-',color=[1.,0.4,0])
    dii = 250
    ilim1 = max([i1-dii,0])
    ilim2 = min([i2+dii,xdm])
    jlim1 = max([j1-dii,0])
    jlim2 = min([j1+dii,ydm])

    ax2.set_xlim([ilim1,ilim2])
    ax2.set_ylim([jlim1,jlim2])
    ax2.set_xticks([])
    ax2.set_yticks([])

  if btx:
    bottom_text(btx,pos=[0.02, 0.03])

  return 

class TRANSP():
  def __init__(self, dnmb, trnsp1d, XX, YY, Hbtm):
    trnsp1d    = np.expand_dims(trnsp1d, axis=0)
    self.TM    = np.array([dnmb])
    self.trnsp = trnsp1d
    self.xcrd  = XX
    self.ycrd  = YY
    self.Hbtm  = Hbtm

  def add_array(self, dnmb, trnsp1d):
    trnsp1d    = np.expand_dims(trnsp1d, axis=0)
    self.TM    = np.append(self.TM, dnmb)
    self.trnsp = np.append(self.trnsp, trnsp1d, axis=0)
 
class FIELD2D():
  def __init__(self, dnmb, XX, YY, LSgm, Hb, A2d):
    A2d        = np.expand_dims(A2d, axis=(0))
    self.TM    = np.array([dnmb])
    self.Fld2D = A2d
    self.LON   = XX 
    self.LAT   = YY
    self.Lsgm  = Lsgm
    self.Hbtm  = Hb

  def add_array(self, dnmb, A2d):
    A2d        = np.expand_dims(A2d, axis=0)
    self.TM    = np.array([dnmb])
    self.Fld2D = np.append(self.Fld2D, A2d, axis=0)


