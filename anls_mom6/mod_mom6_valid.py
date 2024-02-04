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
      "nlegs" : 1,
      "xl1"   : [3379],
      "xl2"   : [3519],
      "yl1"   : [3007],
      "yl2"   : [3007],
      "ucntr1": [-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35],
      "tcntr" : [x/10 for x in range(-10, 80, 10)],
      "scntr" : [x/10 for x in range(320, 358, 2)],
      "smin"  : 32.5,
      "smax"  : 35.0,
      "tmin"  : -1.5,
      "tmax"  : 5.5,
      "umin"  : -0.1,
      "umax"  : 0.1
    },
    "Fram79"  : {
      "nlegs" : 1,
      "xl1"   : [3353],
      "xl2"   : [3555],
      "yl1"   : [2965],
      "yl2"   : [2965],
      "ucntr1": [-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35],
      "tcntr" : [x/10 for x in range(-10, 80, 10)],
      "scntr" : [x/10 for x in range(320, 358, 2)],
      "smin"  : 32.5,
      "smax"  : 35.0,
      "tmin"  : -1.5,
      "tmax"  : 5.5,
      "umin"  : -0.1,
      "umax"  : 0.1
    },
    "Fram79s2": {
      "nlegs" : 1,
      "xl1"   : [3355],
      "xl2"   : [3530],
      "yl1"   : [2949],
      "yl2"   : [2989],
      "ucntr1": [-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35],
      "tcntr" : [x/10 for x in range(-10, 80, 10)],
      "scntr" : [x/100 for x in range(3480, 3520, 5)],
      "smin"  : 33.5,
      "smax"  : 35.0,
      "tmin"  : -1.5,
      "tmax"  : 5.5,
      "umin"  : -0.1,
      "umax"  : 0.1
    },
    "DavisStr": {
      "nlegs" : 1,
      "xl1"   : [2903],
      "xl2"   : [3006],
      "yl1"   : [2707],
      "yl2"   : [2707],
      "ucntr1": [-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35],
      "scntr" : [x/10 for x in range(348, 368, 2)],
      "smin"  : 32.2,
      "smax"  : 35.0,
      "tmin"  : -1.5,
      "tmax"  : 7.5,
      "umin"  : -0.1,
      "umax"  : 0.1
    },
    "DavisStr2": {
      "nlegs" : 1,
      "xl1"   : [2921],
      "xl2"   : [3009],
      "yl1"   : [2724],
      "yl2"   : [2724],
      "ucntr1": [-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35],
      "scntr" : [x/10 for x in range(348, 368, 2)],
      "smin"  : 32.2,
      "smax"  : 35.0,
      "tmin"  : -1.5,
      "tmax"  : 7.5,
      "umin"  : -0.1,
      "umax"  : 0.1
    },
    "DavisS2" : {
      "nlegs" : 1,
      "xl1"   : [2924],
      "xl2"   : [3004],
      "yl1"   : [2730],
      "yl2"   : [2705],
      "ucntr1": [-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35],
      "tcntr" : [x/10 for x in range(-10, 80, 10)],
      "scntr" : [x/10 for x in range(320, 358, 2)],
      "smin"  : 32.2,
      "smax"  : 35.0,
      "tmin"  : -1.5,
      "tmax"  : 7.5,
      "umin"  : -0.1,
      "umax"  : 0.1
    },
    "Yucatan" : {
      "nlegs" : 2,
      "xl1"   : [2485, 2485],
      "xl2"   : [2485, 2514],
      "yl1"   : [1779, 1785],
      "yl2"   : [1785, 1785],
      "ucntr1": [-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.25,0.5,0.75,1.,1.25,1.5,1.75],
      "scntr" : [x/10 for x in range(348, 368, 2)],
      "tcntr" : [x/10 for x in range(40, 280, 20)],
      "smin"  : 34.8,
      "smax"  : 36.7,
      "tmin"  : 3.,
      "tmax"  : 28.,
      "umin"  : -0.8,
      "umax"  : 0.8
    },
    "Yucatan2": {
      "nlegs" : 1,
      "xl1"   : [2487],
      "xl2"   : [2512],
      "yl1"   : [1778],
      "yl2"   : [1784],
      "ucntr1": [-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.25,0.5,0.75,1.,1.25,1.5,1.75],
      "scntr" : [x/10 for x in range(346, 368, 2)],
      "tcntr" : [x/10 for x in range(40, 280, 20)],
      "smin"  : 34.8,
      "smax"  : 36.7,
      "tmin"  : 3.,
      "tmax"  : 28.,
      "umin"  : -0.8,
      "umax"  : 0.8
    },
    "FlorCabl": {
      "nlegs" : 1,
      "xl1"   : [2571],
      "xl2"   : [2591],
      "yl1"   : [1856],
      "yl2"   : [1856],
      "ucntr1": [x/100 for x in range(-150, 0, 10)],
      "ucntr2": [x/100 for x in range(0, 200, 20)],
      "scntr" : [x/10 for x in range(346, 368, 2)],
      "tcntr" : [x/10 for x in range(40, 280, 20)],
      "smin"  : 34.8,
      "smax"  : 36.7,
      "tmin"  : 3.,
      "tmax"  : 28.,
      "umin"  : -1.5,
      "umax"  : 1.5
    },
    "BarentsS": {
      "nlegs" : 2,
      "xl1"   : [3587, 3634],
      "xl2"   : [3634, 3723],
      "yl1"   : [2938, 2884],
      "yl2"   : [2884, 2820],
      "ucntr1": [-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,.3,0.35,0.4,0.45],
      "scntr" : [x/10 for x in range(340, 360, 2)],
      "tcntr" : [x/10 for x in range(-10, 120, 10)],
      "smin"  : 34.1,
      "smax"  : 35.1,
      "tmin"  : -1,
      "tmax"  : 8.,
      "umin"  : -0.2,
      "umax"  : 0.2
    },
    "BeringS" : {
      "nlegs" : 2,
      "xl1"   : [1393, 1403],
      "xl2"   : [1403, 1415],
      "yl1"   : [2626, 2619],
      "yl2"   : [2619, 2619],
      "ucntr1": [-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [x/100 for x in range(0, 50, 5)],
      "scntr" : [x/10 for x in range(346, 368, 2)],
      "tcntr" : [x/10 for x in range(40, 280, 20)],
      "smin"  : 34.8,
      "smax"  : 36.7,
      "tmin"  : 3.,
      "tmax"  : 28.,
      "umin"  : -0.4,
      "umax"  : 0.4
    },
    "DenmarkS": {
      "nlegs" : 1,
      "xl1"   : [3280],
      "xl2"   : [3301],
      "yl1"   : [2679],
      "yl2"   : [2601],
      "ucntr1": [-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,.3,0.35,0.4,0.45],
      "scntr" : [33.0,33.2,33.4,33.6,33.8,34.0,34.2,\
                 34.4,34.6,34.8,34.9,34.91,34.92,34.93,\
                 34.94,34.95,34.96,34.97,34.98,34.99],
      "tcntr" : [x/100 for x in range(-100, 800, 25)],
      "smin"  : 33.4,
      "smax"  : 35.1,
      "tmin"  : -1.,
      "tmax"  : 7.,
      "umin"  : -0.2,
      "umax"  : 0.2
    },
    "IclShtl" : {
      "nlegs" : 1,
      "xl1"   : [3384],
      "xl2"   : [3541],
      "yl1"   : [2553],
      "yl2"   : [2460],
      "ucntr1": [-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,.3,0.35,0.4,0.45],
      "scntr" : [x/10 for x in range(340, 368, 2)],
      "tcntr" : [x/10 for x in range(-10, 160, 10)],
      "smin"  : 34.0,
      "smax"  : 35.3,
      "tmin"  : -1.,
      "tmax"  : 10.,
      "umin"  : -0.2,
      "umax"  : 0.2
    },
    "ShtlScot": {
      "nlegs" : 1,
      "xl1"   : [3543],
      "xl2"   : [3526],
      "yl1"   : [2454],
      "yl2"   : [2424],
      "ucntr1": [-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,.3,0.35,0.4,0.45],
      "scntr" : [x/10 for x in range(346, 368, 2)],
      "tcntr" : [x/10 for x in range(40, 280, 20)],
      "smin"  : 34.8,
      "smax"  : 36.7,
      "tmin"  : 3.,
      "tmax"  : 28.,
      "umin"  : -0.2,
      "umax"  : 0.2
    },
    "LaManch" : {
      "nlegs" : 1,
      "xl1"   : [3587],
      "xl2"   : [3593],
      "yl1"   : [2255],
      "yl2"   : [2248],
      "ucntr1": [-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,.3,0.35,0.4,0.45],
      "scntr" : [x/100 for x in range(3420, 3530, 2)],
      "tcntr" : [x/10 for x in range(80, 200, 10)],
      "smin"  : 34.6,
      "smax"  : 35.0,
      "tmin"  : 10.,
      "tmax"  : 15.,
      "umin"  : -0.2,
      "umax"  : 0.2
    },
    "NAtl39"  : {
      "nlegs" : 1,
      "xl1"   : [2632],
      "xl2"   : [3457],
      "yl1"   : [2030],
      "yl2"   : [2030],
      "ucntr1": [-0.2,-0.15,-0.1,-0.05],
      "ucntr2": [0,0.05,0.1,0.15,0.2,0.25,.3,0.35,0.4,0.45],
      "scntr" : [x/10 for x in range(346, 368, 2)],
      "tcntr" : [x/10 for x in range(40, 280, 20)],
      "smin"  : 34.8,
      "smax"  : 36.7,
      "tmin"  : 3.,
      "tmax"  : 28.,
      "umin"  : -0.2,
      "umax"  : 0.2
    },
  }

  return STR

def ocean_sections():
  XSCT = {
    "BaffNAFram" : {
      "NP"       : False,
      "II"       : [3052, 2999, 2946, 2929, 2996, 3102, 3217, 3395, \
                    3446, 3399, 3476, 3403],
      "JJ"       : [3086, 3056, 2762, 2546, 2368, 2291, 2287, 2486, \
                    2571, 2675, 2965, 3097],
      "ucntr1"   : [-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05],
      "ucntr2"   : [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35],
      "tcntr"    : [x/10 for x in range(-10, 80, 10)],
      "scntr"    : [x/10 for x in range(348, 358, 1)],
      "smin"     : 33.5,
      "smax"     : 35.5,
      "tmin"     : -1.,
      "tmax"     : 12.,
      "umin"     : -0.1,
      "umax"     : 0.1
    },
    "AlaskaIcld" : {
      "NP"       : True,
      "II"       : [1612, 1056, 3399, 3449, 3486, 3349], 
      "JJ"       : [2878, 3296, 3296, 3004, 2898, 2579], 
      "ucntr1"   : [-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05],
      "ucntr2"   : [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35],
      "tcntr"    : [x/10 for x in range(-20, 80, 10)],
      "scntr"    : [x/10 for x in range(340, 358, 1)],
      "smin"     : 33.1,
      "smax"     : 35.1,
      "tmin"     : -1.,
      "tmax"     : 5.,
      "umin"     : -0.1,
      "umax"     : 0.1
    },
    "GoMCarib"   : {
      "NP"       : False,
      "II"       : [2347, 2493, 2506, 2572, 2594, 2810], 
      "JJ"       : [1810, 1841, 1736, 1737, 1685, 1691], 
      "ucntr1"   : [-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05],
      "ucntr2"   : [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35],
      "tcntr"    : [x/10 for x in range(-20, 80, 10)],
      "scntr"    : [x/10 for x in range(340, 358, 1)],
      "smin"     : 34.8,
      "smax"     : 37.0,
      "tmin"     : 5.,
      "tmax"     : 28.,
      "umin"     : -0.1,
      "umax"     : 0.1
    },
  }

  return XSCT

def ts_prof_regions():
  REGNS = {
    "AmundsAO"  : {
      "lon0"    : 75.,
      "lat0"    : 87,
      "dlon"    : 70.,
      "dlat"    : 3.
    },
    "NansenAO"  : {
      "lon0"    : 85.,
      "lat0"    : 81.,
      "dlon"    : 60.,
      "dlat"    : 4.
    },
    "MakarovAO" : {
      "lon0"    : 225.,
      "lat0"    : 87.,
      "dlon"    : 90.,
      "dlat"    : 3.
    },
    "CanadaAO"  : {
      "lon0"    : 200.,
      "lat0"    : 77.,
      "dlon"    : 20.,
      "dlat"    : 7.
    },
  }
  
  return REGNS 

def interp_zlevels():
  # Interpolate onto fixed z-levels:
  ZZi = np.concatenate((np.arange(0.,     -6.,    -1),
                        np.arange(-6.,    -20.,   -2),
                        np.arange(-20.,   -100.,  -5),
                        np.arange(-100.,  -250.,  -10),
                        np.arange(-250.,  -500.,  -25),
                        np.arange(-500.,  -700.,  -50),
                        np.arange(-700.,  -2000., -100),
                        np.arange(-2000., -5250., -250)))
  return ZZi

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
               btm_midpnt=False, \
               IJs=[], f_intrf=True, shad='flat', cntr1=[], cntr2=[]):
  """
  Plot vertical section of a scalar field
  IJs = I and J indices of the segment vertices or all points
  XX - arrays of along-section coordinates (distances, or lon/lat, ...)
  shad - shading: flat, gouraud
  cntr1, cntr2 - plotting contours
  btm_midpnt - plot bottom at the hlaf grid point to better
               represent the true depth in the grid point
  """
  class nf(float):
    def __repr__(self):
      s = f'{self:.2f}'
      return f'{self:.0f}' if s[-1] == '0' else s

  try: 
    dmm = clrmp.colors
  except:
    print(f"!!! Colormap provided doesnt have attribute colors, using default")
#  clrmp = copy(plt.cm.rainbow)
    clrmp = copy(plt.cm.Spectral_r)

  clrmp.set_bad(color=[0.98,0.98,0.98])

  nx = A2d.shape[1]
  ll = A2d.shape[0]

  if not rmin:
    rmin, rmax = minmax_clrmap(A2d)

  plt.ion()
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])
  #ax1.plot(XX,Hb)
  im1 = ax1.pcolormesh(XX, ZZ, A2d, shading=shad, \
                 cmap=clrmp,\
                 vmin=rmin, \
                 vmax=rmax)
  if len(cntr1) > 0:
    CS1 = ax1.contour(XX, ZZ, A2d, cntr1, colors=[(0.,0.6,0.9)], \
                linestyles='solid',linewidths=1)
    # Recast levels to new class
#    CS1.levels = [nf(val) for val in CS1.levels]
    fmt='%.2f'
    ax1.clabel(CS1, CS1.levels, inline=True, fmt=fmt, fontsize=10)


  if len(cntr2) > 0:
    CS2 = ax1.contour(XX, ZZ, A2d, cntr2, colors=[(0.6,0.6,0.6)], \
                linestyles='solid',linewidths=1)
    # Recast levels to new class
#    CS2.levels = [nf(val) for val in CS2.levels]
    fmt='%.2f'
    ax1.clabel(CS2, CS2.levels, inline=True, fmt=fmt, fontsize=10)

# Compute half-grid point coords:
  if btm_midpnt:
    XXhf = XX.copy()
    dXX  = np.diff(XX)
    dXX  = np.append(dXX,dXX[-1])
    XXhf = XX + 0.5*dXX

  # Patch bottom:
  Bmin = np.floor(1.05*np.min(Hb))
  if btm_midpnt:
    verts = [(np.min(XX),-8000),(np.min(XX),Hb[0]),*zip(XXhf,Hb),\
             (np.max(XX),Hb[-1]),(np.max(XX),-8000)]
  else:
    verts = [(np.min(XX),-8000),*zip(XX,Hb),(np.max(XX),-8000)]
  poly = Polygon(verts, facecolor='0.6', edgecolor='0.6', zorder=10)
  ax1.add_patch(poly)

  # Plot interfaces only if zz is 2D
  if len(ZZ.shape) == 1:
    f_intrf = False

  if f_intrf:
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
          if cc%dcntr == 0 and dzz > 0.1:
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
  if len(IJs)>0:
    xdm = HH.shape[1]
    ydm = HH.shape[0]
    i1 = np.min(IJs[:,0])
    i2 = np.max(IJs[:,0])
    j1 = np.min(IJs[:,1])
    j2 = np.max(IJs[:,1])
    ax2 = plt.axes([0.74, 0.010, 0.18, 0.18])
    ax2.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
    if np.min(HH) <= -1000.:
      ax2.contour(HH,[-3000,-2000,-1000],colors=[(0.7,0.7,0.7)],\
                linestyles='solid',linewidths=1)
    ax2.plot(IJs[:,0],IJs[:,1],'.',ms=8,color=[1.,0.4,0])
    dii = 60
    ilim1 = max([i1-dii,0])
    ilim2 = min([i2+dii,xdm])
    jlim1 = max([j1-dii,0])
    jlim2 = min([j2+dii,ydm])

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
  def __init__(self, dnmb, XX, YY, LSgm, ZZi, Hb, A2d):
    A2d        = np.expand_dims(A2d, axis=(0))
    self.TM    = np.array([dnmb])
    self.Fld2D = A2d
    self.LON   = XX 
    self.LAT   = YY
    self.Lsgm  = LSgm
    self.Hbtm  = Hb
    self.ZZi   = ZZi

  def add_array(self, dnmb, A2d):
    A2d        = np.expand_dims(A2d, axis=0)
    self.TM    = np.append(self.TM, dnmb)
    self.Fld2D = np.append(self.Fld2D, A2d, axis=0)

class UTS2D():
  def __init__(self, dnmb,II,JJ, XX, YY, LSgm, ZZi, Hb, A2d):
    A2d        = np.expand_dims(A2d, axis=(0))
    self.TM    = np.array([dnmb])
    self.Fld2D = A2d
    self.Iindx = II
    self.Jindx = JJ
    self.LON   = XX 
    self.LAT   = YY
    self.Lsgm  = LSgm
    self.Hbtm  = Hb
    self.ZZi   = ZZi

  def add_array(self, dnmb, A2d):
    A2d        = np.expand_dims(A2d, axis=0)
    self.TM    = np.append(self.TM, dnmb)
    self.Fld2D = np.append(self.Fld2D, A2d, axis=0)

def interp_2Dsect_segmhalf(A2d, ZZi, ZZ2d, ZM2d, Hb):
  """
  Interpolate 2D section from hybrid to z-fixed levels
   for half segments
  """
  import mod_interp1D as minterp

  nsgm  = A2d.shape[1]
  nlvls = A2d.shape[0]
  nintp = ZZi.shape[0]
  A2di  = np.zeros((nintp, nsgm)) + 1.e30
  for ihf in range(nsgm):
    isgm = int(np.floor(ihf/2))
    hb0 = Hb[isgm]
    if hb0 >= 0.:
#      print(f"isgm={isgm} ihf={ihf}")
      A2di[:,ihf] = np.nan
      continue

    zz = ZZ2d[:,isgm].squeeze().copy()
    zm = ZM2d[:,isgm].squeeze().copy()
    zm = minterp.make_monotonic(zm)
    aa = A2d[:, ihf].squeeze().copy()
# To avoid problems with interpolation at the end-points (Gibbs effect
# for polynom > 1 degree): keep same values below bottom
# Or use linear interp 
#    print(f"isgm={isgm} zz={zz} hb0={hb0}")
    dbtm = abs(zz-hb0)
#    dbtm = np.diff(zz)
    try:
      izb = min(np.where(dbtm <= 0.01)[0]) # interf. depths!
    except:
      print(f"interp_2Dsect: couldnt find bottom: isgm={isgm} " +\
            f"Hbottom={hb0} ZZ={zz[-1]}")
      raise Exception("Failed locating bottom")

    if izb <= nlvls:
      aa[izb:] = aa[izb-1]

# Add extra point for 0 interp at the surface:
    zm = np.insert(zm, 0, 5.)
    aa = np.insert(aa, 0, aa[0])

# Add extra point at the bottom for interp deep z levels:
    zm = np.insert(zm, -1, zm[-1])
    zm[-1] = -10000.
    aa = np.insert(aa, -1, aa[-1])
#

# Depth levels:
#    print(f"isgm={isgm}")
    Aint = np.zeros((nintp))
    for kk in range(nintp):
# 2nd degree polynom gives errors near the bottom when profile makes jumps
# to random bottom values, use 1st degree
#      aai = minterp.pcws_lagr2(zm, aa, ZZi[kk])
      aai = minterp.pcws_lagr1(zm, aa, ZZi[kk])
      Aint[kk] = aai
      A2di[kk,ihf] = aai

  return A2di


def interp_2Dsect(A2d, ZZi, ZZ2d, ZM2d, Hb):
  """
  Interpolate 2D section from hybrid to z-fixed levels
  """
  import mod_interp1D as minterp

  nsgm  = A2d.shape[1]
  nlvls = A2d.shape[0]
  nintp = ZZi.shape[0]
  A2di  = np.zeros((nintp, nsgm)) + 1.e30
  for isgm in range(nsgm):
    hb0 = Hb[isgm]
    if hb0 >= 0.:
      A2di[:,isgm] = np.nan
      continue

    zz = ZZ2d[:,isgm].squeeze().copy()
    zm = ZM2d[:,isgm].squeeze().copy()
    zm = minterp.make_monotonic(zm)
    aa = A2d[:, isgm].squeeze().copy()
# To avoid problems with interpolation at the end-points (Gibbs effect
# for polynom > 1 degree): keep same values below bottom
# Or use linear interp 
#    print(f"isgm={isgm} zz={zz} hb0={hb0}")
    dbtm = abs(zz-hb0)
#    dbtm = np.diff(zz)
    try:
      izb = min(np.where(dbtm <= 0.01)[0]) # interf. depths!
    except:
      print(f"interp_2Dsect: couldnt find bottom: isgm={isgm} " +\
            f"Hbottom={hb0} ZZ={zz[-1]}")
      raise Exception("Failed locating bottom")

    if izb <= nlvls:
      aa[izb:] = aa[izb-1]

# Add extra point for 0 interp at the surface:
    zm = np.insert(zm, 0, 5.)
    aa = np.insert(aa, 0, aa[0])

# Add extra point at the bottom for interp deep z levels:
    zm = np.insert(zm, -1, zm[-1])
    zm[-1] = -10000.
    aa = np.insert(aa, -1, aa[-1])

# Depth levels:
#    print(f"isgm={isgm}")
    Aint = np.zeros((nintp))
    for kk in range(nintp):
# 2nd degree polynom gives errors near the bottom when profile makes jumps
# to random bottom values, use 1st degree
#      aai = minterp.pcws_lagr2(zz, aa, ZZi[kk])
      aai = minterp.pcws_lagr1(zm, aa, ZZi[kk])
      Aint[kk] = aai
      A2di[kk,isgm] = aai

  return A2di

def segm_half_zintrf(II, JJ, ZZ2d):
  nI      = len(II)
  kdm     = ZZ2d.shape[0]
  ZZ_hf   = np.zeros((kdm,2*nI))-999.
  for ik in range(nI):
# indices for 1st half segment
    ix1 = ik*2
    ix2 = ix1+1
# Interf depths:
    if ik < nI-1:
      ZZ_hf[:,ix1] = ZZ2d[:,ik] + 0.25*(ZZ2d[:,ik+1]-ZZ2d[:,ik])      
      ZZ_hf[:,ix2] = ZZ2d[:,ik] + 0.75*(ZZ2d[:,ik+1]-ZZ2d[:,ik])      
    else:
      ZZ_hf[:,ix1] = ZZ2d[:,ik]
      ZZ_hf[:,ix2] = ZZ2d[:,ik]

  return ZZ_hf

def segm_half_coord(II, JJ, hLsgm1, hLsgm2, XX, YY, Hb):
  """
    Define coordinates, segment legnths, layer depths for 
    half segments

    The code should work for either left/right up/down direction 
    of the segments (i.e. start from Istart < Iend, or Istart > Iend, same 
    for J)

    Has been tested for Istart < Iend (i.e. going from West to East) both
    for J north-south and south-north 
  """
# Half-segment indices and coordinates for plotting:
  nI      = len(II)
  II_hf   = np.zeros((2*nI))-999.
  JJ_hf   = np.zeros((2*nI))-999.
  XX_hf   = np.zeros((2*nI))-999.
  YY_hf   = np.zeros((2*nI))-999.
  LSgm_hf = np.zeros((2*nI))-999.
  Hb_hf   = np.zeros((2*nI))-999.
  for ik in range(nI):
# Define direction of the half-segments: >0 or <0:
# dI = 0 if vertical half-segm, dJ = 0 for horizontal half-segments
# Left half-segments wrt to grid cell 
    if ik > 0:
      dI1 = II[ik] - II[ik-1]
      dJ1 = JJ[ik] - JJ[ik-1]
    else:
      dI1 = 0
      dJ1 = 0
# Right half-segments:
    if ik < nI-1:
      dI2 = II[ik+1] - II[ik]
      dJ2 = JJ[ik+1] - JJ[ik]
    else:
      dI2 = 0
      dJ2 = 0 

# indices for 1st half segment
    ix1 = ik*2
    ix2 = ix1+1
    II_hf[ix1] = II[ik] - 0.25*np.sign(dI1)
    JJ_hf[ix1] = JJ[ik] - 0.25*np.sign(dJ1)
    II_hf[ix2] = II[ik] + 0.25*np.sign(dI2)
    JJ_hf[ix2] = JJ[ik] + 0.25*np.sign(dJ2)

#
#    II_hf[ix1] = II[ik] - 0.25*abs(Vnrm1[ik,1]) # 0 if vertical half-segm
#    JJ_hf[ix1] = JJ[ik] - 0.25*abs(Vnrm1[ik,0]) # 0 if horiz half-segm
#    II_hf[ix2] = II[ik] + 0.25*abs(Vnrm2[ik,1]) # 0 if vertical half-segm
#    JJ_hf[ix2] = JJ[ik] + 0.25*abs(Vnrm2[ik,0]) # 0 if horiz half-segm
#
# lengths:
    LSgm_hf[ix1] = hLsgm1[ik]
    LSgm_hf[ix2] = hLsgm2[ik]

# coordinates:
    if ik > 0:
      dX1 = XX[ik] - XX[ik-1]
      dY1 = YY[ik] - YY[ik-1]
    else:
      dX1 = 0.
      dY1 = 0.

    if ik < nI-1:
      dX2 = XX[ik+1] - XX[ik]
      dY2 = YY[ik+1] - YY[ik]
    else:
      dX2 = 0.
      dY2 = 0.

    XX_hf[ix1] = XX[ik] - 0.25*dX1
    XX_hf[ix2] = XX[ik] + 0.25*dX2
    YY_hf[ix1] = YY[ik] - 0.25*dY1
    YY_hf[ix2] = YY[ik] + 0.25*dY2

# Bottom:
    if ik < nI-1:
      Hb_hf[ix1] = Hb[ik] + 0.25*(Hb[ik+1]-Hb[ik])
      Hb_hf[ix2] = Hb[ik] + 0.75*(Hb[ik+1]-Hb[ik])
    else:
      Hb_hf[ix1] = Hb[ik]
      Hb_hf[ix2] = Hb[ik]

    if Hb[ik] > 0.:
      Hb_hf[ix1] = Hb[ik]
      Hb_hf[ix2] = Hb[ik]

  return II_hf, JJ_hf, XX_hf, YY_hf, Hb_hf, LSgm_hf

def project2X(II0, JJ0, Fav):
  """
  For plotting to avoid discontinuities caused by
  zigzaging segments of the sections
  project 2D section on X-axis reducing Y-oriented segments
  Fill gaps of Y-segments by interpolating between X-segm
  for smooth plotting
  """
# Find start-end indices for interpolation into gaps
# Want to eliminate Y-oriented segments that cause discontinuities
# in the U fields (here: nodes 2, 3, 4, as they have u-components
# normal to the section or part of the grid cell at the corners)
# interpolate 1 point before and after the zigzag point
#            4         5 
#            *---------*-------
#            |
#            |
#            * 3
#            |
#            |
# --*--------*  
#   1        2
#
#  Interpolate (1) - (5) to nodes (2,3,4)
#
  import mod_interp1D as mintrp
  Favi = Fav.copy()
  nI   = len(II0)
  dJ   = np.diff(JJ0)
  dI   = np.diff(II0)
  ndj  = len(dJ)
#  Irmv = np.where( (abs(dJ)>0) & (dI==0.) )[0]   
  Irmv = np.where( (abs(dJ)>0) )[0]   
  nIrmv= len(Irmv)
  if nIrmv == 0:
    print('Section is on X axis, no projection needed')
    return Favi

  Indx1   = []
  Indx2   = []
  f_gap = False
  for irr in range(ndj-1):
    if abs(dJ[irr]) == 0:
      continue
    else:
      if not f_gap:
# Start of the gap: point to use for interpolation over the gap
        Indx1.append(irr)
        f_gap = True

      if f_gap:
        if abs(dJ[irr+1]) == 0:
# End of the gap: 2nd pnt for interpolation:
          Indx2.append(irr+1)
          f_gap = False

# Case when Y-segm at the beginning
# Interpolate from next X-segment
  if Indx1[0] <= 0: Indx1[0] = Indx2[0]

# Now check the last Y-segment - if Y segm is the last
# interpolate from the previous X-segment
#  if Indx2[-1] > nI:
#    Indx2[-1] = Indx1[-1]
  if f_gap and abs(dJ[-1]) > 0:
    Indx2.append(Indx1[-1])
    f_gap = False

# Interpolate into gaps:
  ngp = len(Indx1)
  print(f"len Indx1={ngp} Indx2={len(Indx2)}")
  for ii in range(ngp):
    ix1 = Indx1[ii]
    ix2 = Indx2[ii]
    Xp  = np.array([JJ0[ix1],JJ0[ix2]])
    Yp  = np.array([Fav[:,ix1], Fav[:,ix2]]).transpose()

    for igp in range(ix1+1,ix2):
      xx0 = JJ0[igp]
      Ai  = mintrp.lagr_polynom1D(Xp, Yp, xx0)
      Favi[:,igp] = Ai

  return Favi

def runmn1D(AA, npnts=3, axis=0):
  """
    Running mean along axis = 0 - by rows along X axis, 
    = 1 - by columns along Y axis
 
    assumed equidistant nodes
    averaging is over npnts nodes, if npnts is even +1 to make odd

    AA is 2D array
  """
  if axis == 1:
    AA = AA.transpose()

  idm = AA.shape[1]
  jdm = AA.shape[0]

  if npnts%2 == 0:
    print(f"Specified # points for averaging even {npnts} ---> {npnts+1}")
    npnts = npnts+1

  dii = int((npnts-1)/2)

  Aav = AA.copy()
  for jj in range(jdm):
    dmm = AA[jj,:]
 
    for ii in range(idm):
      i1 = ii-dii
      i2 = ii+dii
      i1 = max([i1,0])
      i2 = min([i2,idm])
 
      Aav[jj,ii] = np.nanmean(dmm[i1:i2+1])

  if axis == 1:
    Aav = Aav.transpose()

  return Aav

def box_fltr(AA, npnts=3):
  """
    Running box-fltr 
    Nans are ignored
 
    assumed equidistant nodes
    averaging is over (npnts x npnts), if npnts is even +1 to make odd

    AA is 2D array
  """
  idm = AA.shape[1]
  jdm = AA.shape[0]

  if npnts%2 == 0:
    print(f"Specified # points for averaging even {npnts} ---> {npnts+1}")
    npnts = npnts+1

  dii = int((npnts-1)/2)

  Aav = AA.copy()
  for jj in range(jdm):
    for ii in range(idm):
      if np.isnan(AA[jj,ii]):
        continue
      i1 = ii-dii
      i2 = ii+dii
      i1 = max([i1,0])
      i2 = min([i2,idm])
 
      j1 = jj-dii
      j2 = jj+dii
      j1 = max([j1,0])
      j2 = min([j2,jdm])

      Aav[jj,ii] = np.nanmean(AA[j1:j2+1,i1:i2+1])

  return Aav

def plot_section_orthomap(II, JJ, IJ, LON, LAT, HH,\
                     lon0=-10, lat0=60, res='l',\
                     XI=[], dX=500., fgnmb=1, btx=[], \
                     sttl = 'section'):
  """
    Plot section line on orthonormal projection to show
    sections in the polar regions
    To change projection angle: set lon0, lat0
  """
  from mpl_toolkits.basemap import Basemap, cm
  import matplotlib.colors as colors
  import matplotlib.mlab as mlab
#  from matplotlib.colors import ListedColormap
  m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)

  xh, yh = m(LON,LAT)  # modl grid coordinates on the projections

  fig1 = plt.figure(fgnmb,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

  m.drawcoastlines()
  m.fillcontinents(color=(0.2,0.2,0.2))
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))
  m.contour(xh, yh, HH, [-4000,-3000,-2000,-1000], 
            colors=[(0.8,0.8,0.8)], 
            linestyles='solid')

  m.plot(xh[JJ,II],yh[JJ,II],'.')
  ax1.set_title(sttl)

# Show distance or other coordinates along the section:
  if len(XI) > 0:
    nX = 0
    while nX <= XI[-1]:
      dd  = np.abs(XI-nX)
      ix0 = np.argmin(dd)
      ii1 = II[ix0]
      jj1 = JJ[ix0]
      dh  = np.sqrt((xh[jj1,ii1] - xh[jj1+1,ii1])**2 \
                     + (xh[jj1,ii1] - xh[jj1,ii1+1])**2)
      xt  = xh[jj1,ii1] + 10.*dh 
      yt  = yh[jj1,ii1] + 10.*dh

      m.plot(xh[jj1,ii1], yh[jj1,ii1], '.', ms=11, color=(1.,0.,0.))
      ax1.text(xt, yt, f"{XI[ix0]:.0f}", fontsize=10, color=(1.,0.,0.))
      nX += dX

  if len(btx) > 0:
    bottom_text(btx,pos=[0.08, 0.08])

  return ax1

def plot_section_map(II, JJ, IJ, Vnrm1, Vnrm2, II_hf, JJ_hf, \
                     HH, XI=[], dX=500., fgnmb=1, btx='mod_mom6_valid.py', \
                     sttl = 'section'):
  """
    Show section line on the map
    with half-segments and norms
 
    If XI is not empty - array of along-section coordinates,
    indicate every dX point on the section

  """
  xl1 = np.min(IJ[:,0])
  xl2 = np.max(IJ[:,0])
  yl1 = np.min(IJ[:,1])
  yl2 = np.max(IJ[:,1])
  dxy = 50

  plt.ion()
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
  ax1.contour(HH,[0], 
              colors=[(0,0.,0)], 
              linestyles='solid')
  ax1.contour(HH,list(range(-5000,0,500)),
              colors=[(0.8,0.8,0.8)], 
              linestyles='solid')
  ax1.axis('scaled')
  ax1.set_xlim([xl1-dxy, xl2+dxy])
  ax1.set_ylim([yl1-dxy, yl2+dxy])
  ax1.plot(IJ[:,0],IJ[:,1],'b-')
  ax1.plot(IJ[:,0],IJ[:,1],'.', ms=17, color=(0.,0.5,0.8))

  if len(XI) > 0:
#    ax1.plot(II[0], JJ[0],'.', ms=18, color=(1.,0.,0.))
#    ax1.text(II[0]-40,JJ[0],f"0", fontsize=12, color=(1.,0.,0.))
    nX = 0
    while nX <= XI[-1]:
      dd  = np.abs(XI-nX)
      ii0 = np.argmin(dd)
      im1 = ii0 - 5
      ip1 = ii0 + 5
      im1 = max([0, im1])
      ip1 = min([ip1, len(II)-1])
      di  = II[ip1] - II[im1]
      dj  = JJ[ip1] - JJ[im1]
# Perpendicular direction:  
      pdi = -dj
      pdj = di
      lpd = np.sqrt(pdi*pdi + pdj*pdj)
      pdi = pdi/lpd
      pdj = pdj/lpd
      xt  = II[ii0] - pdi*50.
      yt  = JJ[ii0] - pdj*50.     
 
#      print(f"nX={nX} di={di} dj={dj} pdi={pdi} pdj={pdj}") 

      ax1.plot(II[ii0], JJ[ii0],'.', ms=18, color=(1.,0.,0.))
      ax1.text(xt, yt ,f"{XI[ii0]:.0f}", fontsize=10, color=(1.,0.,0.),
               horizontalalignment='center',
               verticalalignment='center')
      nX += dX 

  ax1.set_title(sttl)
#
# Show half-segments locations - should coincide with the norm vectors:
  if len(II_hf) > 0:
    ax1.plot(II_hf, JJ_hf, 'd', ms=5, color=(0.6,0.6,0.6))

  nsgm = len(II)
  for isgm in range(nsgm):
    ii0 = II[isgm]
    jj0 = JJ[isgm]
    if isgm < nsgm-1:
      ip1 = II[isgm+1]
      jp1 = JJ[isgm+1]
    else:
      ip1 = ii0
      jp1 = jj0
    if isgm > 0:
      im1 = II[isgm-1]
      jm1 = JJ[isgm-1]
    else:
      im1 = ii0
      jm1 = jj0
    ih1 = 0.5*(im1+ii0)
    jh1 = 0.5*(jm1+jj0)
    ih2 = 0.5*(ip1+ii0)
    jh2 = 0.5*(jp1+jj0)

    ax1.plot(ii0, jj0, '.', ms=12, color=(0., 0.2, 0.5))
    if isgm > 0:
      ax1.plot([ih1,ii0],[jh1,jj0],'r-') # 1st half of the segment
    if isgm < nsgm-1:
      ax1.plot([ii0,ih2],[jj0,jh2],'g-') # 2nd half
# Plot norm
    scl = 0.2
    v1  = scl*Vnrm1[isgm,:]
    v2  = scl*Vnrm2[isgm,:]
    in1 = 0.5*(ii0-ih1)+ih1
    jn1 = 0.5*(jj0-jh1)+jh1
    in2 = 0.5*(ih2-ii0)+ii0
    jn2 = 0.5*(jh2-jj0)+jj0
    if isgm > 0:
      ax1.plot([in1, in1+v1[0]],[jn1, jn1+v1[1]],'-', color=(0.7,0.2,0))
    if isgm < nsgm-1:
      ax1.plot([in2, in2+v2[0]],[jn2, jn2+v2[1]],'-', color=(0.,1.,0.4))

  bottom_text(btx,pos=[0.08, 0.08])

  return ax1

def read_flcable(YR, dflobs):
  """
    Read Fl transport from cable observations
    Data from https://www.aoml.noaa.gov/phod/floridacurrent/data_access.php
  """
  import mod_time as mtime

  fid = open(dflobs,'r')
# Find length of the input file:
  fid.seek(0,2)
  fend = fid.tell()
  dmm = 0
# Read header then data
  fid.seek(0)
  fpos = fid.tell()

  FLX = []
  TM0 = []
  while fpos < fend:
    dmm = fid.readline().split()
    if len(dmm)==0:
      continue
    elif len(dmm[0].strip()) == 0:
      print('Empty string, end ...')
      break
    elif dmm[0] == '%':
      continue

    yr0 = int(dmm[0])
    mm0 = int(dmm[1])
    dd0 = int(dmm[2])
    flx = float(dmm[3])
    flg = int(dmm[4])

    dnmb0 = mtime.datenum([yr0,mm0,dd0])
    TM0.append(dnmb0)
    if flg == 2:
      FLX.append(1.e30)
    else:
      FLX.append(flx)

    fpos = fid.tell() 

  fid.close()

  FLX = np.array(FLX)
  TM0 = np.array(TM0)

  return FLX, TM0

def read_flcableZ(dflin):
  """
    Read Zulema's transport estimates
  """
  fid = open(dflin,'r')
# Find length of the input file:
  fid.seek(0,2)
  fend = fid.tell()
  dmm = 0
  fid.seek(0)
  fpos = fid.tell()

  FLX = []
  while fpos < fend:
    dmm = fid.readline().split()
    if len(dmm)==0:
      continue
    elif len(dmm[0].strip()) == 0:
      print('Empty string, EOF done ...')
      break
    elif dmm[0] == '%':
      continue

    flx = float(dmm[0])
    FLX.append(flx)
    fpos = fid.tell()
 
  fid.close()
  FLX = np.array(FLX)

  return FLX


def plot_points_orthomap(X,Y, LON, LAT, HH, clr=[0,0,0.8],\
                     lon0=-10, lat0=60, res='l',\
                     fgnmb=1, btx=[], \
                     sttl = 'Points'):
  """
    Plot data points on orthonormal projection to show
    observations or something else in the polar regions
    X, Y - are point coordinates (lon, lat)
    To change projection angle: set lon0, lat0
  """
  from mpl_toolkits.basemap import Basemap, cm
  import matplotlib.colors as colors
  import matplotlib.mlab as mlab
#  from matplotlib.colors import ListedColormap
  m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)

  xh, yh = m(LON,LAT)  # modl grid coordinates on the projections

  plt.ion()
  fig1 = plt.figure(fgnmb,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

  m.drawcoastlines()
  m.fillcontinents(color=(0.2,0.2,0.2))
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))
  m.contour(xh, yh, HH, [-4000,-3000,-2000,-1000],
            colors=[(0.8,0.8,0.8)],
            linestyles='solid')

  xp, yp = m(X,Y)
  m.plot(xp,yp,'.',color=clr)
  ax1.set_title(sttl)

  if len(btx) > 0:
    bottom_text(btx,pos=[0.08, 0.03])

  return ax1


