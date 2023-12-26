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
      "smin"  : 32.5,
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
      "smin"  : 34.8,
      "smax"  : 36.7,
      "tmin"  : 3.,
      "tmax"  : 28.,
      "umin"  : -0.8,
      "umax"  : 0.8
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
               IJs=[], f_intrf=True, shad='flat', cntr1=[], cntr2=[]):
  """
  Plot vertical section of a scalar field
  IJs = I and J indices of the segment vertices or all points
  XX - arrays of along-section coordinates (distances, or lon/lat, ...)
  shad - shading: flat, gouraud
  cntr1, cntr2 - plotting contours
  """
  class nf(float):
    def __repr__(self):
      s = f'{self:.2f}'
      return f'{self:.0f}' if s[-1] == '0' else s

  try: 
    dmm = clrmp.colors
  except:
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


  # Patch bottom:
  Bmin = np.floor(np.min(Hb))-100.
  verts = [(np.min(XX),-8000),*zip(XX,Hb),(np.max(XX),-8000)]
  poly = Polygon(verts, facecolor='0.6', edgecolor='0.6')
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
  xdm = HH.shape[1]
  ydm = HH.shape[0]
  if len(IJs)>0:
    i1 = np.min(IJs[:,0])
    i2 = np.max(IJs[:,0])
    j1 = np.min(IJs[:,1])
    j2 = np.max(IJs[:,1])
    ax2 = plt.axes([0.74, 0.010, 0.18, 0.18])
    ax2.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
    ax2.plot(IJs[:,0],IJs[:,1],'-',color=[1.,0.4,0])
    dii = 100
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


def interp_2Dsect(A2d, ZZi, ZM2d, Hb):
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

    zz = ZM2d[:,isgm].squeeze().copy()
    zz = minterp.make_monotonic(zz)
    aa = A2d[:, isgm].squeeze().copy()
# To avoid problems with interpolation at the end-points (Gibbs effect
# for polynom > 1 degree): keep same values below bottom
# Or use linear interp 
#    print(f"isgm={isgm} zz={zz} hb0={hb0}")
    dbtm = abs(zz-hb0)
#    dbtm = np.diff(zz) 
    izb = min(np.where(dbtm <= 0.01)[0])
    aa[izb+1:] = aa[izb]

# Add extra point for 0 interp at the surface:
    zz = np.insert(zz, 0, 5.)
    aa = np.insert(aa, 0, aa[0])

# Add extra point at the bottom for interp deep z levels:
    zz = np.insert(zz, -1, zz[-1])
    zz[-1] = -10000.
    aa = np.insert(aa, -1, aa[-1])
#

# Depth levels:
#    print(f"isgm={isgm}")
    Aint = np.zeros((nintp))
    for kk in range(nintp):
# 2nd degree polynom gives errors near the bottom when profile makes jumps
# to random bottom values, use 1st degree
#      aai = minterp.pcws_lagr2(zz, aa, ZZi[kk])
      aai = minterp.pcws_lagr1(zz, aa, ZZi[kk])
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

def segm_half_coord(II, JJ, Vnrm1, Vnrm2, hLsgm1, hLsgm2, XX, YY, Hb):
  """
    Define coordinates, segment legnths, layer depths for 
    half segments
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
# indices for 1st half segment
    ix1 = ik*2
    ix2 = ix1+1
    II_hf[ix1] = II[ik] - 0.25*abs(Vnrm1[ik,1]) # 0 if vertical half-segm
    JJ_hf[ix1] = JJ[ik] - 0.25*abs(Vnrm1[ik,0]) # 0 if horiz half-segm
    II_hf[ix2] = II[ik] + 0.25*abs(Vnrm2[ik,1]) # 0 if vertical half-segm
    JJ_hf[ix2] = JJ[ik] + 0.25*abs(Vnrm2[ik,0]) # 0 if horiz half-segm

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
  Irmv = np.where( (dJ>0) & (dI==0.) )[0]   
  nIrmv= len(Irmv)
  if nIrmv == 0:
    print('Section is on X axis, no projection needed')
    return Favi

  Indx1   = []
  Indx2   = []
  f_gap = False
  for irr in range(nIrmv):
    ir1 = Irmv[irr]
    if ~f_gap:
      Indx1.append(ir1-1)
      f_gap = True
    
    if irr == nIrmv-1:  # last Y segm
      Indx2.append(ir1+2)
      f_gap = False
    else:
      ir2 = Irmv[irr+1]
      if (ir2-ir1) > 1:  # not the same Y-segment
        Indx2.append(ir1+2)
        f_gap = False      

# Case when Y-segm at the beginning
# Interpolate from next X-segment
  if Indx1[0] < 0: Indx1[0] = Indx2[0]

# Now check the last Y-segment - if Y segm is the last
# interpolate from the previous X-segment
  if Indx2[-1] > nI:
    Indx2[-1] = Indx1[-1]

# Interpolate into gaps:
  ngp = len(Indx1)
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



