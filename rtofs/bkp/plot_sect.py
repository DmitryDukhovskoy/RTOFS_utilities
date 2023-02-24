"""
  Plot sections, etc.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
#import matplotlib as mpl
#import struct
#from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
#import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
#from mpl_toolkits.basemap import Basemap, shiftgrid

sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython')
import mod_read_hycom
#importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo, read_hycom, read_topo
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text

def Jsections():
  """
    Dictionary of sections
    in W-E directions
    with i_start (west),i_end (east), j indices
  """
  SCT = {
    "SPac1" : [1512,2580,1000],
    "SPac2" : [2000,2600,745],
    "GoM1"  : [2340,2560,1798],
    "NAtl1" : [2900,3000,1897]
  }

  return SCT

def plot_Jsection(xsct,ZZ,HH,LON,LAT, drfnm, sct_show=[], \
                  fgnmb=5,dcntr=1, zstart=-1000., jlrs=-1):
# Plot J-section
# Write contours at the deepest point with dcntr intervals
# make dcntr=0 if not needed
# jrls >=0 indicate location where layer thickness is printed out
#
# Sections: i1,i2,j1 - W-E sections
#
# Input: xsct     - section name from dictionary SCT
#        ZZ       - interface depths 3D        
#        HH       - bathymetry 2D
#        LON/LAT  - model grid 2D
#        fgnmb    - figure # to plot
#        dcntr    - print layer # on xsection every dcntr layer
#        zstart   - if dcntr>0, where to start layer notation
#        jlrs     - show location (where layers are printed out, e.g.)
#        sct_show - show section location on map in figure(sct_show)
#
  SCT = Jsections()

  aa = SCT.get(xsct)
  i1 = aa[0]
  i2 = aa[1]
  j1 = aa[2]

  ZZs = ZZ[:,j1,i1:i2]
  Hb  = HH[j1,i1:i2]
  XX  = LON[j1,i1:i2]
  y0  = LAT[j1,i1]
  Bmin = np.floor(np.min(Hb))-100.
  dZZs = np.abs(np.diff(ZZs, axis=0))

  ll = ZZ.shape[0]-1 

  plt.ion()
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
  #ax1.plot(XX,Hb)
  # Patch bottom:
  verts = [(np.min(XX),-8000),*zip(XX,Hb),(np.max(XX),-8000)]
  poly = Polygon(verts, facecolor='0.6', edgecolor='0.6')
  ax1.add_patch(poly)

  # Plot interfaces
  for kk in range(ll+1):
    dmm = ZZs[kk,:]
    if kk > 0:
      dzz = dZZs[kk-1,:]
#      if all (dzz == 0):  # all layers at the bottom
#        break
    ax1.plot(XX,dmm,color=[0,0,0])

  iB = np.where(Hb == np.min(Hb))[0][0]
  if dcntr>0:
    cc=-1
    for kk in range(ll+1):
      zlr = ZZs[kk,iB]
      if kk > 0:
        dzz = dZZs[kk-1,iB]
      else:
        dzz = 0.

      if zlr <= zstart:
        cc += 1
        if cc%dcntr == 0 and dzz > 1.e-18:
          zm = zlr+0.5*dzz
          txt = '{0}'.format(kk)
          ax1.text(XX[iB],zm,txt,fontsize=10)

# Inidcate where layers are printed out:
  if jlrs >= 0:
    xlrs = XX[jlrs]
    ax1.plot((xlrs,xlrs),(Bmin,0),'--',color=[0.8,0,0])
  #ax1.autoscale(enable=True, axis='both', tight='True')
  ax1.set_xlim([np.min(XX), np.max(XX)])
  ax1.set_ylim([Bmin, 0])


  stl = ('{5}, {4}, Hybrd layers, lat={0:5.2f}, j={1}, i={2}:{3}'.\
         format(y0,j1+1,i1+1,i2+1,drfnm,xsct))
  ax1.set_title(stl)

  btx = 'find_collapsed_lrs.py'
  bottom_text(btx)

#  if jp0 >= 0:
#    print_1D(dZZs[:,jp0])

#
# Draw a section on the map:
  if sct_show:
    sct_show.plot((i1,i2),(j1,j1),'-',color=[0.,0.4,1])

  return ZZs, Hb, XX, dZZs

def plot_rho_Jsct(xsct,Rho,ZZ0,HH,LON,LAT, stl=[], btx=[], sct_show=[], \
                  fgnmb=5, rmin=[], rmax=[], clrmp=[]):
# pcolormesh in python:
# https://matplotlib.org/stable/gallery/images_contours_and_fields/pcolormesh_grids.html
#
# Plot rho fields along J-section
# Write contours at the deepest point with dcntr intervals
# make dcntr=0 if not needed
# jrls >=0 indicate location where layer thickness is printed out
#
# Sections: i1,i2,j1 - W-E sections
#
# Input: xsct     - section name from dictionary SCT
#        Rho      - 3D field to plot (density)
#        ZZ0      - z-level interface depths 1D array or 3D
#                   if 1D constant Z is assumed for the section
#                   3D - Z layers extracted for the section
#        HH       - bathymetry 2D
#        LON/LAT  - model grid 2D
#        fgnmb    - figure # to plot
#        dcntr    - print layer # on xsection every dcntr layer
#        zstart   - if dcntr>0, where to start layer notation
#        jlrs     - show location (where layers are printed out, e.g.)
#        sct_show - show section location on map inset
#
  print('Plotting section '+xsct)

  SCT = Jsections()
  aa = SCT.get(xsct)
  i1 = aa[0]
  i2 = aa[1]
  j1 = aa[2]
#
# For flat shading need to have X,Z dimension +1
# i.e. vertices of the grid cells
# cheat with horizontal dir - take the mid-point Longitudes
  rhox= np.squeeze(Rho[:,j1,i1:i2])
  Hb  = HH[j1,i1:i2+1]
  XX  = LON[j1,i1:i2+1]
  y0  = LAT[j1,i1]
  Bmin = np.floor(np.min(Hb))-100.

  dmm = ZZ0.shape
  if len(dmm)==1:
    ZZ = ZZ0
  else:
    ZZ = ZZ0[:,j1,i1:i2+1]

#
# Make XX and ZZ same dimensions:
# It is either same dim or ZZ(depths) > XX dim
  lXX = len(XX.shape)
  lZZ = len(ZZ.shape)

  if lZZ > lXX and lXX == 1:
    kzz = ZZ.shape[0]
    nzz = ZZ.shape[1]
    dmm = XX.copy()
    dmm = np.expand_dims(dmm, axis=0) 
    XX  = np.expand_dims(XX, axis=0)
    for kk in range(1,kzz):
      XX = np.append(XX, dmm, axis=0)



# Colormap
  if not clrmp: 
#  clrmp = copy(plt.cm.rainbow)
    clrmp = copy(plt.cm.gist_ncar_r)

  clrmp.set_bad(color=[0.7,0.7,0.7])

# For positive/negative make colorbar symmetric for comparison
  if (not rmax) or (not rmin):
    rmin = np.floor(np.nanmin(rhox*10))/10
    rmax = np.ceil(np.nanmax(rhox*10))/10
    if rmin < 0. and rmax >0:
      dmm = max([abs(rmin), rmax])
      rmin = -dmm
      rmax = dmm

  plt.ion()
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.25, 0.8, 0.7])
  im1 = ax1.pcolormesh(XX,ZZ,rhox, shading='flat', \
                 cmap=clrmp,\
                 vmin=rmin, \
                 vmax=rmax)

  #ax1.plot(XX,Hb)
  # Patch bottom:
  lXX = len(XX.shape)
  if lXX == 1:
    xx1 = XX
  else:
    xx1 = XX[0,:]

  verts = [(np.min(xx1),-8000),*zip(xx1,Hb),(np.max(xx1),-8000)]
  poly = Polygon(verts, facecolor='0.6', edgecolor='0.6')
  ax1.add_patch(poly)

  #ax1.autoscale(enable=True, axis='both', tight='True')
  ax1.set_xlim([np.min(xx1), np.max(xx1)])
  ax1.set_ylim([Bmin, 0])

  if stl:
    ax1.set_title(stl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='both')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.tick_params(direction='in', length=12)
  plt.sca(ax1)


# Draw a section on the map:
  xdm = HH.shape[1]
  ydm = HH.shape[0] 
  if sct_show:
    ax2 = plt.axes([0.75, 0.025, 0.18, 0.18])
    ax2.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
    ax2.plot((i1,i2),(j1,j1),'-',color=[1.,0.4,0])
    dii = 350
    ilim1 = max([i1-dii,0])
    ilim2 = min([i2+dii,xdm])
    jlim1 = max([j1-dii,0])
    jlim2 = min([j1+dii,ydm])

    ax2.set_xlim([ilim1,ilim2])
    ax2.set_ylim([jlim1,jlim2])

  if btx:
    bottom_text(btx,pos=[0.02, 0.03])

  return Hb, XX


def plot_A3d_Jsct(xsct,Rho,ZZ0,HH,LON,LAT, stl=[], btx=[], sct_show=[], \
                  fgnmb=5, rmin=[], rmax=[], clrmp=[], dcntr=1, zstart=-1000., jlrs=-1):
#
# Plot any section from 3D array - similar to _rho_Jsct but more options
#
# Write contours at the deepest point with dcntr intervals
# make dcntr=0 if not needed
# jrls >=0 indicate location where layer thickness is printed out
#
# Sections: i1,i2,j1 - W-E sections
#
# Input: xsct     - section name from dictionary SCT
#        Rho      - 3D field to plot (density)
#        ZZ0      - z-level interface depths 1D array or 3D
#                   if 1D constant Z is assumed for the section
#                   3D - Z layers extracted for the section
#        HH       - bathymetry 2D
#        LON/LAT  - model grid 2D
#        fgnmb    - figure # to plot
#        dcntr    - print layer # on xsection every dcntr layer
#        zstart   - if dcntr>0, where to start layer notation
#        jlrs     - show location (where layers are printed out, e.g.)
#        sct_show - show section location on map inset
#
  print('Plotting section '+xsct)

  SCT = Jsections()
  aa = SCT.get(xsct)
  i1 = aa[0]
  i2 = aa[1]
  j1 = aa[2]
#
# For flat shading need to have X,Z dimension +1
# i.e. vertices of the grid cells
# cheat with horizontal dir - take the mid-point Longitudes
  rhox= np.squeeze(Rho[:,j1,i1:i2])
  Hb  = HH[j1,i1:i2+1]
  XX  = LON[j1,i1:i2+1]
  y0  = LAT[j1,i1]
  Bmin = np.floor(np.min(Hb))-100.
  X1  = XX.copy()

  dmm = ZZ0.shape
  if len(dmm)==1:
    ZZ = ZZ0
  else:
    ZZ = ZZ0[:,j1,i1:i2+1]

#
# Make XX and ZZ same dimensions:
# It is either same dim or ZZ(depths) > XX dim
  lXX = len(XX.shape)
  lZZ = len(ZZ.shape)

  if lZZ > lXX and lXX == 1:
    kzz = ZZ.shape[0]
    nzz = ZZ.shape[1]
    dmm = XX.copy()
    dmm = np.expand_dims(dmm, axis=0) 
    XX  = np.expand_dims(XX, axis=0)
    for kk in range(1,kzz):
      XX = np.append(XX, dmm, axis=0)



# Colormap
  if not clrmp: 
#  clrmp = copy(plt.cm.rainbow)
    clrmp = copy(plt.cm.gist_ncar_r)

  clrmp.set_bad(color=[0.7,0.7,0.7])

# For positive/negative make colorbar symmetric for comparison
  if (not rmax) or (not rmin):
    rmin = np.floor(np.nanmin(rhox*10))/10
    rmax = np.ceil(np.nanmax(rhox*10))/10
    if rmin < 0. and rmax >0:
      dmm = max([abs(rmin), rmax])
      rmin = -dmm
      rmax = dmm

  plt.ion()
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.25, 0.8, 0.7])
  im1 = ax1.pcolormesh(XX,ZZ,rhox, shading='flat', \
                 cmap=clrmp,\
                 vmin=rmin, \
                 vmax=rmax)

  #ax1.plot(XX,Hb)
  # Patch bottom:
  lXX = len(XX.shape)
  if lXX == 1:
    xx1 = XX
  else:
    xx1 = XX[0,:]

  verts = [(np.min(xx1),-8000),*zip(xx1,Hb),(np.max(xx1),-8000)]
  poly = Polygon(verts, facecolor='0.6', edgecolor='0.6')
  ax1.add_patch(poly)

  iB = np.where(Hb == np.min(Hb))[0][0]
  if dcntr>0:
    cc=-1
    for kk in range(kzz):
      zlr = ZZ[kk,iB]
      if kk > 0:  
        dzz = abs(ZZ[kk,iB]-ZZ[kk-1,iB])
      else:
        dzz = 0.

      if zlr <= zstart:
        cc += 1
        if cc%dcntr == 0 and dzz > 0.1:
          zm = zlr+0.5*dzz
          txt = '{0}'.format(kk)
          ax1.text(X1[iB],zm,txt,fontsize=10)

# Inidcate where layers are printed out:
  if jlrs >= 0:
    xlrs = X1[jlrs]
    ax1.plot((xlrs,xlrs),(Bmin,0),'--',color=[0.8,0,0])

  #ax1.autoscale(enable=True, axis='both', tight='True')
  ax1.set_xlim([np.min(xx1), np.max(xx1)])
  ax1.set_ylim([Bmin, 0])

  if stl:
    ax1.set_title(stl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='both')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.tick_params(direction='in', length=12)
  plt.sca(ax1)


# Draw a section on the map:
  xdm = HH.shape[1]
  ydm = HH.shape[0] 
  if sct_show:
    ax2 = plt.axes([0.75, 0.025, 0.18, 0.18])
    ax2.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
    ax2.plot((i1,i2),(j1,j1),'-',color=[1.,0.4,0])
    dii = 350
    ilim1 = max([i1-dii,0])
    ilim2 = min([i2+dii,xdm])
    jlim1 = max([j1-dii,0])
    jlim2 = min([j1+dii,ydm])

    ax2.set_xlim([ilim1,ilim2])
    ax2.set_ylim([jlim1,jlim2])

  if btx:
    bottom_text(btx,pos=[0.02, 0.03])

  return Hb, X1, ZZ, rhox



def print_1D(A,wd=8,prc=2):
  """
    Print out 1 colume of data
  """
  ndim1 = A.shape[0]
  for k in range (ndim1):
    print('{0}: {1:{width}.{precis}f}'.format(k+1,A[k],width=wd,precis=prc))

  return

def print_2D(A1,A2,wd=8,prc=2,kend=[]):
  """
    Print out 2 columns of data, same size
  """
  ndim1 = A1.shape[0]
  if not kend:
    kend = ndim1

  for k in range (kend):
    print('{0}: {1:{width}.{precis}f}  {2:{width}.{precis}f}'.\
          format(k+1,A1[k],A2[k],width=wd,precis=prc))

  return

def clrmp_BlGrWhYlOrRd(nclrs=200):
  """
    Create colormap blue - Green - White - Yellow - Orange - Red
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  nhlf = round(nclrs/2)
  top = cm.get_cmap('YlOrRd',nhlf)
  btm = cm.get_cmap('GnBu_r',nhlf)
#
# Make white transition between half-colormaps
  clrtop = top(range(nhlf))
  clrbtm = btm(range(nhlf))

  ixtop  = round(nhlf*0.1)-1
  cxtop  = clrtop[ixtop,:]
  ixbtm  = nhlf-ixtop-1
  cxbtm  = clrbtm[ixbtm,:]
  chtop  = np.zeros((ixtop,4))
  chbtm  = np.zeros((ixtop,4))

  chtop[:,3] = cxtop[3]  # alfa
  chbtm[:,3] = cxbtm[3]

  for ik in range(3):
    chtop[:,ik]  = np.linspace(cxtop[ik],1,ixtop)

  for ik in range(3):
    chbtm[:,ik]  = np.linspace(cxbtm[ik],1,ixtop)

  chtop = np.flip(chtop,0)

  clrtop[0:ixtop,:]  = chtop
  clrbtm[ixbtm+1:nhlf,:] = chbtm

  newclrs = np.concatenate((clrbtm,clrtop))
  newcmp  = ListedColormap(newclrs)

  return newcmp


