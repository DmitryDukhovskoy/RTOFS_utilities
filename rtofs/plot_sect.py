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

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
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
    "SPac1"  : [1512,2580,1000],
    "SPac2"  : [2000,2600,745],
    "GoM1"   : [2340,2560,1798],
    "NAtl1"  : [2900,3000,1897],
    "SeaJpn" : [650,820,2049],
    "Test1"  : [3000,3250,819],
    "Test2"  : [1300,1550,1324]
  }

  return SCT

def find_indx_lonlat(x0,y0,LON,LAT,xsct="none"):
  """
  Find closest grid point to lon/lat coordinate
  """
  if x0 > 180.:
    x0 = x0-360.

  dmm = np.sqrt((LON-x0)**2+(LAT-y0)**2)
  jj0, ii0 = np.where(dmm == np.min(dmm)) # global indices
#
# If section provided, indices for this section:
  ip0 = -1
  jp0 = -1

  if not xsct=="none":
    SCT = Jsections()
    aa  = SCT.get(xsct)
    i1  = aa[0]
    i2  = aa[1]
    j1  = aa[2]
    ip0 = ii0[0]-i1
    jp0 = jj0[0]-j1 

  return ii0[0], jj0[0], ip0, jp0


def plot_Jsection(xsct,ZZ,HH,LON,LAT, drfnm, sct_show=False, \
                  fgnmb=5,dcntr=1, zstart=-1000., jlrs=-1, \
                  btx=[]):
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

# Check is lon is at singular longitude:
# Should happen only when X(1, ...) >0 and after 180E, X<0 (-179,-178,...)
  if XX[0] > XX[-1]:
    dmm = XX-360.
    XX = np.where(XX >= 0.,dmm,XX)

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

# indices in Fortran convention:
  stl = ('{5}, {4}, Hybrd layers, lat={0:5.2f}, jF={1}, iF={2}:{3}'.\
         format(y0,j1+1,i1+1,i2+1,drfnm,xsct))
  ax1.set_title(stl)

#  btx = 'find_collapsed_lrs.py'
#  bottom_text(btx)

#  if jp0 >= 0:
#    print_1D(dZZs[:,jp0])

#
# Draw a section on the map:
#  if sct_show:
#    sct_show.plot((i1,i2),(j1,j1),'-',color=[0.,0.4,1])
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

  return ZZs, Hb, XX, dZZs


def plot_rho_Jsct(xsct,Rho,ZZ0,HH,LON,LAT, stl=[], btx=[], sct_show=[], \
                  fgnmb=5, rmin=[], rmax=[], clrmp=[], jlrs=-1):
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
#        jlrs     - show location (where field is printed out by layers)
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

# Check is lon is at singular longitude:
# Should happen only when X(1, ...) >0 and after 180E, X<0 (-179,-178,...)
  if XX[0] > XX[-1]:
    dmm = XX-360.
    XX = np.where(XX >= 0.,dmm,XX)

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

# Inidcate where info by layers printed out:
  if jlrs >= 0:
    xlrs = XX[0,jlrs]
    ax1.plot((xlrs,xlrs),(Bmin,0),'--',color=[0.8,0,0])
    iprf = i1+jlrs-1
    jprf = j1 
    ax1.text(xlrs,0.9*Bmin,'  i={0}\n  j={1}'.format(iprf,jprf))


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

  return Hb, XX, rhox, ZZ 

#
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


def plot_A3d_Jsct(xsct,Rho,ZZ0,HH,LON,LAT, stl=[], btx=[], sct_show=[], \
                  fgnmb=5, rmin=[], rmax=[], clrmp=[], dcntr=1, \
                  zstart=-1000., jlrs=-1):
#
# Plot any section from 3D array - similar to _rho_Jsct but more options
#
# Write contours at the deepest point with dcntr intervals
# make dcntr=0 if not needed
# jrls >=0 indicate location where A3d is printed out by layers
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


def plot_2Dmap(A2d,HH,fgnmb=1,clrmp=[],rmin=0.,rmax=50.,\
               ctl="2D map",Ipnt=[],Jpnt=[],btx=" "):
  """
  Plot 2 D map
  Ipnt,Jpnt - locations to plot on the map
  """
  print('Plotting  '+ctl+'...')
  if not clrmp:
    clrmp = copy(plt.cm.afmhot_r)

  clrmp.set_bad(color=[0.7,0.7,0.7])
  plt.ion()
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im1 = plt.pcolormesh(A2d,shading='flat',\
                       vmin=rmin, vmax=rmax, cmap=clrmp)
  plt.contour(HH,[0.0],colors=[(0.5,0.5,0.5)],linewidths=1)

  if len(Ipnt)>0:
    plt.plot(Ipnt,Jpnt,linestyle='',marker=".",markersize=5,color="m")

  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  if rmin == 0:
    clb = plt.colorbar(im1, cax=ax2, extend='max')
  elif rmax == 0:
    clb = plt.colorbar(im1, cax=ax2, extend='min')
  else:
    clb = plt.colorbar(im1, cax=ax2, extend='both')

  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(ticklabs,fontsize=12)

  plt.sca(ax1)
  ax1.set_title(ctl)

  bottom_text(btx)

# Print several locations of points:
  if len(Ipnt):
    ax3 = plt.axes([0.1,0.3,0.25,0.4])
    nmax = min(20,Ipnt.shape[0])
    ss1  = 'Some locations > eps \n  I      J\n'
    for kk in range(nmax):
      ss1 = ss1+'{0}  {1}\n'.format(Ipnt[kk],Jpnt[kk])
    ax3.text(0.1,0.1,ss1)
    ax3.axis('off')

  return

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

def print_3D(A1,A2,A3,wd=8,prc=2,kend=[]):
  """
    Print out 3 columns of data, same size
  """
  ndim1 = A1.shape[0]
  if not kend:
    kend = ndim1

  for k in range (kend):
    print('{0:3d}: {1:{width}.{precis}f} {2:{width}.{precis}f} {3:{width}.{precis}f}'.\
          format(k+1,A1[k],A2[k],A3[k],width=wd,precis=prc))

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

def clrmp_WhYlOrRd(nclrs=200):
  """
    Create colormap White - Yellow - Orange - Red
    showing positive values e.g. from 0 to ...
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  nhlf = nclrs
  top = cm.get_cmap('YlOrRd',nhlf)
#
# Make white transition between half-colormaps
  clrtop = top(range(nhlf))

  ixtop  = round(nhlf*0.1)-1
  cxtop  = clrtop[ixtop,:]
  chtop  = np.zeros((ixtop,4))

  chtop[:,3] = cxtop[3]  # alfa

  for ik in range(3):
    chtop[:,ik]  = np.linspace(cxtop[ik],1,ixtop)

  chtop              = np.flip(chtop,0)
  clrtop[0:ixtop,:]  = chtop

  newclrs = clrtop
  newcmp  = ListedColormap(newclrs)

  return newcmp


def clrmp_BlGrWh(nclrs=200):
  """
    Create colormap blue - Green - White
    for negatives
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  nhlf = nclrs
  btm = cm.get_cmap('GnBu_r',nhlf)
#
# Make white transition between half-colormaps
  clrbtm = btm(range(nhlf))

  ixtop  = round(nhlf*0.1)-1
  ixbtm  = nhlf-ixtop-1
  cxbtm  = clrbtm[ixbtm,:]
  chbtm  = np.zeros((ixtop,4))

  chbtm[:,3] = cxbtm[3]

  for ik in range(3):
    chbtm[:,ik]  = np.linspace(cxbtm[ik],1,ixtop)

  clrbtm[ixbtm+1:nhlf,:] = chbtm

  newclrs = clrbtm
  newcmp  = ListedColormap(newclrs)

  return newcmp

def plot_1prof(B,Z,ctl='Profile'):
  plt.ion()
  print('Plotting '+ctl)
  fig1 = plt.figure(1,figsize=(7,8),constrained_layout=False)
  plt.clf()
  plt.plot(B,Z)
  plt.title(ctl)

  return

def plot_1prof_map(B,Z,HH,zlim=-1.e6,xlim1=[],xlim2=[],\
                   ctl='Porfile',fgn=1,i1=[],j1=[]):
  """
  Plot profiles and show location on a map
  limit to the upper zlim meters
  """
  from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

  plt.ion()
  print('Plotting '+ctl)
  fig1 = plt.figure(fgn,figsize=(7,8), constrained_layout=False)
  fig1.clf()
  ax1 = plt.axes([0.1, 0.1, 0.35, 0.8])
  plt.plot(B,Z)
  
  klim = max(np.where(Z>=zlim)[0])
  zmin = max([zlim,np.min(Z)])
  if not xlim1 or not xlim2:
    xlim1 = np.min(B[0:klim+1])
    xlim2 = np.max(B[0:klim+1])
    dx = abs(xlim2-xlim1)
    if dx > 1:
      xlim1 = np.floor(xlim1)
      xlim2 = np.ceil(xlim2)
      ax1.xaxis.set_major_locator(MultipleLocator(4))
      ax1.xaxis.set_minor_locator(MultipleLocator(2))
    else:
      xlim1 = 0.1*(np.floor(xlim1*10))
      xlim2 = 0.1*(np.ceil(xlim2*10))
      ax1.xaxis.set_major_locator(MultipleLocator(0.2))
      ax1.xaxis.set_minor_locator(MultipleLocator(0.1))

  plt.ylim(ymax = 0, ymin = zmin) 
  plt.xlim(xmax = xlim2, xmin = xlim1)

  ax1.grid(True)
#  ax1.set_xlabel('cycle/hr')
  ax1.set_title(ctl)

  xdm = HH.shape[1]
  ydm = HH.shape[0]

  ax2 = plt.axes([0.55, 0.5, 0.35, 0.35])
  ax2.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
  if i1:
    ax2.plot((i1),(j1),marker='.',markersize=8,color=[1.,0.4,0])
    dii = 400
    ilim1 = max([i1-dii,0])
    ilim2 = min([i1+dii,xdm])
    jlim1 = max([j1-dii,0])
    jlim2 = min([j1+dii,ydm])

    ax2.set_xlim([ilim1,ilim2])
    ax2.set_ylim([jlim1,jlim2])

  ax2.grid(True)

  return


def plot_Nprof_map(A3d,ZM,i0,j0,HH,dI=1,zlim=-1.e6,xlim1=[],xlim2=[],\
                   ctl='Porfile',fgn=1):
  """
  Plot N profiles around i0,j0
  and show location on a map
  limit to the upper zlim meters
  """
  from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
  plt.ion()
  print('Plotting '+ctl)
  fig1 = plt.figure(fgn,figsize=(7,8), constrained_layout=False)
  fig1.clf()
  ax1 = plt.axes([0.1, 0.1, 0.35, 0.8])
  xdm = HH.shape[1]
  ydm = HH.shape[0]

  ii1 = i0-dI
  ii1 = max([0,ii1])
  ii2 = i0+dI
  ii2 = min([xdm,ii2])
  jj1 = j0-dI
  jj1 = max([0,jj1])
  jj2 = j0+dI
  jj2 = min([ydm,jj2])

  bmax = -1.e6
  bmin = 1.e6
  for ii in range(ii1,ii2+1):
    for jj in range(jj1,jj2+1):
      if HH[jj,ii] >= 0.:
        continue
      B = A3d[:,jj,ii]
      Z = ZM[:,jj,ii]
      plt.plot(B,Z)
      klim = max(np.where(Z>=zlim)[0])
      bmax = max([bmax,np.max(B[0:klim+1])])
      bmin = min([bmin,np.min(B[0:klim+1])])
  
  zmin = max([zlim,np.min(Z)])
  if not xlim1 or not xlim2:
    xlim1 = bmin
    xlim2 = bmax 
    dx = abs(xlim2-xlim1)
    if dx >= 3.:
      xlim1 = np.floor(xlim1-dx/10.)
      xlim2 = np.ceil(xlim2+dx/10.)
#      plt.xticks(np.arange(xlim1-1,xlim2+1, step=2.))
      ax1.xaxis.set_major_locator(MultipleLocator(4))
      ax1.xaxis.set_minor_locator(MultipleLocator(2))
    else:
      mjtck = 0.2
      mntck = 0.1
      xlim1 = 0.1*(np.floor(xlim1*10-dx))
      xlim2 = 0.1*(np.ceil(xlim2*10+dx))
      ax1.xaxis.set_major_locator(MultipleLocator(mjtck))
      ax1.xaxis.set_minor_locator(MultipleLocator(mntck))


  plt.ylim(ymax = 0, ymin = zmin) 
  plt.xlim(xmax = xlim2, xmin = xlim1)

  ax1.grid(True)
#  ax1.set_xlabel('cycle/hr')
  ax1.set_title(ctl)


  ax2 = plt.axes([0.55, 0.5, 0.35, 0.35])
  ax2.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
  for ii in range(ii1,ii2+1):
    for jj in range(jj1,jj2+1):
      ax2.plot((ii),(jj),marker='.',markersize=8,color=[1.,0.4,0])

  dii = 400
  ilim1 = max([ii1-dii,0])
  ilim2 = min([ii2+dii,xdm])
  jlim1 = max([jj1-dii,0])
  jlim2 = min([jj2+dii,ydm])

  ax2.set_xlim([ilim1,ilim2])
  ax2.set_ylim([jlim1,jlim2])

  ax2.grid(True)

  return

def plot_2prof(Z1,T1,Z2,T2,ctl1="",ctl2="",fgn=1):
  plt.ion()
  print('Plotting '+ctl1)
  fig1 = plt.figure(fgn,figsize=(7,8), constrained_layout=False)
  fig1.clf()
  ax1 = plt.axes([0.1, 0.1, 0.35, 0.8])
  plt.plot(T1,Z1)
  ax1.grid(True)
#  ax1.set_xlabel('cycle/hr')
  ax1.set_title(ctl1)

  print('Plotting '+ctl2)
  ax2 = plt.axes([0.55, 0.1, 0.35, 0.8])
  ax2.plot(T2,Z2)
  ax2.grid(True)
  ax2.set_title(ctl2)

  return


def plot_TSprof_stat(Zw,tw_md,tw_lp,tw_up, Zg,tg_md,tg_lp,tg_up, ctl1="",fgn=1,\
                     lb1='WOA18',lb2='GDEM4'):
  """
    Plot T or S profile statistics for WOA and GDEM on the same plot
    tw_md - median or mean profile
    tw_lp/up - lower/upper bounds of the profile (e.g., 10, 90th percentiles)
  """
  plt.ion()
  fig1 = plt.figure(fgn,figsize=(7,8), constrained_layout=False)
  fig1.clf()
  ax1 = plt.axes([0.1, 0.1, 0.5, 0.8])
  clr1=[0.,0.4,0.8]
  plt.plot(tw_md,Zw,color=clr1,linewidth=2.5,label=lb1)
  plt.plot(tw_lp,Zw,color=clr1,linewidth=1)
  plt.plot(tw_up,Zw,color=clr1,linewidth=1)

  clr2 = [0.9,0.3,0]
  plt.plot(tg_md,Zg,color=clr2,linewidth=2.1,label=lb2)
  plt.plot(tg_lp,Zg,color=clr2,linewidth=1)
  plt.plot(tg_up,Zg,color=clr2,linewidth=1)

  izb = np.where(~np.isnan(tw_md))[0][-1]
  hbtm = Zw[izb]-200.

  plt.yticks(np.arange(-10000,0,500))
  plt.ylim(hbtm,0)

  ax1.grid(True)
#  ax1.set_xlabel('cycle/hr')
  ax1.set_title(ctl1)

  legend = ax1.legend(bbox_to_anchor=(1.05,1), loc='upper left')

  return


def plot_EWsection(ix1, ix2, ix0, jx0, ZZ, HH, LON, LAT, \
                  fgnmb=5, dcntr=1, zstart=-1000., lrstart=0,\
                  stl='HYCOM layers', sinfo=[], btx=[], f_intract=True):
# Plot EW-section
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
#        lrstart  - do not plot layers <lrstart
#        jlrs     - show location (where layers are printed out, e.g.)
#        sct_show - show section location on map in figure(sct_show)
#        sinfo    - information text to put outside the plot
#        f_intract -  False - do not show figures, True - show figures
#
  i1 = ix1
  i2 = ix2
  j1 = jx0
  kdm = ZZ.shape[0]
  jdm = ZZ.shape[1]
  idm = ZZ.shape[2]

# i1 < i2
# Check if section is at the boundary
# assuming that desired section is the shortest
  f_xob = False
  if (i1 > i2) and ((i1-i2) < (idm-i1)+i2):  # simply swap i1 and i2
    imm=i2
    i2=i1
    i1=imm
    ZZs  = ZZ[:,j1,i1:i2]
    Hb   = HH[j1,i1:i2]
    XX   = LON[j1,i1:i2]
  elif i1 > i2:      
# Section across the OB:
    f_xob = True
    Z1  = ZZ[:,j1,i1:]
    Hb1 = HH[j1,i1:]
    X1  = LON[j1,i1:]
    Z2  = ZZ[:,j1,:i2]
    Hb2 = HH[j1,:i2]
    X2  = LON[j1,:i2]
    HH1 = HH[:,i1:]
    HH2 = HH[:,:i2]
    Hsub= np.append(HH1,HH2, axis=1)
    ZZs = np.append(Z1,Z2, axis=1)
    Hb  = np.append(Hb1,Hb2)
    XX  = np.append(X1,X2)
  else:
    ZZs  = ZZ[:,j1,i1:i2]
    Hb   = HH[j1,i1:i2]
    XX   = LON[j1,i1:i2]

  lat0 = LAT[jx0,ix0]
  lon0 = LON[jx0,ix0]
  Bmin = np.floor(np.min(Hb))-100.
  dZZs = np.abs(np.diff(ZZs, axis=0))

  ll = ZZ.shape[0]-1 

# Check is lon is at singular longitude:
# Should happen only when X(1, ...) >0 and after 180E, X<0 (-179,-178,...)
  if XX[0] > XX[-1]:
    dmm = XX-360.
    XX = np.where(XX >= 0.,dmm,XX)

  if f_intract:
    plt.ion()
  else:
    plt.ioff()

  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.25, 0.8, 0.7])
  #ax1.plot(XX,Hb)
  # Patch bottom:
  verts = [(np.min(XX),-8000),*zip(XX,Hb),(np.max(XX),-8000)]
  poly = Polygon(verts, facecolor='0.6', edgecolor='0.6')
  ax1.add_patch(poly)

  CCLR = np.zeros((9,3))
  CCLR[0,:] = [0.5,0,0.5]
  CCLR[1,:] = [0.8,0,0.9]
  CCLR[2,:] = [1.0,0,0.5]
  CCLR[3,:] = [1.,0.,0]
  CCLR[4,:] = [0.,0.,1]
  CCLR[5,:] = [0.,0.5,1]
  CCLR[6,:] = [0.,0.7,0.3]
  CCLR[7,:] = [0.,1,0.2]
  CCLR[8,:] = [0.8,0.5,0.]

  # Plot interfaces
  iclr = 0
  for kk in range(ll+1):
    dmm = ZZs[kk,:]
    if kk > 0 and kk < lrstart:
      continue
#    if kk > 0:
#      dzz = dZZs[kk-1,:]
#      if all (dzz == 0):  # all layers at the bottom
#        break
    if (kk+1)%5==0:
      cclr = CCLR[iclr,:]
      iclr = iclr+1
    else:
      cclr=[0,0,0]

    ax1.plot(XX,dmm,color=cclr)

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

  #ax1.autoscale(enable=True, axis='both', tight='True')
  ax1.set_xlim([np.min(XX), np.max(XX)])
  ax1.set_ylim([Bmin, 0])

  ax1.set_title(stl)


# Plot map:
  import mod_utils as mutil
  if not f_xob:
    LMSK = HH.copy()
    LMSK = np.where(LMSK<0.,0.,1.)
    lcmp = mutil.clrmp_lmask(2, clr_land=[0.5,0.5,0.5])
    ax3 = plt.axes([0.75, 0.02, 0.2, 0.2])
    ax3.pcolormesh(LMSK, shading='flat',\
                    cmap=lcmp, \
                    vmin=0, \
                    vmax=1)
    ax3.contour(HH, [-5000,-4000,-3000,-2000,-1000], \
                colors=[(0.9,0.9,0.9)], linestyles='solid',linewidths=1)
    ax3.plot([i1,i2], [j1,j1], '-', color=[1.,0.2,0])

    dmin  = 400
    dh    = int(max([2*dmin-abs(ix2-ix1),60])*0.5)
    dii   = max([abs(ix2-ix1),dmin])
    ilim1 = max([ix1-dh,0])
    ilim2 = min([ix2+dh,idm])
    jlim1 = max([jx0-dii,0])
    jlim2 = min([jx0+dii,jdm])

    ax3.axis('scaled')
  #  ax3.axis('square')
    ax3.set_xlim([ilim1,ilim2])
    ax3.set_ylim([jlim1,jlim2])
  else:    # domain cut by the OB
    LMSK = Hsub.copy()
    LMSK = np.where(LMSK<0.,0.,1.)
    lcmp = mutil.clrmp_lmask(2, clr_land=[0.5,0.5,0.5])
    ax3 = plt.axes([0.75, 0.02, 0.2, 0.2])
    ax3.pcolormesh(LMSK, shading='flat',\
                    cmap=lcmp, \
                    vmin=0, \
                    vmax=1)
    ax3.contour(Hsub, [-5000,-4000,-3000,-2000,-1000], \
                colors=[(0.9,0.9,0.9)], linestyles='solid',linewidths=1)
    i1s = 0
    i2s = Hsub.shape[1]
    ax3.plot([i1s,i2s], [j1,j1], '-', color=[1.,0.2,0])

    dmin  = 400
    dh    = int(max([2*dmin-abs(i2s-i1s),60])*0.5)
    dii   = min([abs(i2s-i1s),dmin])
    ilim1 = max([i1s-dh,0])
    ilim2 = min([i2s+dh,idm])
    jlim1 = max([jx0-dii,0])
    jlim2 = min([jx0+dii,jdm])

    ax3.axis('scaled')
  #  ax3.axis('square')
    ax3.set_xlim([ilim1,ilim2])
    ax3.set_ylim([jlim1,jlim2])


# Set y, x axis not visible
#ahd = plt.gca()
  xax = ax3.axes.get_xaxis()
  xax = xax.set_visible(False)
  yax = ax3.axes.get_yaxis()
  yax = yax.set_visible(False)

# Info text:
  if sinfo:
    ax4 = plt.axes([0.1,0.08,0.5,0.15])
    ax4.text(0.0,0.1,sinfo)
    ax4.axis('off') 

  if btx:
    bottom_text(btx,pos=[0.02, 0.03])

  return ZZs, Hb, XX, dZZs

def plot_SNsection(jx1, jx2, ix0, jx0, ZZ, HH, LON, LAT, \
                  fgnmb=6, dcntr=1, zstart=-1000., lrstart=0,\
                  stl='HYCOM layers', sinfo=[], btx=[], f_intract=True):
# Plot SN-section
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
#        lrstart  - do not plot layers <lrstart
#        jlrs     - show location (where layers are printed out, e.g.)
#        sct_show - show section location on map in figure(sct_show)
#
  i1 = ix0
  j1 = jx1
  j2 = jx2
  kdm = ZZ.shape[0]
  jdm = ZZ.shape[1]
  idm = ZZ.shape[2]


  ZZs  = ZZ[:,j1:j2,i1]
  Hb   = HH[j1:j2,i1]
  XX   = LAT[j1:j2,i1]
  lat0 = LAT[jx0,ix0]
  lon0 = LON[jx0,ix0]
  Bmin = np.floor(np.min(Hb))-100.
  dZZs = np.abs(np.diff(ZZs, axis=0))

  ll = ZZ.shape[0]-1 

# Check is lon is at singular longitude:
# Should happen only when X(1, ...) >0 and after 180E, X<0 (-179,-178,...)
  if XX[0] > XX[-1]:
    dmm = XX-360.
    XX = np.where(XX >= 0.,dmm,XX)

  if f_intract:
    plt.ion()
  else:
    plt.ioff()

  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.25, 0.8, 0.7])
  #ax1.plot(XX,Hb)
  # Patch bottom:
  verts = [(np.min(XX),-8000),*zip(XX,Hb),(np.max(XX),-8000)]
  poly = Polygon(verts, facecolor='0.6', edgecolor='0.6')
  ax1.add_patch(poly)

  CCLR = np.zeros((9,3))
  CCLR[0,:] = [0.5,0,0.5]
  CCLR[1,:] = [0.8,0,0.9]
  CCLR[2,:] = [1.0,0,0.5]
  CCLR[3,:] = [1.,0.,0]
  CCLR[4,:] = [0.,0.,1]
  CCLR[5,:] = [0.,0.5,1]
  CCLR[6,:] = [0.,0.7,0.3]
  CCLR[7,:] = [0.,1,0.2]
  CCLR[8,:] = [0.8,0.5,0.]

  iclr = 0
  # Plot interfaces
  for kk in range(ll+1):
    dmm = ZZs[kk,:]
    if kk > 0 and kk < lrstart:
      continue
#    if kk > 0:
#      dzz = dZZs[kk-1,:]
#      if all (dzz == 0):  # all layers at the bottom
#        break
    if (kk+1)%5==0:
#      cclr=[0.8,0,0]
      cclr = CCLR[iclr,:]
      iclr = iclr+1
      if iclr > len(CCLR):
        iclr=0
    else:         
      cclr=[0,0,0]
    ax1.plot(XX,dmm,color=cclr)

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

  #ax1.autoscale(enable=True, axis='both', tight='True')
  ax1.set_xlim([np.min(XX), np.max(XX)])
  ax1.set_ylim([Bmin, 0])

  ax1.set_title(stl)


# Plot map:
  import mod_utils as mutil
  LMSK = HH.copy()
  LMSK = np.where(LMSK<0.,0.,1.)
  lcmp = mutil.clrmp_lmask(2, clr_land=[0.5,0.5,0.5])
  ax3 = plt.axes([0.75, 0.02, 0.2, 0.2])
  ax3.pcolormesh(LMSK, shading='flat',\
                  cmap=lcmp, \
                  vmin=0, \
                  vmax=1)
  ax3.contour(HH, [-5000,-4000,-3000,-2000,-1000], \
              colors=[(0.9,0.9,0.9)], linestyles='solid',linewidths=1)
  ax3.plot([i1,i1], [j1,j2], '-', color=[1.,0.2,0])

  dmin  = 400
  dh    = int(max([2*dmin-abs(jx2-jx1),60])*0.5)
  dii   = max([abs(jx2-jx1),dmin])
  ilim1 = max([ix0-dii,0])
  ilim2 = min([ix0+dii,idm])
  jlim1 = max([jx1-dh,0])
  jlim2 = min([jx2+dh,jdm])

  ax3.axis('scaled')
  ax3.set_xlim(ilim1,ilim2)
  ax3.set_ylim(jlim1,jlim2)

# Set y, x axis not visible
#ahd = plt.gca()
  xax = ax3.axes.get_xaxis()
  xax = xax.set_visible(False)
  yax = ax3.axes.get_yaxis()
  yax = yax.set_visible(False)

# Info text:
  if sinfo:
    ax4 = plt.axes([0.1,0.08,0.5,0.15])
    ax4.text(0.0,0.1,sinfo)
    ax4.axis('off') 

  if btx:
    bottom_text(btx,pos=[0.02, 0.03])

  return ZZs, Hb, XX, dZZs


# ==========
def plot_EWsectTS(ix1, ix2, ix0, jx0, ZZ, A3d, HH, LON, LAT, \
                 cmpr, rmin=[], rmax=[],\
                 fgnmb=5, dcntr=1, zstart=-1000., lrstart=999,\
                 stl='HYCOM TS', sinfo=[], btx=[], f_intract=True):
# Plot EW-section Salinity
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
#        lrstart  - do not plot layers <lrstart
#        sinfo    - information text to put outside the plot
#
  i1 = ix1
  i2 = ix2
  j1 = jx0
  kdm = ZZ.shape[0]
  jdm = ZZ.shape[1]
  idm = ZZ.shape[2]

# i1 < i2
# Check if section is at the boundary
# assuming that desired section is the shortest
  f_xob = False
  if (i1 > i2) and ((i1-i2) < (idm-i1)+i2):  # simply swap i1 and i2
    imm=i2
    i2=i1
    i1=imm
    ZZs  = ZZ[:,j1,i1:i2]
    Hb   = HH[j1,i1:i2]
    XX   = LON[j1,i1:i2]
    A2d  = A3d[j1,i1:i2]
  elif i1 > i2:      
# Section across the OB:
    f_xob = True
    Z1  = ZZ[:,j1,i1:]
    Hb1 = HH[j1,i1:]
    X1  = LON[j1,i1:]
    Z2  = ZZ[:,j1,:i2]
    Hb2 = HH[j1,:i2]
    X2  = LON[j1,:i2]
    HH1 = HH[:,i1:]
    HH2 = HH[:,:i2]
    Hsub= np.append(HH1,HH2, axis=1)
    ZZs = np.append(Z1,Z2, axis=1)
    Hb  = np.append(Hb1,Hb2)
    XX  = np.append(X1,X2)
    A1  = A3d[:,j1,i1:]
    A2  = A3d[:,j1,:i2]
    A2d = np.append(A1,A2, axis=1)
  else:
    ZZs  = ZZ[:,j1,i1:i2]
    Hb   = HH[j1,i1:i2]
    XX   = LON[j1,i1:i2]
    A2d  = A3d[:,j1,i1:i2]

  lat0 = LAT[jx0,ix0]
  lon0 = LON[jx0,ix0]
  Bmin = np.floor(np.min(Hb))-100.
  dZZs = np.abs(np.diff(ZZs, axis=0))

  ll = ZZ.shape[0]-1 

# Add layer for plotting:
  A2d = np.insert(A2d,0,A2d[0,:],axis=0)

# Check is lon is at singular longitude:
# Should happen only when X(1, ...) >0 and after 180E, X<0 (-179,-178,...)
  if XX[0] > XX[-1]:
    dmm = XX-360.
    XX = np.where(XX >= 0.,dmm,XX)

  if f_intract:
    plt.ion()
  else:
    plt.ioff()

  alf=2. 
  if not rmin:
    rmin = np.nanpercentile(A2d,alf)
    rmin = round(rmin,2)
    rmax = np.nanpercentile(A2d,100.-alf)
    rmax = round(rmax,2)
    print('rmin={0}, rmax={1}'.format(rmin,rmax))

  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.25, 0.8, 0.7])
  #ax1.plot(XX,Hb)
  im1 = ax1.pcolormesh(XX, ZZs, A2d, shading='flat',
                cmap=cmpr, vmin=rmin, vmax=rmax)
# im1.set_clim(34.5,35.5)
  # Patch bottom:
  verts = [(np.min(XX),-8000),*zip(XX,Hb),(np.max(XX),-8000)]
  poly = Polygon(verts, facecolor='0.1', edgecolor='0.1')
  ax1.add_patch(poly)


  # Plot interfaces
  for kk in range(ll+1):
    dmm = ZZs[kk,:]
    if kk < lrstart:
      continue

    cclr=[0,0,0]
#    ax1.plot(XX, dmm, color=cclr, linewidth=1)
    ax1.plot(XX, dmm, color=cclr, linewidth=0.5)


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
        if cc%dcntr == 0 and dzz > 1.0:
          zm = zlr+0.5*dzz
          txt = '{0}'.format(kk)
          ax1.text(XX[iB],zm,txt,fontsize=10)

  #ax1.autoscale(enable=True, axis='both', tight='True')
  ax1.set_xlim([np.min(XX), np.max(XX)])
  ax1.set_ylim([Bmin, 0])

  ax1.set_title(stl)

# Colorbar
#  ax2 = fig1.add_axes([ax1.get_position().x0, ax1.get_position().y0-0.1,
#                     ax1.get_position().width,0.02])
  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')
  

# Plot map:
  import mod_utils as mutil
  if not f_xob:
    LMSK = HH.copy()
    LMSK = np.where(LMSK<0.,0.,1.)
    lcmp = mutil.clrmp_lmask(2, clr_land=[0.5,0.5,0.5])
    ax3 = plt.axes([0.75, 0.02, 0.2, 0.2])
    ax3.pcolormesh(LMSK, shading='flat',\
                    cmap=lcmp, \
                    vmin=0, \
                    vmax=1)
    ax3.contour(HH, [-5000,-4000,-3000,-2000,-1000], \
                colors=[(0.9,0.9,0.9)], linestyles='solid',linewidths=1)
    ax3.plot([i1,i2], [j1,j1], '-', color=[1.,0.2,0])

    dmin  = 300
    dh    = int(max([2*dmin-abs(ix2-ix1),60])*0.5)
    dii   = max([abs(ix2-ix1),dmin])
    ilim1 = max([ix1-dh,0])
    ilim2 = min([ix2+dh,idm])
    jlim1 = max([jx0-dii,0])
    jlim2 = min([jx0+dii,jdm])

    ax3.axis('scaled')
  #  ax3.axis('square')
    ax3.set_xlim([ilim1,ilim2])
    ax3.set_ylim([jlim1,jlim2])
  else:    # domain cut by the OB
    LMSK = Hsub.copy()
    LMSK = np.where(LMSK<0.,0.,1.)
    lcmp = mutil.clrmp_lmask(2, clr_land=[0.5,0.5,0.5])
    ax3 = plt.axes([0.75, 0.02, 0.2, 0.2])
    ax3.pcolormesh(LMSK, shading='flat',\
                    cmap=lcmp, \
                    vmin=0, \
                    vmax=1)
    ax3.contour(Hsub, [-5000,-4000,-3000,-2000,-1000], \
                colors=[(0.9,0.9,0.9)], linestyles='solid',linewidths=1)
    i1s = 0
    i2s = Hsub.shape[1]
    ax3.plot([i1s,i2s], [j1,j1], '-', color=[1.,0.2,0])

    dmin  = 500
    dh    = int(max([2*dmin-abs(i2s-i1s),60])*0.5)
    dii   = min([abs(i2s-i1s),dmin])
    ilim1 = max([i1s-dh,0])
    ilim2 = min([i2s+dh,idm])
    jlim1 = max([jx0-dii,0])
    jlim2 = min([jx0+dii,jdm])

    ax3.axis('scaled')
  #  ax3.axis('square')
    ax3.set_xlim([ilim1,ilim2])
    ax3.set_ylim([jlim1,jlim2])


# Set y, x axis not visible
#ahd = plt.gca()
  xax = ax3.axes.get_xaxis()
  xax = xax.set_visible(False)
  yax = ax3.axes.get_yaxis()
  yax = yax.set_visible(False)

# Info text:
  if sinfo:
    ax4 = plt.axes([0.1,0.08,0.5,0.15])
    ax4.text(0.0,0.1,sinfo)
    ax4.axis('off') 

  if btx:
    bottom_text(btx,pos=[0.02, 0.03])

  return ZZs, Hb, XX, dZZs


def plot_SNsectTS(jx1, jx2, ix0, jx0, ZZ, A3d, HH, LON, LAT, \
                 cmpr, rmin=[], rmax=[],\
                 fgnmb=5, dcntr=1, zstart=-1000., lrstart=999,\
                 stl='HYCOM TS', sinfo=[], btx=[], f_intract=True):
# Plot SN-section Salinity
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
#        lrstart  - do not plot layers <lrstart
#        sinfo    - information text to put outside the plot
#
  i1 = ix0
  j1 = jx1
  j2 = jx2
  kdm = ZZ.shape[0]
  jdm = ZZ.shape[1]
  idm = ZZ.shape[2]

# i1 < i2
# Check if section is at the boundary
# assuming that desired section is the shortest
  ZZs  = ZZ[:,j1:j2,i1]
  Hb   = HH[j1:j2,i1]
  XX   = LAT[j1:j2,i1]
  lat0 = LAT[jx0,ix0]
  lon0 = LON[jx0,ix0]
  A2d  = A3d[:,j1:j2,i1]
  Bmin = np.floor(np.min(Hb))-100.
  dZZs = np.abs(np.diff(ZZs, axis=0))

  lat0 = LAT[jx0,ix0]
  lon0 = LON[jx0,ix0]
  Bmin = np.floor(np.min(Hb))-100.
  dZZs = np.abs(np.diff(ZZs, axis=0))

  ll = ZZ.shape[0]-1 

# Add layer for plotting:
  A2d = np.insert(A2d,0,A2d[0,:],axis=0)

# Check is lon is at singular longitude:
# Should happen only when X(1, ...) >0 and after 180E, X<0 (-179,-178,...)
  if XX[0] > XX[-1]:
    dmm = XX-360.
    XX = np.where(XX >= 0.,dmm,XX)

  if f_intract:
    plt.ion()
  else:
    plt.ioff()

  alf=2.
  if not rmin:
    rmin = np.nanpercentile(A2d,alf)
    rmin = round(rmin,2)
    rmax = np.nanpercentile(A2d,100.-alf)
    rmax = round(rmax,2)
    print('rmin={0}, rmax={1}'.format(rmin,rmax))

  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.25, 0.8, 0.7])
  #ax1.plot(XX,Hb)
  im1 = ax1.pcolormesh(XX, ZZs, A2d, shading='flat',
                cmap=cmpr, vmin=rmin, vmax=rmax)
# im1.set_clim(34.5,35.5)
  # Patch bottom:
  verts = [(np.min(XX),-8000),*zip(XX,Hb),(np.max(XX),-8000)]
  poly = Polygon(verts, facecolor='0.1', edgecolor='0.1')
  ax1.add_patch(poly)


  # Plot interfaces
  for kk in range(ll+1):
    dmm = ZZs[kk,:]
    if kk < lrstart:
      continue

    cclr=[0,0,0]
#    ax1.plot(XX,dmm,color=cclr)
    ax1.plot(XX, dmm, color=cclr, linewidth=0.5)

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
        if cc%dcntr == 0 and dzz > 1.5:
          zm = zlr+0.5*dzz
          txt = '{0}'.format(kk)
          ax1.text(XX[iB],zm,txt,fontsize=10)

  #ax1.autoscale(enable=True, axis='both', tight='True')
  ax1.set_xlim([np.min(XX), np.max(XX)])
  ax1.set_ylim([Bmin, 0])

  ax1.set_title(stl)

# Colorbar
#  ax2 = fig1.add_axes([ax1.get_position().x0, ax1.get_position().y0-0.1,
#                     ax1.get_position().width,0.02])
  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')
  

# Plot map:
  import mod_utils as mutil
  LMSK = HH.copy()
  LMSK = np.where(LMSK<0.,0.,1.)
  lcmp = mutil.clrmp_lmask(2, clr_land=[0.5,0.5,0.5])
  ax3 = plt.axes([0.75, 0.02, 0.2, 0.2])
  ax3.pcolormesh(LMSK, shading='flat',\
                  cmap=lcmp, \
                  vmin=0, \
                  vmax=1)
  ax3.contour(HH, [-5000,-4000,-3000,-2000,-1000], \
              colors=[(0.9,0.9,0.9)], linestyles='solid',linewidths=1)
  ax3.plot([i1,i1], [j1,j2], '-', color=[1.,0.2,0])

  dmin  = 300
  dh    = int(max([2*dmin-abs(jx2-jx1),60])*0.5)
  dii   = max([abs(jx2-jx1),dmin])
  ilim1 = max([ix0-dii,0])
  ilim2 = min([ix0+dii,idm])
  jlim1 = max([jx1-dh,0])
  jlim2 = min([jx2+dh,jdm])

  ax3.axis('scaled')
  ax3.set_xlim([ilim1,ilim2])
  ax3.set_ylim([jlim1,jlim2])
# Set y, x axis not visible
#ahd = plt.gca()
  xax = ax3.axes.get_xaxis()
  xax = xax.set_visible(False)
  yax = ax3.axes.get_yaxis()
  yax = yax.set_visible(False)

# Info text:
  if sinfo:
    ax4 = plt.axes([0.1,0.08,0.5,0.15])
    ax4.text(0.0,0.1,sinfo)
    ax4.axis('off') 

  if btx:
    bottom_text(btx,pos=[0.02, 0.03])

  return ZZs, Hb, XX, dZZs

def plot_profile(ax0, S1, zzm, z0, stl='Prof', zlim=-500., xl1=10., xl2=40.):
  """
    Plot T/S profile
  """
#  ax0 = plt.axes([0.1, 0.08, 0.3, 0.35])
  ax0.plot(S1,zzm,'.-',color=(0.,0.4,0.8))
#  ax0.plot(si,z0,'r.')
  ax0.set_xlim(xl1,xl2)
  ax0.set_ylim(zlim,0)
  ax0.plot([xl1, xl2],[z0, z0],'r--')
  ax0.grid(True)
  ax0.set_title(stl)

  return()

