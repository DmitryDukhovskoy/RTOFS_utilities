"""
  Functions to plot vectors
"""
import numpy as np
import importlib
import sys

sys.path.append('/home/ddmitry/codes/MyPython/hycom_utils')
sys.path.append('/home/ddmitry/codes/MyPython')

def rotate_vector_cart2polar(U, V, xU, yU, lonPol, latPol):
  """
# Project vector from cartesian grid to polar
# given components of the vector U(u,v) are oriented
# along positive Y and positive X of the local grid
# Need to porject this vectors onto new grid
#
# U,V - vector components on cartesian grid
# xU, yU - coordinates of cartesian grid points 
#          on the polar grid
# xPol, yPol - grid points dstiances (coordinates) of the polar grid
# lonPol, latPol - geogr coordinates corresponding to xPol, yPol
# 
  """
# First, locate the North Pole on polar grid
  dm1,dm2 = np.where(latPol == np.max(latPol))
  imx=dm1[0]
  jmx=dm2[0]
# 
# Find 4 points surrounding N. Pole
# For this grid, the N. Pole is in the center of the grid
# if not need to develop search algorithm
# Here, just check that it is in the center of the grid:
  ln1=lonPol[jmx,imx]
  lnm1=lonPol[jmx-1,imx-1]
  lnp1=lonPol[jmx+1,imx+1]

  IJ = [[imx,jmx]]
  if np.abs(ln1-lnm1) > 100.:
#    grx = latPol[jmx+1,imx+1]-latPol[jmx+1,imx]
#    gry = latPol[jmx+1,imx+1]-latPol[jmx,imx+1]
    IJ.append([jmx-1,imx])
    IJ.append([jmx-1,imx-1])
    IJ.append([jmx,imx-1])
  else:
    IJ.append([jmx+1,imx])
    IJ.append([jmx+1,imx+1])
    IJ.append([jmx,imx+1])
#
  """
  Rotate U vectors by finding the local N.P. gradient
  i.e. grad(lat) and finding the angle between this grad
  and N P. direction on the Merc. grid
  """


  ny = U.shape[0]
  nx = U.shape[1]
  Ur = np.zeros((ny,nx))*np.nan
  Vr = np.zeros((ny,nx))*np.nan

  for ii in range(nx):
    for jj in range(ny):
      x0 = xU[jj,ii]
      y0 = yU[jj,ii]
      u0 = U[jj,ii]
      v0 = V[jj,ii]
      if (ii > 0 and jj > 0) and (ii < nx-1 and jj < ny-1):
        dx = (xU[jj,ii+1]-xU[jj,ii-1])*1.e-3
        dphi_dx = (latPol[jj,ii+1]-latPol[jj,ii-1])/dx
        dy = (yU[jj+1,ii]-yU[jj-1,ii])*1.e-3
        dphi_dy = (latPol[jj+1,ii]-latPol[jj-1,ii])/dy
      elif ii == 0 and jj < (ny-1):
        dx = (xU[jj,ii+1]-xU[jj,ii])*1.e-3
        dphi_dx = (latPol[jj,ii+1]-latPol[jj,ii])/dx
        dy = (yU[jj+1,ii]-yU[jj,ii])*1.e-3
        dphi_dy = (latPol[jj+1,ii]-latPol[jj,ii])/dy
      elif ii == (nx-1) and jj < (ny-1):
        dx = (xU[jj,ii]-xU[jj,ii-1])*1.e-3
        dphi_dx = (latPol[jj,ii]-latPol[jj,ii-1])/dx
        dy = (yU[jj+1,ii]-yU[jj,ii])*1.e-3
        dphi_dy = (latPol[jj+1,ii]-latPol[jj,ii])/dy
      elif ii < (nx-1) and jj == 0:
        dx = (xU[jj,ii+1]-xU[jj,ii])*1.e-3
        dphi_dx = (latPol[jj,ii+1]-latPol[jj,ii])/dx
        dy = (yU[jj+1,ii]-yU[jj,ii])*1.e-3
        dphi_dy = (latPol[jj+1,ii]-latPol[jj,ii])/dy
      elif ii < (nx-1) and jj == (ny-1):
        dx = (xU[jj,ii+1]-xU[jj,ii])*1.e-3
        dphi_dx = (latPol[jj,ii+1]-latPol[jj,ii])/dx
        dy = (yU[jj,ii]-yU[jj-1,ii])*1.e-3
        dphi_dy = (latPol[jj,ii]-latPol[jj-1,ii])/dy
      else:
        dx = (xU[jj,ii]-xU[jj,ii-1])*1.e-3
        dphi_dx = (latPol[jj,ii]-latPol[jj,ii-1])/dx
        dy = (yU[jj,ii]-yU[jj-1,ii])*1.e-3
        dphi_dy = (latPol[jj,ii]-latPol[jj-1,ii])/dy

#
# Normalize grad phi - N P direction
      gP = np.sqrt(dphi_dx**2+dphi_dy**2)
      dphi_dx = dphi_dx/gP
      dphi_dy = dphi_dy/gP
#
# Find angle to N. Pole
      alf = np.arctan2(dphi_dy,dphi_dx)
      alf_dgr = alf*180./np.pi
      dalf = alf_dgr-90.

      ur, vr = rotate_vector(u0, v0, dalf)
#     if jj == 0 and ii == 0:
#       print(" jj,ii = 0: ur=", ur," vr=",vr)

      Ur[jj,ii] = ur
      Vr[jj,ii] = vr

  f_chck = False
  if f_chck:
    import mod_draw_vector
#   importlib.reload(mod_draw_vector)
    from mod_draw_vector import compass
    nf = 3
    arrowprops = dict(color='darkorange', linewidth=2)
    ax = compass(dphi_dx, dphi_dy, arrowprops, nf)

  return Ur, Vr


def inpolygon(xq, yq, xv, yv):
  """ 
  Function similar to matlab inpolygon
  from interent stackoverflow
  returns in indicating if the query points specified by xq and yq 
  are inside or on the edge of the polygon area defined by xv and yv.
  """
  from matplotlib import path
  shape = xq.shape
  xq = xq.reshape(-1)
  yq = yq.reshape(-1)
  xv = xv.reshape(-1)
  yv = yv.reshape(-1)
  q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
  p = path.Path([(xv[i], yv[i]) for i in range(xv.shape[0])])

  return p.contains_points(q).reshape(shape)
#
# 
def rotate_vector(uin,vin,thtd):
  """
  Rotate vector U(uin,vin)
  by angle thtd - in degrees
  """
  tht = thtd*np.pi/180.
  R = [np.cos(tht), -np.sin(tht), np.sin(tht), np.cos(tht)] # rotates c/clkws if tht<0
  R = np.array(R).reshape(2,2)
  UV = np.array([uin,vin]).reshape(2,1)
  UVr = R.dot(UV)

  ur = UVr[0].item()
  vr = UVr[1].item()

  """
  nf = 3
  arrowprops = dict(color='darkorange', linewidth=2)
  ax = compass(uin, vin, arrowprops, nf)

  nf = 4
  arrowprops = dict(color='blue', linewidth=2)
  ax = compass(ur, vr, arrowprops, nf)
  """

  return ur, vr


def hycom_dir_gradNorth(lonc,latc,dgr=True):
  """
  Find grad(lat) on HYCOM grid
  return direction of grad(lat) to N.Pole at every grid point
  if dgr=true - in degrees
  """
  nx = lonc.shape[1]
  ny = lonc.shape[0]

  DX, DY = dx_dy_2Dgrid(lonc,latc)

  ALFN = np.zeros((ny,nx))*np.nan
  print("Calculating North direction on HYCOM grid ...")
  for ii in range(nx):
    for jj in range(ny):
        dphi_dx, dphi_dy = dphi(nx,ny,jj,ii,DX,DY,latc,fDX=True)

# Normalize grad phi - N P direction
        gP = np.sqrt(dphi_dx**2+dphi_dy**2)
        dphi_dx = dphi_dx/gP
        dphi_dy = dphi_dy/gP
#
# Find angle to N. Pole on Polar grid:
        alf = np.arctan2(dphi_dy,dphi_dx)
        if dgr: alf_dgr = alf*180./np.pi
        ALFN[jj,ii] = alf_dgr

#       if ii==553 and jj==599:
#         breakpoint()
  return ALFN
#
# 
def dx_dy_2Dgrid(LON,LAT):
  """
  Calculate grid spacing for a 2D grid 
  LON/LAT - grid coordinates
  """
  nn = LON.shape[1]
  mm = LON.shape[0]
  DX = np.zeros((mm,nn))
  DY = np.zeros((mm,nn))

  print('Calculating DX, DY {0}x{1}'.format(mm,nn))

  import mod_misc1
  from mod_misc1 import dist_sphcrd

  for ii in range(nn-1):
    dx = dist_sphcrd(LAT[:,ii],LON[:,ii],LAT[:,ii+1],LON[:,ii+1])
    DX[:,ii] = dx
  DX[:,nn-1] = dx

  for jj in range(mm-1):
    dy = dist_sphcrd(LAT[jj,:],LON[jj,:],LAT[jj+1,:],LON[jj+1,:])
    DY[jj,:] = dy
  DY[mm-1,:] = dy

  return DX, DY


def rotate_vector_hycom2polar(U, V, ALFH, xU, yU, lonPol, latPol):
  """
# Project vector from HYCOM/CICE grid to polar
# given components of the vector U(u,v) are oriented
# along positive Y and positive X of the local grid
# Need to porject this vectors onto new grid
#
# U,V - vector components on cartesian grid
# ALFH - direction to the N.Pole on HYCOM grid, degrees
#        0 dgr - along X axis 
#
# xU, yU - coordinates of cartesian grid points 
#          on the polar grid
# xPol, yPol - grid points dstiances (coordinates) of the polar grid
# lonPol, latPol - geogr coordinates corresponding to xPol, yPol
# 
  """
# First, locate the North Pole on polar grid
  dm1,dm2 = np.where(latPol == np.max(latPol))
  imx=dm1[0]
  jmx=dm2[0]
# 
# Find 4 points surrounding N. Pole
# For this grid, the N. Pole is in the center of the grid
# if not need to develop search algorithm
# Here, just check that it is in the center of the grid:
  ln1=lonPol[jmx,imx]
  lnm1=lonPol[jmx-1,imx-1]
  lnp1=lonPol[jmx+1,imx+1]

  IJ = [[imx,jmx]]
  if np.abs(ln1-lnm1) > 100.:
#    grx = latPol[jmx+1,imx+1]-latPol[jmx+1,imx]
#    gry = latPol[jmx+1,imx+1]-latPol[jmx,imx+1]
    IJ.append([jmx-1,imx])
    IJ.append([jmx-1,imx-1])
    IJ.append([jmx,imx-1])
  else:
    IJ.append([jmx+1,imx])
    IJ.append([jmx+1,imx+1])
    IJ.append([jmx,imx+1])
#
  """
  Rotate U vectors by finding the local N.P. gradient
  i.e. grad(lat) and finding the angle between this grad
  and N P. direction on the Merc. grid
  """


  ny = U.shape[0]
  nx = U.shape[1]
  Ur = np.zeros((ny,nx))*np.nan
  Vr = np.zeros((ny,nx))*np.nan
  ALFP = np.zeros((ny,nx))*np.nan

  for ii in range(nx):
    for jj in range(ny):
      x0 = xU[jj,ii]
      y0 = yU[jj,ii]
      u0 = U[jj,ii]
      v0 = V[jj,ii]
      dphi_dx, dphi_dy = dphi(nx,ny,jj,ii,xU,yU,latPol)
#
# Normalize grad phi - N P direction
      gP = np.sqrt(dphi_dx**2+dphi_dy**2)
      dphi_dx = dphi_dx/gP
      dphi_dy = dphi_dy/gP

# Find angle to N. Pole on Polar grid:
      alf = np.arctan2(dphi_dy,dphi_dx)
      alf_dgr = alf*180./np.pi
      ALFP[jj,ii] = alf_dgr
      alfH_dgr = ALFH[jj,ii]
      dalf = alf_dgr-alfH_dgr

      ur, vr = rotate_vector(u0, v0, dalf)
#     if jj == 0 and ii == 0:
#       print(" jj,ii = 0: ur=", ur," vr=",vr)

      Ur[jj,ii] = ur
      Vr[jj,ii] = vr

# breakpoint()

  f_chck = False
  if f_chck:
    import mod_draw_vector
#   importlib.reload(mod_draw_vector)
    from mod_draw_vector import compass
    nf = 3
    arrowprops = dict(color='darkorange', linewidth=2)
    ax = compass(dphi_dx, dphi_dy, arrowprops, nf)

  return Ur, Vr, ALFP


def dphi(nx,ny,jj,ii,xU,yU,latPol, xscl=1.e-3, fDX=False):
  """
  If fDX true - xU and yU are dx and dy 
  """

  if (ii > 0 and jj > 0) and (ii < nx-1 and jj < ny-1):
    if fDX:
      dx = xU[jj,ii]*xscl
    else:
      dx = (xU[jj,ii+1]-xU[jj,ii-1])*xscl
    if fDX:
      dy = yU[jj,ii]*xscl
    else:
      dy = (yU[jj+1,ii]-yU[jj-1,ii])*xscl
    dphi_dx = (latPol[jj,ii+1]-latPol[jj,ii-1])/dx
    dphi_dy = (latPol[jj+1,ii]-latPol[jj-1,ii])/dy
  elif ii == 0 and jj < (ny-1):
    if fDX:
      dx = xU[jj,ii]*xscl
    else:
      dx = (xU[jj,ii+1]-xU[jj,ii])*xscl
    if fDX:
      dy = yU[jj,ii]*xscl
    else:
      dy = (yU[jj+1,ii]-yU[jj,ii])*xscl
    dphi_dx = (latPol[jj,ii+1]-latPol[jj,ii])/dx
    dphi_dy = (latPol[jj+1,ii]-latPol[jj,ii])/dy
  elif ii == (nx-1) and jj < (ny-1):
    if fDX:
      dx = xU[jj,ii]*xscl
    else:
      dx = (xU[jj,ii]-xU[jj,ii-1])*xscl
    if fDX:
      dy = yU[jj,ii]*xscl
    else:
      dy = (yU[jj+1,ii]-yU[jj,ii])*xscl
    dphi_dx = (latPol[jj,ii]-latPol[jj,ii-1])/dx
    dphi_dy = (latPol[jj+1,ii]-latPol[jj,ii])/dy
  elif ii < (nx-1) and jj == 0:
    if fDX:
      dx = xU[jj,ii]*xscl
    else:
      dx = (xU[jj,ii+1]-xU[jj,ii])*xscl
    if fDX:
      dy = yU[jj,ii]*xscl
    else:
      dy = (yU[jj+1,ii]-yU[jj,ii])*xscl
    dphi_dx = (latPol[jj,ii+1]-latPol[jj,ii])/dx
    dphi_dy = (latPol[jj+1,ii]-latPol[jj,ii])/dy
  elif ii < (nx-1) and jj == (ny-1):
    if fDX:
      dx = xU[jj,ii]*xscl
    else:
      dx = (xU[jj,ii+1]-xU[jj,ii])*xscl
    if fDX:
      dy = yU[jj,ii]*xscl
    else:
      dy = (yU[jj,ii]-yU[jj-1,ii])*xscl
    dphi_dx = (latPol[jj,ii+1]-latPol[jj,ii])/dx
    dphi_dy = (latPol[jj,ii]-latPol[jj-1,ii])/dy
  else:
    if fDX:
      dx = xU[jj,ii]*xscl
    else:
      dx = (xU[jj,ii]-xU[jj,ii-1])*xscl
    if fDX:
      dy = yU[jj,ii]*xscl
    else:
      dy = (yU[jj,ii]-yU[jj-1,ii])*xscl
    dphi_dx = (latPol[jj,ii]-latPol[jj,ii-1])/dx
    dphi_dy = (latPol[jj,ii]-latPol[jj-1,ii])/dy

  return dphi_dx, dphi_dy


