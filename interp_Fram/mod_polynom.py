"""
  Poynomial algorithms and utilities
"""
import os
import numpy as np
from copy import copy
import importlib
import matplotlib.pyplot as plt
import sys

def read_vflux(yr, expt, strnm, pthout, res='008'):
# pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/'
# Saved in matlab code: save_VolFWTflux.m
# /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/anls_Fram2
#
  if res == "008":
    flxfile = pthout+'hycom008_{0:03d}_Fluxes_'.format(expt)+strnm+'_{0}.dat'.format(yr)
  else:
    flxfile = pthout + 'hycom004_{0:03d}_Fluxes_'.format(expt) + strnm + '_{0}.dat'.format(yr)

  print('Reading ',flxfile)

  fid = open(flxfile,'rb')
  fid.seek(0)
  dmm = np.fromfile(fid,dtype='>i',count=1)
  nflds = dmm[0]
  dmm = np.fromfile(fid,dtype='>i',count=1)
  npnts = dmm[0]
  dL = np.fromfile(fid,dtype='>f',count=npnts)
  Vflx = np.fromfile(fid,dtype='>f',count=npnts)
  Tflx = np.fromfile(fid,dtype='>f',count=npnts)
  FWflx = np.fromfile(fid,dtype='>f',count=npnts)
  fid.close()

  Vflx = 1.e-6*Vflx  # Sv
  return dL, Vflx, Tflx, FWflx

def read_coord(pthout,flcoord):
  flpcrd = pthout + flcoord
  print('Reading ',flpcrd)

  fid = open(flpcrd,'rb')
  fid.seek(0)
  dmm   = np.fromfile(fid, dtype='>i', count=1)
  nflds = dmm[0]
  dmm   = np.fromfile(fid, dtype='>i', count=1)
  npnts = dmm[0]
  dL    = np.fromfile(fid, dtype='>f', count=npnts)
  Lonx  = np.fromfile(fid, dtype='>f', count=npnts)
  Latx  = np.fromfile(fid, dtype='>f', count=npnts)
  Indx  = np.fromfile(fid, dtype='>f', count=npnts)
  Jndx  = np.fromfile(fid, dtype='>f', count=npnts)
  fid.close()

  return dL, Lonx, Latx, Indx, Jndx

def read_sect_grid(fsect_grid):
  fid = open(fsect_grid,'rb')
  fid.seek(0)
  dmm = np.fromfile(fid,dtype='>i',count=1)
  nflds = dmm[0]
  dmm = np.fromfile(fid,dtype='>i',count=1)
  npnts = dmm[0]
  Long  = np.fromfile(fid, dtype='>f', count=npnts)
  Latd  = np.fromfile(fid, dtype='>f', count=npnts)
  Hbtm  = np.fromfile(fid, dtype='>f', count=npnts)
  fid.close()

  return Hbtm, Long, Latd

def anlfn(nn):
  """
    Create analytical function for testing 
    best polynomila node searching algorithms
  """
  II   = np.arange(0, nn, dtype=int)
  dx   = 1.e3
  dL   = np.zeros((nn))+dx
  XX   = np.arange(0,nn)*dx
  i1   = int(np.floor(0.3*nn))
  x1   = XX[i1]
  i2   = int(np.floor(0.65*nn))
  x2   = XX[i2]
  V0   = 1.
  sgm  = 10.*dx
  Vflx = V0*np.exp(-(XX-x1)**2/(sgm**2)) - \
         0.55*V0*np.exp(-(XX-x1)**2/(3.*sgm**2))

  Tflx  = Vflx
  FWflx = Vflx

  return dL, Vflx, Tflx, FWflx 

def pcws_lagr1(Xp,Yp,xx):
  """
   Piecewise linear Lagr. Polynomial degree 1
  """
 # Find interval that contains point xx
  dmm = abs(Xp-xx)
  imin = np.argmin(dmm)
  if Xp[imin] < xx :
    i1 = imin
    i2 = i1+1
  elif Xp[imin] == xx:  # node point, to avoid left/right boundary i-1/i+1 problem
    i1 = max([0,imin-1])
    i2 = i1+1
  else:
    i1 = imin-1
    i2 = imin

 # Check dim if i1/i2 >/< max/min size Xp
  if i1 < 0 or i2 > Xp.shape[0]:
    print('ERROR: i1 {0} or i2 {1} is out of bound'.format(i1,i2))
    sys.exit(1)

  Pn = Yp[i1]*(xx-Xp[i2])/(Xp[i1]-Xp[i2]) + \
       Yp[i2]*(xx-Xp[i1])/(Xp[i2]-Xp[i1])

  return Pn

def pcws_lagr2(Xp,Yp,xx):
  """
   Piecewise Lagr. Polynomial degree 2
   needs 3 nodes, near the boundary if < 3 nodes - use linear polynom
  """
# for 1st and last nodes - return the node value:
  if xx == Xp[0]:
    Pn = Yp[0]
    return Pn
  elif xx == Xp[-1]:
    Pn = Yp[-1]
    return Pn

# Define non-overlapping intervals containing 3 nodes:
  nx = Xp.shape[0]
  Is = []
  Xs = []
  for ii in range(0,nx,2):
    if ii == nx-1:   # exact # of segments with 3 nodes
      break

    if ii+2 < nx: 
      Is.append([ii,ii+1,ii+2])
      Xs.append([Xp[ii],Xp[ii+1],Xp[ii+2]]) 
    else:
#  Last segment has less than 3 nodes
      Is.append([ii,ii+1,ii+1])
      Xs.append([Xp[ii],Xp[ii+1],Xp[ii+1]])

  Is  = np.array(Is)
  Xs  = np.array(Xs)  
  nIs = Xs.shape[0]
# Find interval that contains point xx
  for jj in range(nIs):
    if xx > Xs[jj,0] and xx <= Xs[jj,2]:
      iIs = jj
      break
  i1 = Is[iIs,0]
  i2 = Is[iIs,1]
  i3 = Is[iIs,2]
  if not i2 == i3:
    phi1 = ((xx-Xp[i2])*(xx-Xp[i3]))/((Xp[i1]-Xp[i2])*(Xp[i1]-Xp[i3]))
    phi2 = ((xx-Xp[i1])*(xx-Xp[i3]))/((Xp[i2]-Xp[i1])*(Xp[i2]-Xp[i3]))
    phi3 = ((xx-Xp[i1])*(xx-Xp[i2]))/((Xp[i3]-Xp[i1])*(Xp[i3]-Xp[i2]))
    Pn   = Yp[i1]*phi1 + Yp[i2]*phi2 + Yp[i3]*phi3 
  else:
# linear interpolation
    Pn = Yp[i1]*(xx-Xp[i2])/(Xp[i1]-Xp[i2]) + \
         Yp[i2]*(xx-Xp[i1])/(Xp[i2]-Xp[i1])

  return Pn

def pcws_hermite(Xp,Yp,xx):
  """
    Piecewise Hermite polynomial Newton form
    taken as a cubic to satisfy the 4 constraints:
    on the Left: derivative at x(i) and f(x(i)) 
    on the right: derivative at x(i+1) and f(x(i+1))
  """
# for 1st and last nodes - return the node value:
  if xx == Xp[0]:
    Pn = Yp[0]
    return Pn
  elif xx == Xp[-1]:
    Pn = Yp[-1]
    return Pn

# Find interval for xx:
  nx = Xp.shape[0]
  iL = max(np.where(Xp <= xx)[0])
  iLm1 = max([iL-1,0])
  iLp1 = min([iL+1,nx-1])
  iR = iL+1
  iRm1 = max([iR-1,0])
  iRp1 = min([iR+1,nx-1])

# Estimate derivative:
  dfL = (Yp[iLp1]-Yp[iLm1])/(Xp[iLp1]-Xp[iLm1])
  dfR = (Yp[iRp1]-Yp[iRm1])/(Xp[iRp1]-Xp[iRm1])
  
# Basis functions:
  psiL = ((xx-Xp[iR])**2)/((Xp[iL]-Xp[iR])**2)*\
         (1.-2.*(xx-Xp[iL])/(Xp[iL]-Xp[iR]))
  psiR = ((xx-Xp[iL])**2)/((Xp[iR]-Xp[iL])**2)*\
         (1.-2.*(xx-Xp[iR])/(Xp[iR]-Xp[iL]))
  xiL  = ((xx-Xp[iR])**2)/((Xp[iL]-Xp[iR])**2)*(xx-Xp[iL])
  xiR  = ((xx-Xp[iL])**2)/((Xp[iR]-Xp[iL])**2)*(xx-Xp[iR])

# Hermite polynom on the interval [i,i+1]:
  Pn = Yp[iL]*psiL + dfL*xiL + Yp[iR]*psiR + dfR*xiR

  return Pn 


def plot_flx(X,F,Xp,Yp,Xest,Peqd,ERR,n1,n2,\
             stl=' title ',xlbl='East Dist, km', ylbl='Flux, Sv'):
  plt.ion()  # enables interactive mode
  fig = plt.figure(figsize=(8,8))
  fig.clf()
  ax = plt.axes([0.1, 0.6,0.85,0.35])
  ax.plot(X,F,label='Obs')
  ax.plot(Xp,Yp,'r*',label='Interp Points')
  ax.plot(Xest,Peqd,label='Pn(x)')
  ax.legend(bbox_to_anchor=(0.98,-0.1))
  ax.set_xlabel(xlbl)
  ax.set_ylabel(ylbl)
  plt.grid()
  ax.set_title(stl)

 # Plot Error
  xr = np.arange(n1,n2+1)
  ax2 = plt.axes([0.1,0.1,0.8,0.35])
  ax2.plot(xr,ERR,'.-')
  ax2.set_yscale("log")
  ax2.set_xlabel('# of intrp points')
  ax2.set_ylabel('||Err||_{inf}')
  ax2.set_title('Error vs N. of interpolation points')
  plt.grid()

  btx = 'polynom_obs.py'
  bottom_text(btx)


def norm2err(Xp,Yobs,iXp,P1):
  """
    Compute L2 error
    Xp  - Polyn nodes
    iXp - global indices of polynom nodes
    
  """
  Err2 = (Yobs-P1)**2
  nsgm = Xp.shape[0]-1
  ErrSgm = np.zeros((nsgm))
  for iis in range(nsgm):
    iL = iXp[iis]
    iR = iXp[iis+1]-1
    ErrSgm[iis] = np.sum(Err2[iL:iR+1])

# Global error
#  err2nrm = np.sqrt(np.dot((Yex-Plagr1),(Yex-Plagr1)))
  err2nrm = np.sum(Err2)

  return ErrSgm, err2nrm, Err2 


def poly1eqd(XX,YY,N):
  """
    1st degree polynomial with equidistant nodes
    Given N - # of nodes (including 2 end-points)
    return estimates at XX locations
  """
  nP = XX.shape[0]

# Equidistant nodes:
  dx = (nP-1)/(N-1)
  iXp = np.round(np.arange(0,nP-0.1,dx)).astype(int)
  iXp = list(iXp)

  Yp = YY[iXp]
  Xp = XX[iXp]
  Plagr = []  # estimated flux for polynomial points
  for ip in range(nP):
    yip1 = pcws_lagr1(Xp,Yp,XX[ip])
    Plagr.append((yip1))
  Plagr = np.array(Plagr)

  return Plagr, Xp, iXp


def poly1maxerr(XX,YY,N):
  """
    Find best node locations by placing the nodes
    at the largest errors
  """
  nP = XX.shape[0]
# Distributing points to reduce error
# Start with a straight line and 2 end points:
  iXp = [0, nP-1]
  Yp = YY[iXp]
  Xp = XX[iXp]

  Plagr1 = []  # estimated flux for polynomial points
  for ip in range(nP):
    yip1 = pcws_lagr1(Xp,Yp,XX[ip])
    Plagr1.append((yip1))
  Plagr1 = np.array(Plagr1)

# 2-norm
# Error by segments:
# Global error
#  err2nrm = np.sqrt(np.dot((Yex-Plagr1),(Yex-Plagr1)))
  ErrSgm, ErrGlb, ERR2 = norm2err(Xp,YY,iXp,Plagr1)

  f_plt = False
  if f_plt:
    fig1 = plt.figure(1)
    plt.clf()
    plt.plot(XX,YY,color=clr1)
    plt.plot(Xp,Yp,'ro')
    plt.plot(XX,Plagr1,color=clr2) #interpolated best
    plt.grid()
    ctl = '{0} ErrGlb0={1:8.6e}'.format(FlxNm,ErrGlb)
    plt.title(ctl)
#    bottom_text(btx,pos=[0.05,0.02])

  ErrS0 = ErrSgm
  ErrG0 = ErrGlb
  Err20 = ERR2
# Position nodes N-2 at the location of max error
  iXpn = iXp
  for kk in range(2,N):
    jMax = np.argmax(ERR2)
    iXpn.append(jMax)
    iXpn.sort()
    Yp = YY[iXpn]
    Xp = XX[iXpn]
    Plagr1 = []  # estimated flux for equidistant polynomial points
    for ip in range(nP):
      yip1 = pcws_lagr1(Xp,Yp,XX[ip])
      Plagr1.append((yip1))
    Plagr1 = np.array(Plagr1)
    ErrSgm, ErrGlb, ERR2 = norm2err(Xp,YY,iXp,Plagr1)

    if f_plt:
      plt.plot(XX,Plagr1) #interpolated
  
  return Plagr1, Xp, iXp

def poly1minmaxErrGlb(XX,YY,N):
  """
    Find best node iC by minimizing global error
    slide the node from iStart to iEnd
    excluding existing nodes
    Start with all nodes located at the max Error points
    then each point move to the best position 
    that minimizes GlbError
  """
#
  nP = XX.shape[0]
# Distributing points to reduce error
# Start with points at maxError 
  Plagr1,Xp,iXp = poly1maxerr(XX,YY,N)
  ErrSgm, ErrGlb, ERR2 = norm2err(Xp,YY,iXp,Plagr1)
  Yp = YY[iXp]

# Each point (except the end points) 
# move within the segment
# bounded by 2 adjacent points
  for kk in range(1,N-1):
    errG  = []
    Jindx = []
    i1    = iXp[kk-1]
    i2    = iXp[kk+1]
    i0    = iXp[kk]
    for jp in range(i1+1,i2):
      iXpn = iXp.copy()
      iXpn[kk] = jp
      iXpn.sort()
      XpN = XX[iXpn]
      YpN = YY[iXpn]
      P1  = []
      for ip in range(nP):
        yip1 = pcws_lagr1(XpN,YpN,XX[ip])
        P1.append(yip1)
      P1 = np.array(P1)
      ErrSgm, ErrGlb, ERR2 = norm2err(Xp,YY,iXpn,P1)
      errG.append(ErrGlb)
      Jindx.append(jp)

# Pick smallest error:
    errG    = np.array(errG)
    jmin    = np.where(errG == np.min(errG))[0][0]
    iXp[kk] = Jindx[jmin]
    iXp.sort()

# Final selection:
  Yp = YY[iXp]
  Xp = XX[iXp]
  Plagr1 = []  # estimated flux for polynomial points
  for ip in range(nP):
    yip1 = pcws_lagr1(Xp,Yp,XX[ip])
    Plagr1.append((yip1))
  Plagr1 = np.array(Plagr1)
  ErrSgm, ErrGlb, ERR2 = norm2err(Xp,YY,iXp,Plagr1)

  return Plagr1, Xp, iXp


def poly1minErrGlb(XX,YY,N,dlt_eps=0.1, niter0=50, fguess='equidst'):
  """
    Find best node iC by minimizing global error
    iteratively
    Start node distribution as equidistant

    Iteration: 
    for each node, move it between adjacent nodes finding
    the location that minimizes GlbErr

    Repeat until GlbErr does not change
  """
#
# Start with equidistant nodes
  nP = XX.shape[0]
  if fguess == 'equidist':
    Plagr1,XpE,iXp = poly1eqd(XX,YY,N)
  elif fguess == 'maxerr':
    Plagr1,Xp,iXp = poly1maxerr(XX,YY,N)
  else:
    raise Exception('First guess method is {0}, not equidist or maxess'.\
                    format(fguess))
 
  Xp = XX[iXp]
  ErrSgm, ErrGlb, ERR2 = norm2err(Xp,YY,iXp,Plagr1)

# Iteration 
  dErr = 1e20
  ErrG_old = ErrGlb
  niter = 0
  while dErr > dlt_eps and niter < niter0:
# Each point (except the end points) 
# move within the segment
# bounded by 2 adjacent points
# to minimize GlbErr
    niter += 1
    for kk in range(1,N-1):
      errG  = []
      Jindx = []
      i1    = iXp[kk-1]
      i2    = iXp[kk+1]
      i0    = iXp[kk]
      for jp in range(i1+1,i2):
        iXpn = iXp.copy()
        iXpn[kk] = jp
        iXpn.sort()
        XpN = XX[iXpn]
        YpN = YY[iXpn]
        P1  = []
        for ip in range(nP):
          yip1 = pcws_lagr1(XpN,YpN,XX[ip])
          P1.append(yip1)
        P1 = np.array(P1)
        ErrSgm, ErrGlb, ERR2 = norm2err(Xp,YY,iXpn,P1)
        errG.append(ErrGlb)
        Jindx.append(jp)

  # Pick smallest error:
      errG    = np.array(errG)
      jmin    = np.where(errG == np.min(errG))[0][0]
      iXp[kk] = Jindx[jmin]
      iXp.sort()

  # Final selection:
    Yp = YY[iXp]
    Xp = XX[iXp]
    Plagr1 = []  # estimated flux for polynomial points
    for ip in range(nP):
      yip1 = pcws_lagr1(Xp,Yp,XX[ip])
      Plagr1.append((yip1))

    Plagr1 = np.array(Plagr1)
    ErrSgm, ErrGlb, ERR2 = norm2err(Xp,YY,iXp,Plagr1)
    dErr = ErrG_old - ErrGlb
    ErrG_old = ErrGlb
    print(' Min Glb Err, iteration {0}, dErr={1:9.5e}'.format(niter,dErr))

  return Plagr1, Xp, iXp



def poly1minE1E2(XX,Flx,N):
  """
    Find best node iC by minimizing error left segm & error right segm
    slide the node from iL to iR
    not finished - 
  """
#
# Start with equidistant nodes
  nP = XX.shape[0]

# Equidistant nodes:
  dx = (nP-1)/(N-1)
  iXp = np.round(np.arange(0,nP-0.1,dx)).astype(int)
  iXp = list(iXp)

  Yp = YY[iXp]
  Xp = XX[iXp]
  Plagr = []  # estimated flux for polynomial points
  for ip in range(nP):
    yip1 = pcws_lagr1(Xp,Yp,XX[ip])
    Plagr.append((yip1))
  Plagr = np.array(Plagr)
  
  for isgm in range (1,N):
    iL = iXp[isgm-1]
    iR = iXp[isgm+1]
#
# Slide iC from iL to iR and compute error for left and right segments
    errL = []
    errR = []
    errG = []
    YLR  = YY[iL:iR+1]
    XLR  = XX[iL:iR+1]
    for icc in range(iL+1, iR):
# Polynom nodes
      iXpS = [iL,icc,iR]
      XpS = XX[iXpS]
      YpS = YY[iXpS]
      P1 = []
      for ip in range(iL,iR+1):
        yip1 = pcws_lagr1(XpS,YpS,XX[ip])
        P1.append((yip1))
      P1 = np.array(P1)
      ErrSgm, ErrGlb, ERR2 = norm2err(XpS,YLR,iXp,P1)
      eL = ErrSgm[0]/(icc-iL+1)
      eR = ErrSgm[1]/(iR-icc+1)

      errG.append(ErrGlb)
      errL.append((eL))
      errR.append((eR))
 
    errL = np.array(errL)
    errR = np.array(errR) 
    errG = np.array(errG)


  return P1

class poly_stat():
  kind = 'Polynomial best nodes'

  def __init__(self):
    self.year     = []
    self.XX       = []
    self.Nnodes   = []
    self.Xnodes   = []
    self.Inodes   = []
    self.ErrGlb   = []

  def add_data(self, year, Nnodes, Xnodes, Inodes, ErrGlb):
    self.year.append(year)
    self.Nnodes.append(Nnodes)
    self.Xnodes.append(Xnodes)
    self.Inodes.append(Inodes)   
    self.ErrGlb.append(ErrGlb)


def extract_indx_nodes(PNOD,Nnds):
  """
    Extract polynom nodes for saved years
    Nnds = # of polynom. nodes
    PNOD - structured array class poly_stat 
  """
  dmm = np.array(PNOD.Nnodes)
  nd1 = dmm[0]
  nd2 = np.max(dmm)

  if Nnds < nd1 or Nnds > nd2:
    raise Exception(" Requested # nodes {2} is not >= {0} and <= {1}".\
                    format(nd1,nd2,Nnds))

  II = np.where(dmm==Nnds)[0]
  nI = II.shape[0]
  

  IX  = np.zeros((nI,Nnds))-999
  icc = -1
  for k in range(nI):
    i0  = II[k]
    nds = PNOD.Inodes[i0] 
    icc += 1
    IX[icc,:] = nds

  mnIX  = np.mean(IX, axis=0)
  stdIX = np.std(IX, axis=0)

  return mnIX, stdIX





