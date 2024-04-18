"""
  Functions/subroutines for bilinear interpolation
  Dmitry Dukhovskoy NOAA NESDIS NCEI 
  June 2022
"""
import numpy as np

def basisFn_RectRef():
  """
    Define basis functions on a reference rectangle (square)
    (-1,-1), (1,-1) (1,1) (-1,1)

    For each basis fn phi1, ..., phi4 solve a system of linear eqns
    such that phi1 = 1 at (x1,y1) and 0 at the other vertices, etc.
    a + b*x1 + c*y1 + d*x1*y1 = 1
    a + b*x2 + c*y2 + d*x2*y2 = 0
    etc

    M*A = X, M is [1 -1 -1 1; 
                   1 1 -1 -1; etc, 
    A is (a; b; c; d), X=[1, 0, 0, 0] for phi1 
    See p. 83-90, Gockenbach, "Understanding Finite-Element Methods" textbook
    inverse of M is known
  """
  Minv = 1./4.*np.array([[1, 1, 1, 1],[-1,1,1,-1],[-1,-1,1,1],[1,-1,1,-1]])
  X1 = np.transpose(np.array([[1,0,0,0]]))
  X2 = np.transpose(np.array([[0,1,0,0]]))
  X3 = np.transpose(np.array([[0,0,1,0]]))
  X4 = np.transpose(np.array([[0,0,0,1]]))

  phi1 = np.dot(Minv,X1).squeeze()
  phi2 = np.dot(Minv,X2).squeeze()
  phi3 = np.dot(Minv,X3).squeeze()
  phi4 = np.dot(Minv,X4).squeeze()

  return phi1, phi2, phi3, phi4

def sort_gridcell_indx(II, JJ, fsens='positive'):
  """
    In case when rectangular indices are in random order
    need to put them in c/clckwize (positive) or clockwize (negative)
    order

    assumed: indices represent vertices of a grid box:
    i = 0,1,2,3

     II[i]JJ[i]     II[..]JJ[..]
     *----------------*
     |                |
     |                |
     |                |
     *----------------*
     II[..]JJ[..]     II[..]JJ[..]

   e.g.: II = [10,11,10,11] J=[23, 22, 22, 23]
     want:  II=[10,11,11,10] J=[22,22,23,23]
  """
  if isinstance(II, list):
    II = np.array(II)
  if isinstance(JJ, list):
    JJ = np.array(JJ)
    
# Start from smallest i,j:
  D    = np.sqrt(II**2 + JJ**2)
  isrt = np.argsort(D)
  I0   = II.copy()
  J0   = JJ.copy()  

  ix1         = I0[isrt[0]]
  jx1         = J0[isrt[0]]
  I0[isrt[0]] = -999
  J0[isrt[0]] = -999

  if fsens == 'positive':
# c/clockwise sense:
    imm     = np.argwhere(J0 == jx1)[0][0]
    isrt[1] = imm
    ix2     = I0[imm]
    jx2     = J0[imm]
    I0[imm] = -999
    J0[imm] = -999

    imm     = np.argwhere(I0 == ix2)[0][0]
    isrt[2] = imm
    ix3     = I0[imm]
    jx3     = J0[imm]
    I0[imm] = -999
    J0[imm] = -999
 
    imm     = np.argwhere(J0 == jx3)[0][0]
    isrt[3] = imm

  else:
# clockwise sense:
    imm     = np.argwhere(I0 == ix1)[0][0]
    isrt[1] = imm
    ix2     = I0[imm]
    jx2     = J0[imm]
    I0[imm] = -999
    J0[imm] = -999

    imm     = np.argwhere(J0 == jx2)[0][0]
    isrt[2] = imm
    ix3     = I0[imm]
    jx3     = J0[imm]
    I0[imm] = -999
    J0[imm] = -999

    imm     = np.argwhere(I0 == ix3)[0][0]
    isrt[3] = imm

  IIsrt = II[isrt]
  JJsrt = JJ[isrt]

  return IIsrt, JJsrt

def phi_x0y0(phi,x1,y1):
  """
    Calculate basis fn phi at pnt(x1,y1)
  """
  bb1 = phi[0] + phi[1]*x1 + phi[2]*y1 + phi[3]*x1*y1

  return bb1


def map_x2xhat(XX,YY,x0,y0):
  """
    Map a point from a quadrilateral (x1,y1), ..., (x4,y4) 
    to a reference square
    (-1,-1), ..., (-1,1)
    The mapping is of the form:
    xhat = a1+a2*x+a3*y+a4*x*y = Alf*XY0
    yhat = b1+b2*x+b3*y+b4*x*y = Bet*XY0
    Solve 4 equations to find coefficients:
    a1+a2*x1+a3*y1+a4*x1*y1 = -1 
    a1+a2*x2+a3*y2+a4*x2*y2 = 1 
    ...

    AA*Alf = Cx

   same for b coefficients

   Input: 
     XX, YY - X, Y coordinates of 4 nodes, as 1D vectors of a grid cell
              encompassing the interpolation point
     x0,y0 - coordinates of the interpolation point (same unites as X,Y)
  """
  AA=np.array([[1, XX[0], YY[0], XX[0]*YY[0]],\
               [1, XX[1], YY[1], XX[1]*YY[1]],\
               [1, XX[2], YY[2], XX[2]*YY[2]],\
               [1, XX[3], YY[3], XX[3]*YY[3]]])

#  AA = np.float64(AA)
  AAinv = np.linalg.inv(AA)
  Cx = np.transpose(np.array([[-1.,1.,1.,-1.]]))
  Alf = np.dot(AAinv,Cx)

  Cy = np.transpose(np.array([[-1.,-1.,1.,1.,]]))
  Bet = np.dot(AAinv,Cy)


  XY0 = np.transpose(np.array([[1,x0,y0,x0*y0]]))
  xhat = np.dot(Alf.transpose(),XY0)[0][0]
  yhat = np.dot(np.transpose(Bet),XY0)[0][0]

  return xhat, yhat

def map_xhat2x(XX,YY,xht,yht):
  """
    Inversing mapping from a reference quadrilateral ---> original
    quadrilateral with vertices XX, YY
    The mapping is:
    a1+a2*xhat+a3*yhat+a4*xhat*yhat = x
    b1+b2*xhat+b3*yhat+b4*xhat*yhat = y

    Solve a system of lin eq:
    a1 + a2*s1 + a3*t1 + a4*s1*t1 = x1
    a1 + a2*s2 + a3*t2 + a4*s2*t2 = x2
    ...

    b1 + b2*s1 + b3*t1 + b4*s1*t1 = y1
    b1 + b2*s2 + b3*t2 + b4*s2*t2 = y2
    ...

  """
  XX    = np.expand_dims(XX, axis=0)
  YY    = np.expand_dims(YY, axis=0)

  AA = np.array([[1, -1, -1,  1],\
                 [1,  1, -1, -1],\
                 [1,  1,  1,  1],\
                 [1, -1,  1, -1]])

  AAinv = np.linalg.inv(AA)
  Alf   = np.matmul(AAinv,XX.transpose())
  Bet   = np.matmul(AAinv,YY.transpose())
  ST0   = np.transpose(np.array([[1, xht, yht, xht*yht]]))
  xx0   = np.dot(Alf.transpose(), ST0)[0][0]
  yy0   = np.dot(Bet.transpose(), ST0)[0][0]

  return xx0, yy0

def bilin_interp(phi1,phi2,phi3,phi4,xht,yht,HT):
  """
  Bilinear interpolation H(x,y) = sum(H(x1,y1)*psi_i(x,y)), i=1,...,4
  Input: 
    phi1,phi2,phi3,phi4 - linear basis functions
    xht, yht - interpolated point mapped onto a reference square
               use map_x2xhat   
  """
  p1x = phi_x0y0(phi1,xht,yht) # Phi1 at pnt xht,yht
  p2x = phi_x0y0(phi2,xht,yht)
  p3x = phi_x0y0(phi3,xht,yht)
  p4x = phi_x0y0(phi4,xht,yht)

  Hi = HT[0]*p1x + HT[1]*p2x + HT[2]*p3x + HT[3]*p4x

  return Hi

def bilin_interp1D(phi1,phi2,phi3,phi4,xht,yht,HT):
  """
  Bilinear interpolation H(x,y) = sum(H(x1,y1)*psi_i(x,y)), i=1,...,4
  For 1 D inpyt arrays combined in HT = [[1:n],[1:n],[1:n],[1:n]]
  where n - number of observations at each box vertex, e.g.
  n = number of vertical layers in a model   

  Input: 
    phi1,phi2,phi3,phi4 - linear basis functions
    xht, yht - interpolated point mapped onto a reference square
               use map_x2xhat   
  """
  p1x = phi_x0y0(phi1,xht,yht) # Phi1 at pnt xht,yht
  p2x = phi_x0y0(phi2,xht,yht)
  p3x = phi_x0y0(phi3,xht,yht)
  p4x = phi_x0y0(phi4,xht,yht)

  if not HT.ndim == 2:
    raise Exception('HT should be N x 4 array')

  npnt = HT.shape[1]
  nlv  = HT.shape[0]
  PX  = np.array([[p1x, p2x, p3x, p4x]]).transpose()
  Hi = np.matmul(HT,PX)
#  Hi = (HT[:,0]*p1x + HT[1]*p2x + HT[2]*p3x + HT[3]*p4x)[0]
  Hi = np.squeeze(Hi)

  return Hi

def blinrint(xht,yht,HT):
  """
  Bilinear interpolation H(x,y) = sum(H(x1,y1)*psi_i(x,y)), i=1,...,4
  Input: 
    xht, yht - interpolated point mapped onto a reference square
               use map_x2xhat   
    HT       - 1D array of node values
  """
# Find basis functions for a reference square:
  phi1,phi2,phi3,phi4 = basisFn_RectRef()

  p1x = phi_x0y0(phi1,xht,yht) # Phi1 at pnt xht,yht
  p2x = phi_x0y0(phi2,xht,yht)
  p3x = phi_x0y0(phi3,xht,yht)
  p4x = phi_x0y0(phi4,xht,yht)

  Hi = (HT[0]*p1x + HT[1]*p2x + HT[2]*p3x + HT[3]*p4x)[0]

  return Hi

def lonlat2xy_wrtX0(XX, YY, x0, y0, dltx=1.):
  """
    Convert lon/lat coordinates of 4 quadrilateral vertices XX, YY
    into cartisian coordinates (meters) wrt to a point (x0,y0)

    so that x0,y0 --> 0,0
    and XX, YY ---> +/- distance from 0,0
  """
  import mod_misc1 as mmisc1

  xrf = x0-1.
  if xrf < -180.:
    xrf = xrf+360.
  yrf = y0-1.
  xd0 = mmisc1.dist_sphcrd(y0,xrf,y0,x0)
  yd0 = mmisc1.dist_sphcrd(yrf,x0,y0,x0)

  XV  = np.zeros((4))
  YV  = np.zeros((4))
  for ipp in range(0,4):
    dlx     = mmisc1.dist_sphcrd(YY[ipp],xrf,YY[ipp],XX[ipp]) - xd0
    dly     = mmisc1.dist_sphcrd(yrf,XX[ipp],YY[ipp],XX[ipp]) - yd0
    XV[ipp] = dlx
    YV[ipp] = dly

  return XV, YV, 0., 0. 







