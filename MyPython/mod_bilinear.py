"""
  Functions/subroutines for bilinear interpolation
  Dmitry Dukhovskoy NOAA NESDIS NCEI 
  June 2022
"""
import numpy as np

def basisFn_RectRef():
  """
    Define basis functions on a reference rectangle
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

  phi1 = np.dot(Minv,X1)
  phi2 = np.dot(Minv,X2)
  phi3 = np.dot(Minv,X3)
  phi4 = np.dot(Minv,X4)

  return phi1, phi2, phi3, phi4


def phi_x0y0(phi,x1,y1):
  """
    Calculate basis fn phi at pnt(x1,y1)
  """
  bb1 = phi[0] + phi[1]*x1 + phi[2]*y1 + phi[3]*x1*y1

  return bb1


def map_x2xhat(XX,YY,x0,y0):
  """
    Map a point from a rectangle (x1,y1), ..., (x4,y4) to a reference recatngle
    (-1,-1), ..., (-1,1)
    The mapping is of the form:
    xhat = a1+a2*x+a3*y+a4*x*y
    yhat = b1+b2*x+b3*y+b4*x*y
    Solve 4 equations to find coefficients:
    a1+a2*x1+a3*y1+a4*x1*y1 = -1 
    a1+a2*x2+a3*y2+a4*x2*y2 = 1 
    ...

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

  AAinv = np.linalg.inv(AA)
  Cx = np.transpose(np.array([[-1.,1.,1.,-1.]]))
  Alf = np.dot(AAinv,Cx)

  Cy = np.transpose(np.array([[-1.,-1.,1.,1.,]]))
  Bet = np.dot(AAinv,Cy)


  XY0 = np.transpose(np.array([[1,x0,y0,x0*y0]]))
  xhat = np.dot(np.transpose(Alf),XY0)[0][0]
  yhat = np.dot(np.transpose(Bet),XY0)[0][0]

  return xhat, yhat

def bilin_interp(phi1,phi2,phi3,phi4,xht,yht,HT):
  """
  Bilinear interpolation H(x,y) = sum(H(x1,y1)*psi_i(x,y)), i=1,...,4
  Input: 
    phi1,phi2,phi3,phi4 - linear basis functions
    xht, yht - interpolated point mapped onto a reference rectangle 
               use map_x2xhat   
  """
  p1x = phi_x0y0(phi1,xht,yht) # Phi1 at pnt xht,yht
  p2x = phi_x0y0(phi2,xht,yht)
  p3x = phi_x0y0(phi3,xht,yht)
  p4x = phi_x0y0(phi4,xht,yht)

  Hi = (HT[0]*p1x + HT[1]*p2x + HT[2]*p3x + HT[3]*p4x)[0]

  return Hi

def blinrint(xht,yht,HT):
  """
  Bilinear interpolation H(x,y) = sum(H(x1,y1)*psi_i(x,y)), i=1,...,4
  Input: 
    xht, yht - interpolated point mapped onto a reference rectangle 
               use map_x2xhat   
    HT       - 1D array of node values
  """
# Find basis functions for a reference rectangle:
  phi1,phi2,phi3,phi4 = basisFn_RectRef()

  p1x = phi_x0y0(phi1,xht,yht) # Phi1 at pnt xht,yht
  p2x = phi_x0y0(phi2,xht,yht)
  p3x = phi_x0y0(phi3,xht,yht)
  p4x = phi_x0y0(phi4,xht,yht)

  Hi = (HT[0]*p1x + HT[1]*p2x + HT[2]*p3x + HT[3]*p4x)[0]

  return Hi



