"""
  1d interpolation algorithms global and picewise 

  Dmitry Dukhovskoy, NOAA NCEP EMC
"""
def lagr_polynom(Xp,Yp,xx):
  """
  Lagrange polynomial
  estimate function at xx

  Pn(x) = sum(i to n): y(i)*Li(x)
  """
  Np = Xp.shape[0]
  Pn = 0.
  for ii in range(Np):
    xi = Xp[ii]
    mi_x = 1.
    mi_xi = 1.
    for jj in range(Np):
      if jj == ii:
        continue
      mi_xi = mi_xi*(xi-Xp[jj])  # mi(xi)
#     print('jj={0}, xi={1}, Xp={2}, mi_xi={3}'.format(jj,xi,Xp[jj],mi_xi))
      mi_x = mi_x*(xx-Xp[jj])  # mi(x)
#     print('xx={0}, mi_x={1}'.format(xx,mi_x))

    Li_x = mi_x/mi_xi
    Pn = Pn+Yp[ii]*Li_x

  return Pn


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











