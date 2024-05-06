"""
  1d interpolation algorithms global and picewise 

  Dmitry Dukhovskoy, NOAA NCEP EMC
"""
import numpy as np
import sys

def lagr_polynom(Xp,Yp,xx):
  """
  Lagrange polynomial
  estimate function at xx
  Given values Yp at interpolation nodes Xp
  find value at xx

  The number of interpolation nodes Xp (n)
  defines the polynomial (degree = n-1)

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
#    print(f"ii={ii} mi_x={mi_x} mi_xi={mi_xi} yi*phi={Yp[ii]*Li_x}")
#    print(f"ii={ii} Y={Yp[ii]} phi={Li_x}")

  return Pn

def lagr_polynom2D(Xp, Yp, xx, axis_interp=0):
  """
  Lagrange polynomial
  estimate function at point xx
  Given values Yp(N,M) at interpolation nodes Xp(N,M)
  find value at xx
  xx is either 1 value or a vector
  len(xx) = N or M depending on axis_interp

  axis_interp = 0: interpolation along Y-axis (N nodes,
  M time series)
  then len(xx) = M

  axis_interp = 1: interpolation along X-axis (M nodes,
  N time series)
  len(xx) = N 

  The number of interpolation nodes Xp (n)
  defines the polynomial (degree = n-1)

  Pn(x) = sum(i to n): y(i)*Li(x)
  """
  try:
    lx = len(xx)
  except:
    lx = 0

  f1d = False
  if lx > 0: f1d = True

  if axis_interp > 0: 
    Xp = Xp.transpose()
    Yp = Yp.transpose()

  Np = Xp.shape[0]
  M  = Xp.shape[1]
  if f1d and not lx == M:
    raise Exception(f"len of xx array {lx} should be {M}")

  
  Pn = np.zeros(M)
  for ii in range(Np):
    xi = Xp[ii,:]
    mi_x = np.ones((M))
    mi_xi = np.ones((M))
    for jj in range(Np):
      if jj == ii:
        continue
      mi_xi = mi_xi*(xi-Xp[jj,:])  # mi(xi)
#     print('jj={0}, xi={1}, Xp={2}, mi_xi={3}'.format(jj,xi,Xp[jj],mi_xi))
      mi_x = mi_x*(xx-Xp[jj,:])  # mi(x)
#     print('xx={0}, mi_x={1}'.format(xx,mi_x))

    Li_x = mi_x/mi_xi
    Pn = Pn+Yp[ii,:]*Li_x

  return Pn

def lagr_polynom1D(Xp,Yp,xx):
  """
  Lagrange polynomial
  estimate function at xx
  Given values Yp[N,M] at interpolation nodes Xp[M]
  find value at xx

  Yp is a 2D array with data at each node in columns, rows - locations
  e.g., Yp = values at depths along a section with nodes at Xp

  Pn(x) = sum(i to n): y(i)*Li(x)
  """
  Np = Xp.shape[0]
  Nz = Yp.shape[0]
  Pn = np.zeros((Nz))
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
    Pn = Pn+Yp[:,ii]*Li_x

  return Pn


def pcws_lagr1(Xp,Yp,xx):
  """
   Piecewise linear Lagr. Polynomial degree 1
   Xp is assumed in Ascending Order
  """
  if np.any(np.isnan(Xp)) or np.any(np.isnan(Yp)):
    raise Exception ('No nans are allowed in the input Xp or Yp')

  if not Xp.shape[0] == Yp.shape[0]:
    raise Exception ('Both Xp and Yp arrays should be same length')

  if not np.all(np.diff(Xp)):
    raise Exception ('Xp should be monotonically increasing, no repeating values')

# Check if Xp is increasing:
  if Xp[0] > Xp[1]:
    Xp = -Xp
    xx = -xx

#  if Xp[1] < Xp[0]:
#    raise Exception('Xp should be increasing')
#    print('Xp is not in ascending order, flipping ...')
#    Xp = np.flip(Xp)
#    Yp = np.flip(Yp)

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

def pcws_lagr2(Xp, Yp, xx):
  """
   Piecewise Lagr. Polynomial degree 2
   needs 3 nodes, near the boundary if < 3 nodes - use linear polynom
   Xp -  nodes
   Yp - values at Xp
   xx - interpolation point
   Pn - interpolated value at point xx

   Assumed: X axis is increasing, if X is revearsed (decreasing order, 
            like depths) X is changed sign to make it increasing
   X should not have repeating values (i.e. monotonically increasing)
   No nan values allowed 
   Xp and Yp should be same length

  """
# for 1st and last nodes - return the node value:
  if xx == Xp[0]:
    Pn = Yp[0]
    return Pn
  elif xx == Xp[-1]:
    Pn = Yp[-1]
    return Pn

  if np.any(np.isnan(Xp)) or np.any(np.isnan(Yp)):
    raise Exception ('No nans are allowed in the input Xp or Yp')

  if not Xp.shape[0] == Yp.shape[0]:
    raise Exception ('Both Xp and Yp arrays should be same length')

  if not np.all(np.diff(Xp)):
    raise Exception ('Xp should be monotonically increasing, no repeating values')

# Check if Xp is increasing:
  if Xp[0] > Xp[1]:
    Xp = -Xp
    xx = -xx
#    Xp = np.flip(Xp)
#    Yp = np.flip(Yp)

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

def pcws_lagr3(Xp, Yp, xx):
  """
   Piecewise Lagr. Polynomial degree 3
   needs 4 nodes, near the boundary if < 4 nodes - use lower-degree polynom
   Xp -  nodes
   Yp - values at Xp
   xx - interpolation point
   Pn - interpolated value at point xx

   Assumed: X axis is increasing, if X is revearsed (decreasing order, 
            like depths) X is changed sign to make it increasing
   X should not have repeating values (i.e. monotonically increasing)
   No nan values allowed 
   Xp and Yp should be same length

  """
  nnd = 4
# for 1st and last nodes - return the node value:
  if xx == Xp[0]:
    Pn = Yp[0]
    return Pn
  elif xx == Xp[-1]:
    Pn = Yp[-1]
    return Pn

 # Check if xx is inside the global interval:
  if abs(xx) > np.max(abs(Xp)):
    raise Exception(f'xx {xx} outisde range: {np.min(Xp)}/{np.max(Xp)}')

  if np.any(np.isnan(Xp)) or np.any(np.isnan(Yp)):
    raise Exception ('No nans are allowed in the input Xp or Yp')

  if not Xp.shape[0] == Yp.shape[0]:
    raise Exception ('Both Xp and Yp arrays should be same length')

  if not np.all(np.diff(Xp)):
    raise Exception ('Xp should be monotonically increasing, no repeating values')

# Check if Xp is increasing:
  if Xp[0] > Xp[1]:
    Xp = -Xp
    xx = -xx

# Define non-overlapping intervals containing 4 nodes:
  nx = Xp.shape[0]
  Is = []
  Xs = []
#  for ii in range(0,nx,nnd-1):
  for ii in range(0,nx):
    if ii == nx-1:   # exact # of segments with 4 nodes
      break
    if ii+nnd-1 > nx:  # end 
      break

    tmp = list([k for k in range(ii,ii+nnd)])
    if ii+nnd-1 < nx:
      Is.append(tmp)
      Xs.append(list(Xp[tmp]))
    else:
#  When x is near the end - use different cardinal basis
# constructed on nodes backward, i.e.
# for y(x0): (x-x(i-1)*(x-x(i-2))*(x-x(i-3))/
# ((x(i)-x(i-1)*(x(i)-x(i-2))*(x(i)-x(i-3)), etc
      tmp = np.arange(nx-nnd,nx)
      tmp = np.flip(tmp)
      tmp = list(tmp)
      Is.append(tmp)
      Xs.append(list(Xp[tmp]))

  Is  = np.array(Is)
  Xs  = np.array(Xs)
  nIs = Xs.shape[0]
# Find interval that contains point xx
# Start with middle interval - avoid end intervals due to large errors
# if xx is in the last interval - check next group, it should be in the middle
# if xx is in the 1st interval - accept (this should be xx close to the 1st pnt)
  for jj in range(nIs):
#print(f'jj={jj} nIs={nIs} xx={xx} min/max = {np.min(Xs[jj,:])}/{np.max(Xs[jj,:])}')
    if xx >= Xs[jj,1] and xx <= Xs[jj,2]:
      iIs = jj
      break
    elif xx > Xs[jj,0] and xx <= Xs[jj,1]:
      iIs = jj
      break
    elif jj == nIs-1:
      if xx > np.min(Xs[jj,:]) and xx <= np.max(Xs[jj,:]):
        iIs = jj
        break

#  print(f"{Xs[jj,1]} > {xx} < {Xs[jj,2]}")
#  print(f"iIs={iIs} Xs={Xs[jj,:]}") 
 
  i1 = Is[iIs,0]
  i2 = Is[iIs,1]
  i3 = Is[iIs,2]
  i4 = Is[iIs,3]
# 3rd degree polynom:
  x0 = Xp[i1]
  x1 = Xp[i2]
  x2 = Xp[i3]
  x3 = Xp[i4]
# Cardinal basis functions:
  phi0 = ((xx-x1)*(xx-x2)*(xx-x3))/((x0-x1)*(x0-x2)*(x0-x3))
  phi1 = ((xx-x0)*(xx-x2)*(xx-x3))/((x1-x0)*(x1-x2)*(x1-x3))
  phi2 = ((xx-x0)*(xx-x1)*(xx-x3))/((x2-x0)*(x2-x1)*(x2-x3))
  phi3 = ((xx-x0)*(xx-x1)*(xx-x2))/((x3-x0)*(x3-x1)*(x3-x2))
  Pn   = Yp[i1]*phi0 + Yp[i2]*phi1 + Yp[i3]*phi2 + Yp[i4]*phi3

#  print(f"x0={x0} x1={x1} x2={x2} x3={x3} xx={xx}")
#  print(f"y0={Yp[i1]} y1={Yp[i2]} y2={Yp[i3]} y3={Yp[i4]}")
#  print(f"phi0={phi0} phi1={phi1} phi2={phi2} phi3={phi3}")
#  print(f"y*p0={Yp[i1]*phi0} y*p1={Yp[i2]*phi1} " +\
#        f"y*p2={Yp[i3]*phi2} y*p3={Yp[i4]*phi3} ")

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


def make_monotonic(Xp, eps0=0.001):
  """
    For interpolation, polynomial nodes have to be
    monotonically increasing / decreasing
    In some cases (below botom), Xp can be not monotonic function
    Change to monotonic

    Xp = 1D array
  """
  if Xp[-1] > Xp[0]:
    MIncr = True
    eps0 = abs(eps0)
  else:
    MIncr = False
    eps0 = -abs(eps0)

  zrr = 1.e-23
  dX  = abs(np.diff(Xp))
  if min(abs(dX)) >= zrr:
    return Xp

  ix0 = min(np.where(dX <= zrr)[0]) 
  nxx = len(Xp)
  for ixx in range(ix0, nxx-1):
    Xp[ixx+1] = Xp[ixx] + eps0

  return Xp


def interp_bins(H1d, zzH, zmH, zzM, zmM, eps_dZ=1.e-3, \
                fill_btm=False, fill_val=0.):
  """
    For interpolation / remapping values from one set of vertical layers
    into another (e.g., HYCOM --> MOM6)

    zzH - input interface depths, m <0
    zzM - interface depths (m, <0) where the values are being interpolated to

    eps_dZ - min layer thickness, if < eps_dZ - bottom layer
    For layers below bottom depth, default - copy last value
    if fill_btm - then bottom layers filled with fill_value 


    Remapping preserves total quantity integrated over all layers
    Remapping done by binning
  

     V. grid 1                    V. grid 2

    k=5 ------  z=1      k=8 --------  z=1
                              lr=7
                         k=7 -------  z=0.857
    k=4 -----  z=0.75         lr=6
                         k=6 -------  z=0.714
                              lr=5
                         k=5 -------  z=0.653
    k=3 -----  z=0.5


  """
  nintH = len(zmH)
  nintM = len(zmM)
 
  zzH  = -abs(zzH)
  zmH  = -abs(zmH)
  zzM  = -abs(zzM)
  zmM  = -abs(zmM)
  dhH  = abs(np.diff(zzH))
  hbtH = np.min(zzH)
  hbtM = np.min(zzM)
#
# To avoid errors in the bottom cells,
# make HYCOM layers at least the same depth as MOM:
  zzH[-1] = np.min([hbtH, hbtM-0.01])
 
  Htot     = np.sum(dhH*H1d)
  Mtot     = 0.
  dz_tot   = 0.
  M1d      = np.zeros((nintM))
  for k in range(nintM):
#    print(f"k={k}")
    eta_tot  = 0.
    ztop = zzM[k]
    zbtm = zzM[k+1]
    dzM  = abs(ztop - zbtm)

    if dzM < eps_dZ:
      if fill_btm:
        M1d[k] = fill_val
      else:
        M1d[k] = M1d[k-1]
      continue

# Find HYCOM layer interfaces for MOM layer k bounded by:
# MOM bottom interface
    iHbtm_btm = np.argwhere(zzH <= zbtm)[0][0]   #  btm HYCOM intrf for MOM btm
    iHbtm_top = iHbtm_btm - 1                    # top HYCOM intrf for MOM btm
    zHbtm_btm = zzH[iHbtm_btm]                   # depth of HYCOM btm intrf
    zHbtm_top = zzH[iHbtm_top]                   # depth of HYCOM srf intrf

# MOM top interface:
    iHtop_btm = np.argwhere(zzH < ztop)[0][0]   # btm HYCOM intrf for MOM top
    if iHtop_btm == 0:
      iHtop_top = 0
    else:
      iHtop_top = iHtop_btm - 1                    # top HYCOM for cice6 top
    zHtop_btm = zzH[iHtop_btm]
    zHtop_top = zzH[iHtop_top]
  
# Find weights:
# eta1 = thickn from zbtm MOM to HYCOM top interface
# if MOM lr is inside 1 HYCOM lr eta1 = dhM (MOM lr thickness)
    eta1   = min([zHbtm_top,ztop]) - zbtm
# eta2 = thickn from ztop MOM to zbtm HYCOM
# Check if MOM ice layer is within cice H layer
    if iHbtm_btm == iHtop_btm and iHbtm_top == iHtop_top:
      eta2   = 0.
    else:
      eta2 = ztop - max([zHtop_btm,zbtm])

# Compute value in MOM layer:
    kH_btm = iHbtm_top    # HYCOM layer # where MOM bottom interface
    kH_top = iHtop_top    # HYCOM layer # where MOM top interf
    qvM = (H1d[kH_btm]*eta1 + H1d[kH_top]*eta2)/dzM

    M1d[k]   = qvM
    Mtot     = Mtot + qvM*dzM
    eta_tot  = eta_tot + eta1 + eta2
    dz_tot   = dz_tot + dzM
    rdz_eta  = abs(1.-eta_tot/dzM)

#    print(f"eta_tot={eta_tot}, dzM={dzM}, eta/dzM={rdz_eta}")
#  print(f"Input tot={Htot}, Output tot={Mtot}, dz_tot={dz_tot}")
  return M1d


