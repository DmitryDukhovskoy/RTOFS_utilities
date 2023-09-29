"""
  Draw vectors, arrows, compass diagrams etc
"""
import matplotlib.pyplot as plt
import numpy as np

def compass(u, v, nf=10, v_col=[0.,0.,0.], lwd=1, cf_ahd=0.1, beta=25.,
            larr_min=None, larr_max=None, ax=None):
# arrowprops = dict(color='darkorange', linewidth=2)
  """
  Program similar to Matlab compass
  draw a vecotr in polar coordinates
  using example from internet
  """

# Derive polar coordinates from U,V
  tht, Rr = cart2polar(u, v)

# breakpoint()
  if ax is not None:
#   plt.figure(nf)
    Ylim = np.max(ax.get_ylim())
  else:
    plt.figure(nf).clf()
    ax  = plt.axes(projection='polar')
    Ylim = -1.

  if larr_min is None:
    larr_min = 0.1*Rr
    larr_max = 0.2*Rr

  vec_props={
    'v_col': v_col,
    'lwd': lwd,
    'cf_ahd': cf_ahd,
    'beta': beta,
    'larr_min': larr_min,
    'larr_max': larr_max
  }

  draw_arrowF_polar(tht,Rr,ax=ax,**vec_props)

  Ylim = max(Ylim, 1.05*Rr)
  ax.set_ylim(0, Ylim)

  return ax

#
def cart2polar(x, y):
  """ 
  Convert cartisian coordinates --> polar (theta, radius)
  angle theta is in math convention: =0 on x axis, positive rotation - c/clckwise
  """

  RR = np.sqrt(x**2+y**2)
  THT = np.arctan2(y, x)
  RR = RR.tolist()
  THT = THT.tolist()

  return THT, RR
#
#
#
def draw_arrowF(x1,x2,y1,y2,**vec_props):
  """
  # draws arrow with filled "head"
  on a current figure/axes
  Note: first X-coordinates x1, x2, then Y-coordinates, Y1,y2, ...
  Draws an arrow between pnt (x1,y1) and pnt (x2,y2)
  Optional vector properties:
  cf_ahd - scaling coefficient of the arrowhead, if [0 to 1)
       arrow head is smaller than the vector
       otherwise, the vector will be closed by the arrowhead
  beta - angle between vector and arrow head beams (degrees)
  v_col - color ([R G B]) # default color is black
  lwd - line width, default = 1
  larr_min / larr_max - min/max size of the arrow heads
 
  """
  v_col    = vec_props.get('v_col',[0, 0, 0])
  lwd      = vec_props.get('lwd',1)
  cf_ahd   = vec_props.get('cf_vec',0.3)
  beta_dgr = vec_props.get('beta',25.)
  larr_min = vec_props.get('larr_min',0.)    # min arrow-head size
  larr_max = vec_props.get('larr_max',0.)

  xlim1, xlim2 = plt.gca().get_xlim()
  ylim1, ylim2 = plt.gca().get_ylim()

  uu = x2-x1
  vv = y2-y1
  sp = np.sqrt(uu*uu+vv*vv)
  alfa = np.arctan2(uu, vv)  # vector angle from Y
  beta = beta_dgr*np.pi/180.
  var = cf_ahd*sp
  if var < larr_min:
    var = larr_min

  if larr_max > larr_min and var > larr_max:
    var = larr_max

  dX2 = var*np.sin(alfa-beta)     # arrow head coordinates
  dX3 = var*np.sin(alfa+beta)
  dY2 = var*np.cos(alfa-beta)
  dY3 = var*np.cos(alfa+beta)
  dL = np.sqrt(dX2**2+dY2**2)

# Length of the vector with the arrow-head
  Lv=sp+dL-0.1*dL  # to avoid gap btw stem and arrowhead
  un=uu/sp
  vn=vv/sp
  x0=x1+un*Lv  # scale to adjust for the arrowhead
  y0=y1+vn*Lv

  ax2=x0-dX2
  ax3=x0-dX3
  ay2=y0-dY2
  ay3=y0-dY3
  plt.plot([x1, x2],[y1, y2], color=v_col, linewidth=lwd)       #vector

  X=[x0,ax2,ax3,x0]
  Y=[y0,ay2,ay3,y0]
# XY = list(zip(X,Y))
# XY = np.array(XY)
  plt.fill(X,Y, color=v_col)
  plt.axis('equal')
  plt.xlim([xlim1, xlim2])
  plt.ylim([ylim1, ylim2])

  return

def draw_arrowF_polar(tht,R,ax=None,**vec_props):
  """
  Plot an arrow in polar coordinates
  similar to compass in matlab
  tht - in radians
  """
  v_col    = vec_props.get('v_col',[0, 0, 0])
  lwd      = vec_props.get('lwd',1)
  cf_ahd   = vec_props.get('cf_vec',0.1)
  beta_dgr = vec_props.get('beta',25.)
  larr_min = vec_props.get('larr_min',0.)    # min arrow-head size
  larr_max = vec_props.get('larr_max',0.)

  xlim1, xlim2 = plt.gca().get_xlim()
  ylim1, ylim2 = plt.gca().get_ylim()

  beta = beta_dgr*np.pi/180.
  var = cf_ahd*R
  if var < larr_min:
    var = larr_min

  if larr_max > larr_min and var > larr_max:
    var = larr_max

#
# Find points for arrow head
  Vhead = R*cf_ahd
  dH = Vhead*np.tan(beta)
  Shft = R-Vhead  # arrow shaft
  Lp1 = np.sqrt(dH*dH+Shft*Shft)  # distance to arrow point 1
  dtht = np.arccos(Shft/Lp1)
  tht_p1 = tht-dtht
  tht_p2 = tht+dtht

  Pnts = [R,Lp1,Lp1]
  Angls = [tht,tht_p1,tht_p2]
  if ax is None:
    ax  = plt.axes(projection='polar')

  ax.plot([tht, tht],[0,R],color=v_col)
#   ax.plot(tht_p1,Lp1,'k*')
#   ax.plot(tht_p2,Lp1,'g*')
  ax.fill(Angls,Pnts, color=v_col)

  return


def compass_old(u, v, arrowprops=None, nf=1, ax=None):
# arrowprops = dict(color='darkorange', linewidth=2)
  """
  Program similar to Matlab compass
  draw a vecotr in polar coordinates
  using example from internet
  """

# Derive polar coordinates from U,V
  tht, Rr = cart2polar(u, v)

# breakpoint()
  if ax is not None:
    plt.figure(nf)
  else:
    plt.figure(nf).clf()
    ax  = plt.figure(nf).subplots(subplot_kw=dict(polar=True))

  kw = dict(arrowstyle="->", color='k')
  if arrowprops:
    kw.update(arrowprops)

  [ax.annotate("", xy=(tht, Rr), xytext=(0, 0),
               arrowprops=kw)]

# breakpoint()
#
# [ax.annotate("", xy=(tht, Rr), xytext=(0, 0),
#               arrowprops=kw) 
#              for tht,Rr in zip(THT, RR)]

  ax.set_ylim(0, 1.05*Rr)

  return ax

