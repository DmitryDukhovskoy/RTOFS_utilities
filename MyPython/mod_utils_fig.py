"""
  Functions for plotting and figure settings
"""
import os
import numpy as np
import matplotlib.pyplot as plt

def bottom_text(btx, gtdir='RTOFS_utilities', \
                ipwd=1, pos=[0.05, 0.05], fsz=8, f_short=False):
  # add machine name:
  syst_info = os.uname()
  machine = syst_info.nodename

  #sgit = 'github.com/DmitryDukhovskoy/' + gtdir
  drr = os.getcwd()
  bnm = os.path.basename(drr)
  if f_short:
    btnm = os.path.join(bnm, f"@{machine}: {btx}")
  else:
    btnm = os.path.join(drr, f"@{machine}: {btx}")

  if ipwd==0:
    btnm = (f"@{machine}: {btx}")

  plt.text(pos[0],pos[1],btnm,horizontalalignment='left',
         transform=plt.gcf().transFigure, fontsize=fsz)
  return

def colorbar_horiz(fig1, ax1, img, rmin=None, rmax=None, 
        cb_height=0.016, cb_gap=0.016, nticks=11, fontsz=12, 
        decim=2, extd='both'):
  """
    Plot colorbar under the bottom axis = ax1
    cb_height - cbar height
    cb_gap - gap between cbar and bottom axis
    nticks - number of ticks on cbar
    fontsz - fontsize on cbar
    decim - number of decimals in cbar tick labels: 0, 1, 2, 3, ...
    extd - draw an arrow at the end: min, max, both, neither
  """
  from matplotlib.ticker import StrMethodFormatter

  fig1.canvas.draw()

  pos = ax1.get_position()

  axcb = fig1.add_axes([
    pos.x0,
    pos.y0 - cb_gap - cb_height,
    pos.x1 - pos.x0,
    cb_height
  ])

  clb = fig1.colorbar(
    img,
    cax=axcb,
    extend=extd,
    orientation='horizontal'
  )

  if rmin is None:
    rmin = img.norm.vmin
  if rmax is None:
    rmax = img.norm.vmax

  ticks = np.linspace(rmin, rmax, nticks)
  clb.set_ticks(ticks)

  clb.ax.tick_params(
    direction='in',
    length=10,
    labelsize=fontsz
  )

  fmt = f'{{x:.{decim}f}}'
  clb.ax.xaxis.set_major_formatter(StrMethodFormatter(fmt))

  return clb

def correct_evensp_grid(x, y):
  xcrct = x.copy()
  ycrct = y.copy()
  if x.ndim == 1:
    pass
  elif x.ndim == 2:
    x_row = x[0, :]
    if not np.allclose(x_row, x):
      raise ValueError("The rows of 'x' must be equal")
  else:
    raise ValueError("'x' can have at maximum 2 dimensions")
#   
  if y.ndim == 1:
    pass
  elif y.ndim == 2:
    y_col = y[:, 0]
    if not np.allclose(y_col, y.T):
      raise ValueError("The columns of 'y' must be equal")
  else:
    raise ValueError("'y' can have at maximum 2 dimensions")

  nx = len(x_row)
  ny = len(y_col)

  width = x_row[-1] - x_row[0]
  height = y_col[-1] - y_col[0]

  rtol = 1.e-6
  atol = 0.0
  dx = width/(nx-1)
  if not np.allclose(np.diff(x_row), dx, atol, rtol):
    x1r = np.arange(x_row[0],x_row[-1]+0.1*dx,dx)

    for ii in range(ny):
      xcrct[ii] = x1r

# breakpoint() 
  dy = height/(ny-1)
  if not np.allclose(np.diff(y_col), dy, atol, rtol):
    y1r = np.arange(y_col[0],y_col[-1]+0.1*dy,dy)

    for ii in range(nx):
      ycrct[:,ii] = y1r

# breakpoint()
  x_row2 = xcrct[0, :]
  np.allclose(np.diff(x_row2), dx, atol, rtol)

  return xcrct, ycrct

def plot_fld2D(aa,ctitle,cl1,cl2,X=np.array([]),Y=np.array([]),clrmp='turbo',
               btx='mod_utils_fig.py', figout='plot_field.jpg'):
  """
    Create a quick figure of 2D map and 
    put it on web to view via browser
    otherwise python get slow
    open figure in the browser:
        https://www-dev.star1.nesdis.noaa.gov/ncei_wod/pthn_fig.jpg
  """
  plt.ion()
  cmpS = plt.cm.get_cmap(clrmp)
  fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  jdm=aa.shape[0]
  idm=aa.shape[1]

#  if np.isfinite(X) and np.isfinite(Y):
  if X.size == 0 | Y.size == 0:
    im = ax1.pcolormesh(X,Y,aa,cmap=cmpS,shading='flat')
    ax1.axis('scaled')
    ax1.set(xlim=(np.min(X),np.max(X)),
            ylim=(np.min(Y),np.max(Y)))
  else:
    im = ax1.pcolormesh(aa,cmap=cmpS,shading='flat')
    ax1.axis('scaled')
    ax1.set(xlim=(0,idm),ylim=(0,jdm))

  im.set_clim(cl1,cl2)
  ax1.autoscale(enable=True, axis='both', tight=True)

  ax1.set_title(ctitle)
#
# Colorbar
  divider = make_axes_locatable(ax1)
  cax = divider.append_axes("bottom",size="5%",pad=0.25)
  fig1.colorbar(im,ax=ax1,orientation='horizontal',cax=cax)


  bottom_text(btx)

  plt.show()
  plt.savefig(figout)
  plt.close(fig1)

  wwwpth='/net/www-dev/www/ncei_wod/'
  wwwfig='pthn_fig.jpg'
  cmd = ('mv {0} {1}{2}'.format(figout,wwwpth,wwwfig))
  os.system(cmd)


  return

def minmax_clrmap(dmm, pmin=5, pmax=95,
                  cpnt=0.01, fsym=False):

  a1 = np.nanpercentile(dmm, pmin)
  a2 = np.nanpercentile(dmm, pmax)

  rmin = cpnt * np.floor(a1 / cpnt)
  rmax = cpnt * np.ceil (a2 / cpnt)

  if fsym and (rmin < 0.) and (rmax > 0.):
    mx = max(abs(rmin), abs(rmax))
    rmin = -mx
    rmax =  mx

  return rmin, rmax


