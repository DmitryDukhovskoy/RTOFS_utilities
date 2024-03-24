"""
  Functions for plotting and figure settings
"""

import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.pyplot as plt
from copy import copy
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_prof(Z,T,ctl="",fgn=1):
  plt.ion()
  fig1 = plt.figure(fgn,figsize=(7,8), constrained_layout=False)
  fig1.clf()
  plt.plot(T,Z)
  plt.title(ctl)

  return

def plot_2prof(Z1,T1,Z2,T2,ctl1="",ctl2="",fgn=1):
  plt.ion()
  fig1 = plt.figure(fgn,figsize=(7,8), constrained_layout=False)
  fig1.clf() 
  ax1 = plt.axes([0.1, 0.1, 0.35, 0.8])
  plt.plot(T1,Z1)
  ax1.grid(True)
#  ax1.set_xlabel('cycle/hr')
  ax1.set_title(ctl1)

  ax2 = plt.axes([0.55, 0.1, 0.35, 0.8])
  ax2.plot(T2,Z2)
  ax2.grid(True)
  ax2.set_title(ctl2)

def plot_checkN2Phi(N2,ZZ,Z_phi,Phi,PRB,zbtm,prbnm,krcrd,Rrsb,Cphase,fgn=4):
  """
    Plot N2, N2 smoothed, 1st eigenvector
  """
  import numpy as np
  import mod_utils_rossby as ursb
  import mod_solver as uslv

  N2zF = uslv.runmn(N2,ZZ,mnwnd=3)  # Filter
  N_hrz = np.sqrt(N2)
  N_chr = 3600.*N_hrz   # cycles per hr
  NF_hrz = np.sqrt(N2zF)
  NF_chr = 3600.*NF_hrz

  plt.ion()
  fig4 = plt.figure(fgn,figsize=(8,8), constrained_layout=False)
  fig4.clf()
  ax1 = plt.axes([0.1, 0.1, 0.35, 0.8])
  plt.plot(N_chr,Z_phi)
  plt.plot(NF_chr,Z_phi)
  ax1.grid(True)
  ax1.set_xlabel('cycle/hr')
  uniqnmb = PRB.uniqnmb
  ctl = 'N^2, {1}, U#{2} Zbtm={0:.1f}'.format(zbtm,prbnm,uniqnmb)
  ax1.set_title(ctl)

  ax2 = plt.axes([0.55, 0.1, 0.35, 0.8])
  ax2.plot(Phi,Z_phi)
  ax2.grid(True)
  ctl2 = '$\Phi$(z)'
  ax2.set_title(r'$\Phi$(z)')

  CC  = PRB.CC
  lat = PRB.lat
  lon = PRB.lon
  yr  = PRB.yr
  mo  = PRB.mo
  md  = PRB.md
  gtm = PRB.gtm
  cid = PRB.cid

  Txt = ['Cruise# {1}{0}'.format(cid,CC)]
  Txt.append('lat={0:.2f}N, lon={0:.2f}E'.format(lat,lon))
  Txt.append('{0:02d}/{1:02d}/{2}:{3:.2f}'.format(md,mo,yr,gtm))
  Txt.append('UnqNmb={0}'.format(uniqnmb))
  Txt.append('Seq.Recrd={0}'.format(krcrd))
  Txt.append('Rossby(km)={0:.1f}'.format(Rrsb))
  Txt.append('PhaseSp(m/s)={0:.2f}'.format(Cphase))
  nt=len(Txt)

  x0 = 1
  y0 = 1
  ax3 = plt.axes([0.56,0.12,0.25,0.1])
  ax3.set(xlim=(x0,4.5),ylim=(y0,6))
  for itt in range(nt):
    ssl=Txt[itt]
    x1 = x0
    y1 = y0+itt
    ax3.text(x1,y1,ssl)
  ax3.axis('off')

  return 

def bottom_text(btx, pos=[0.03, 0.05], fsz=10):
  drr = os.getcwd()
  btnm = drr+'/'+btx
  plt.text(pos[0],pos[1],btnm,horizontalalignment='left',
         transform=plt.gcf().transFigure, fontsize=fsz)


def plot_fld2D(aa,ctitle,cl1,cl2,X=np.array([]),Y=np.array([]),\
               clrmp='turbo', fig_save=False, \
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

  if fig_save:
    plt.show()
    plt.savefig(figout)
    plt.close(fig1)

    wwwpth='/net/www-dev/www/ncei_wod/'
    wwwfig='pthn_fig.jpg'
    cmd = ('mv {0} {1}{2}'.format(figout,wwwpth,wwwfig))
    os.system(cmd)


  return

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


