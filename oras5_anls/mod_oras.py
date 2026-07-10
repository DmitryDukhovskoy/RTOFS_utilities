import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import mod_misc1 as mmisc

def demean_ssh(ssh):
  """
    Demean SSH in ORAS5 fields
  """
  iB1 = 200
  iB2 = 693
  jB1 = 860
  jB2 = 1020

  amn = np.nanmean(ssh[jB1:jB2+1,iB1:iB2+1])
  ssh_dmn = ssh-amn

  return ssh_dmn, amn

def calc_grad_sizeBG(ABG, lonh, lath, Acell, MSK_BG, dltC=0.005):
  """
    Find ssh grad in the BG and size BG based on the
    last closed contour
    ABG is a subset region that includes BG
  """
  plt.ioff()
  figA = plt.figure(10,figsize=(8,8))
  plt.clf()

  axA1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  # BG Domain excluding shallow regions where max SSH can > BG 
  Adeep = np.where(MSK_BG==0,np.nan,ABG) 
  #amx = np.nanmax(Adeep)
  ssh_max = np.nanmax(Adeep)
  jmax, imax = np.unravel_index(np.nanargmax(Adeep), Adeep.shape)
  #tol = 1e-6
  #IJNDX = np.argwhere(np.isclose(ABG, ssh_max, atol=tol))
  #if IJNDX.size > 0:
  #    jmax, imax = IJNDX[0]  # First approximate match
  #else:
  #    print(f"No value in ABG close to ssh_max={ssh_max} found")
  jdma, idma = ABG.shape
  IIA, JJA = np.meshgrid(np.arange(idma), np.arange(jdma))

  cff = 0
  if 1. <= dltC*10 < 10.:
    cff = 10
  elif 1. <= dltC*100 < 10.:
    cff = 100
  elif 1. <= dltC*1000 < 10.:
    cff = 1000
  else:
    raise Exception(f'Could not determine coefficient for {dltC}')

  Cstart = np.floor(ssh_max*cff)/float(cff)

  # Mask all land:
  ABG = np.where(np.isnan(ABG),-0.1,ABG)
  cc0 = Cstart
  CNTR = []
  SSHval = []
  while cc0 >= 0:
    CS = axA1.contour(ABG,[cc0])
    SGS  = CS.allsegs[0]  # should be only 1 contoured value
    nsgs = len(SGS)
#    print(f'Contour {cc0:.3f}')
  # Check all contours that includ ssh_max and check if they are closed
    for isg in range(nsgs):
      XY = SGS[isg]
      X  = XY[:,0]
      Y  = XY[:,1]

      if not mmisc.inpolygon_1pnt(imax,jmax, X,Y):
        continue
      #if len(X) <= nmin: continue
    
      dEnd = np.sqrt((X[0]-X[-1])**2+(Y[0]-Y[-1])**2)
      if dEnd < 1.:
        CNTR.append(XY)
        SSHval.append(cc0)
        continue

    cc0 -= dltC

  if len(CNTR) == 0: 
    raise Exception("Could not find closed SSH contours")

  ssh_min = SSHval[-1]  # last closed contour in the BG
  dltH = ssh_max - ssh_min
  xymin = CNTR[-1]
  xx = xymin[:,0]
  yy = xymin[:,1]

  di = 5
  xindx = xx[0::5].astype(int)
  yindx = yy[0::5].astype(int)
  XC = lonh[yindx,xindx]
  YC = lath[yindx,xindx]
  xmax = lonh[jmax,imax]
  ymax = lath[jmax,imax]
  DST = mmisc.dist_sphcrd(YC,XC,ymax,xmax)
  assert(np.min(DST)>0.), f'min DST in grad_SSH == 0'
  grdH = dltH / DST

  # Find BG area:
  MS, _, _ = mmisc.inpolygon_v2(IIA, JJA, xindx, yindx)  # 
  JBG, IBG = np.where(MS == 1) 
  area_bg = np.sum(Acell[JBG,IBG])

  grdH_mn = np.mean(grdH)
  grdH_md = np.median(grdH)
  grdH_lprc = np.percentile(grdH,10)
  grdH_uprc = np.percentile(grdH,90)
  

  f_show = False
  if f_show:
    import mod_colormaps as mclrmps
    clrmp = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBuGn_r')
    rmin = -0.5
    rmax = 0.5
    axA1.pcolormesh(ABG,cmap=clrmp, vmin=rmin, vmax=rmax)
    axA1.plot(xx,yy,'r-')
    axA1.plot(imax,jmax,'o')
    axA1.axis('scaled')

  plt.close(figA)

  plt.ion()

  return grdH_mn, grdH_md, grdH_lprc, grdH_uprc, area_bg, ssh_max





