"""
  Functions for Gulfstream analysis
"""
import os
import numpy as np
import sys
import importlib
import matplotlib.pyplot as plt
import datetime

def derive_contour(AA, tz0=12., xFS=2565, yFS=1809):
  """
    AA is T field with regions that are outside the 
   study region being masked out
   index space - for 0.08 global RTOFS tripolar grid 4500 x 3298 
  """

  plt.ion()
  figA = plt.figure(10,figsize=(8,8))
  plt.clf()

  axA1 = plt.axes([0.1, 0.2, 0.7, 0.7])
  CS   = axA1.contour(AA,[tz0])

  axA1.axis('equal')
  axA1.set_xlim([2520,3050])
  axA1.set_ylim([1800,2250])

  SGS  = CS.allsegs[0]  # should be only 1 contoured value
  nsgs = len(SGS)

# Delete all closed contours
  CNTR = []  
  for isg in range(nsgs):
    XY = SGS[isg]
    X  = XY[:,0]
    Y  = XY[:,1]

    dEnd = np.sqrt((X[0]-X[-1])**2+(Y[0]-Y[-1])**2)
    if dEnd < 1.:
      continue

    CNTR.append(XY)

# Aranage all segments in order
# First find segment that starts in Fl. Str
  nC  = len(CNTR)
  TCNT = []

  for ii in range(nC):
    if ii == 0:
      x0 = xFS
      y0 = yFS
    else:
      x0 = TCNT[-1,0]
      y0 = TCNT[-1,1]

    xsgm, ysgm, imin, jmin = arange_1segm(CNTR,x0,y0)  

# Remove selected segment:
    CNTR.pop(imin)

    if ii == 0:
      TCNT = np.transpose(np.array((xsgm,ysgm)))
    else:
      aa   = np.transpose(np.array((xsgm,ysgm)))
      TCNT = np.append(TCNT, aa, axis=0)


#X = TCNT[:,0]
#Y = TCNT[:,1]
#axA1.plot(X,Y,'.-') 

  return TCNT



def arange_1segm(CNTR,x0,y0):
  """
    Find segment closest to x0, y0 
    arange the orientation of the segment
    so that it starts from the pnt closest to x0, y0
    CNTR - list with segments X,Y as np arrays
  """
  nC  = len(CNTR)
  DFS = np.zeros((nC,2))*np.nan
  for isg in range(nC):
    XY = CNTR[isg]
    X  = XY[:,0]
    Y  = XY[:,1]

    d1 = np.sqrt((X[0]-x0)**2+(Y[0]-y0)**2)
    d2 = np.sqrt((X[-1]-x0)**2+(Y[-1]-y0)**2)

    DFS[isg,0] = d1
    DFS[isg,1] = d2   
 
  imin = np.argmin(np.min(DFS, axis=1))
  jmin = np.argmin(np.min(DFS, axis=0))

  xsgm = CNTR[imin][:,0]
  ysgm = CNTR[imin][:,1]

  if jmin == 1:
    xsgm = np.flip(xsgm)
    ysgm = np.flip(ysgm)

  return xsgm, ysgm, imin, jmin





