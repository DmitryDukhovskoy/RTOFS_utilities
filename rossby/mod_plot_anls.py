"""
  Functions for plotting/analyzing results
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys

sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
import mod_utils_fig
from mod_utils_fig import bottom_text
from mod_utils_fig import plot_fld2D
plt.ion()

btx2='mod_plot_anls.py'
def zonalavrg(BB,ctitle,lat,LMsk,btx=btx2,ifg=1):
  idm = BB.shape[1]
  jdm = BB.shape[0]
  
# Zonal averaging:
  Bav = np.zeros(jdm)
  for jj in range(jdm):
    aa = BB[jj,:].squeeze()
    Bav[jj] = np.nanmean(aa)

  fig1 = plt.figure(ifg,figsize=(8,8), constrained_layout=False)
  plt.clf()
  ax1 = plt.axes([0.1,0.4,0.8,0.5])
  ax1.plot(lat,Bav)
  Bmx=np.nanmax(Bav)
  ax1.set(xlim=(-89,89),ylim=(0,1.1*Bmx))

# Make a plot with major ticks that are multiples of 20 
# and minor ticks that
# are multiples of 5.  Label major ticks with '%d' 
# formatting but don't label
# minor ticks.
  ax1.xaxis.set_major_locator(MultipleLocator(20))
  ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  ax1.yaxis.set_major_locator(MultipleLocator(50))
  ax1.yaxis.set_major_formatter(FormatStrFormatter('%d'))

# For the minor ticks, use no labels; default NullFormatter.
  ax1.xaxis.set_minor_locator(MultipleLocator(5))
  ax1.yaxis.set_minor_locator(MultipleLocator(10))

  plt.grid(axis='both', linestyle='--')
  ax1.set_title(ctitle)

  bottom_text(btx,pos=[0.1,0.3]) 

  plt.show()

  return


