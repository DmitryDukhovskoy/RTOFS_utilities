"""
  Plot T/S profiles at a given locations
  fomr GDEM and WOA18 climatologies

  Cara Manning
  Baffin Bay T/S anomaly in central Baffing 2017-2019

The station locations are:
BB2: 72.75 N, 67 W
224: 70.43 N, 62.96 W

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
from copy import copy

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
import mod_read_ncoda as rncoda
from mod_utils_fig import bottom_text

MM = 8
latPrf = 72.75
lonPrf = -67.

pthgdem = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/ncoda/fix/gdem/'
pthwoa  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/WOA18/'
flgdem  = f"gdem_salt{MM:02d}.short"
fld_gdem = os.path.join(pthgdem, flgdem)

import mod_gdem as mgdem
Zgdem, Tgdem = mgdem.gdem_profile(pthgdem, MM, 'temp', lonPrf, latPrf)
Zgdem, Sgdem = mgdem.gdem_profile(pthgdem, MM, 'salt', lonPrf, latPrf)

import plot_sect as psct
importlib.reload(psct)

fgnmb=1
plt.ion()
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.35, 0.8])
ax2 = plt.axes([0.55, 0.1, 0.35, 0.8])

stl = f"GDEM T MM={MM} Lon={lonPrf:.2f} Lat={latPrf:.2f}"
psct.plot_profile(ax1, Tgdem, Zgdem, 10., stl=stl, zlim=-2100, \
                  xl1=-2, xl2=5.5)

stl = f"GDEM S MM={MM} Lon={lonPrf:.2f} Lat={latPrf:.2f}"
psct.plot_profile(ax2, Sgdem, Zgdem, 10., stl=stl, zlim=-2100, \
                  xl1=29.5, xl2=34.95)

btx = 'plot_TSprof_gdem.py'
bottom_text(btx)
