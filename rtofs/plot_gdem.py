"""
  Plot GDEM from binary 
  see code:
  /scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/ncoda/sorc/ncoda_qc/libsrc/prfobs
  gdem_mod.f
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text

import mod_utils as mutil
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_read_ncoda as rncoda
#importlib.reload(rncoda)
import mod_gdem as mgdem
importlib.reload(mgdem)

imo = 4
fld = 'salt' # salt or temp

pthgdem = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/ncoda/fix/gdem/'

LON, LAT, ZM = mgdem.gdem_grid()
NI = LON.shape[0]
NJ = LAT.shape[0]
NK = ZM.shape[0]

A3d = mgdem.read_gdem3d(pthgdem, imo, fld, NI, NJ, NK)





