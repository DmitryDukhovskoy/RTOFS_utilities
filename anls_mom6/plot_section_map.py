"""
  Show section on the map
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pdb
import importlib
import time
#import struct
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/mom6_utils')

import mod_read_hycom as mhycom
import mod_time as mtime
from mod_utils_fig import bottom_text

import mod_misc1 as mmisc
import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)

plt.ion()

# Check orientation of the line/norms 
# curve_ornt
f_cont = False     # load saved and start from last saved
#sctnm = 'DavisS2'
#sctnm = 'Yucatan'  # 2-leg right angle section
#sctnm = 'Yucatan2'  # slented section
#sctnm = 'DavisS2'
#sctnm = 'Fram79s2'
#sctnm = 'BarentsS'
sctnm = 'BeringS'
#sctnm = 'DenmarkS'
#sctnm = 'IclShtl'
#sctnm = 'ShtlScot'
#sctnm = 'LaManch'
#sctnm = 'NAtl39'
#======= Ocean Sections =====
#sctnm = 'BaffNAFram'
#sctnm = 'AlaskaIcld' 
#sctnm = 'GoMCarib'

pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
DX, DY = mhycom.dx_dy(LON,LAT)
idm  = HH.shape[1]
jdm  = HH.shape[0]

# Either strait or ocean section:
NPsect = False
STR = mom6vld.ocean_straits()
if sctnm in STR:
  oc_strait = True

  nlegs = STR[sctnm]["nlegs"]
  I1    = STR[sctnm]["xl1"]
  I2    = STR[sctnm]["xl2"]
  J1    = STR[sctnm]["yl1"]
  J2    = STR[sctnm]["yl2"]
  Ni    = np.zeros((nlegs))
  Nj    = np.zeros((nlegs))
  IJ    = np.zeros((nlegs+1,2))
  for kk in range(nlegs):
    Ni[kk]     = I2[kk]+1
    Nj[kk]     = J2[kk]+1
    IJ[kk,0]   = I1[kk]
    IJ[kk,1]   = J1[kk]
    IJ[kk+1,0] = I2[kk]
    IJ[kk+1,1] = J2[kk]

else:
  STR = mom6vld.ocean_sections()
  if sctnm in STR:
    oc_strait = False
    NPsct   = STR[sctnm]["NP"]
    Is      = STR[sctnm]["II"]
    Js      = STR[sctnm]["JJ"]
    Is      = np.array(Is)
    Js      = np.array(Js)
# Replace with N. Pole values if needed:
# MOM6 and RTOFS have differen # of jdm
    if NPsct:
      dJ = abs(Js - jdm)
      Js = np.where(dJ < 3, jdm-1, Js)

    nlegs   = len(Is) - 1
    IJ      = np.zeros((nlegs+1,2))
    IJ[:,0] = Is
    IJ[:,1] = Js

  else:
    raise Exception(f"Name {sctnm} is not defined as a strait or section")


SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive')
#    SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='negative')
II     = SGMT.I_indx
JJ     = SGMT.J_indx
hLsgm1 = SGMT.half_Lsgm1
hLsgm2 = SGMT.half_Lsgm2
Vnrm1  = SGMT.half_norm1
Vnrm2  = SGMT.half_norm2
nLeg   = SGMT.Leg_number
XX     = LON[JJ,II]
YY     = LAT[JJ,II]
Hb     = HH[JJ,II]
LSgm   = np.zeros((len(II)))  # total segment length = half1 + half2
for ik in range(len(II)):
   LSgm[ik] = hLsgm1[ik] + hLsgm2[ik]

Hb    = HH[JJ,II].squeeze()


IIhf = []
JJhf = []
nI = len(Hb)
XI = np.arange(0, nI, 1, dtype=int)  # index # along the seciont
XI = (np.cumsum(LSgm)-LSgm[0])*1.e-3           # distance along section, km

if NPsect:
  ax1 = mom6vld.plot_section_orthomap(II, JJ, IJ, LON, LAT, HH, \
                lon0=-10, lat0=60, XI=XI, sttl=sctnm, btx='plot_section_map.py')
else:
  ax1 = mom6vld.plot_section_map(II, JJ, IJ, Vnrm1, Vnrm2, IIhf, JJhf, \
             HH, fgnmb=1, XI=XI, dX=500., sttl=sctnm, btx='plot_section_map.py')

  f_lon = True
  if f_lon:
    ax1.contour(LON,list(range(-180,180,5)), 
                colors=[(0.9,0.9,0.9)], 
                linestyles='solid')
    ax1.contour(LAT,list(range(-80,89,2)), 
                colors=[(0.9,0.9,0.9)], 
                linestyles='solid')
  





