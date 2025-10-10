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

PPTHN = '/home/Dmitry.Dukhovskoy/python'
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

import mod_read_hycom as mhycom
import mod_time as mtime
from mod_utils_fig import bottom_text

import mod_misc1 as mmisc
import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)

from mpl_toolkits.basemap import Basemap, cm
import matplotlib.colors as colors
import matplotlib.mlab as mlab


plt.ion()

SCTNMS = ['DavisS2','Fram79s2','BarentsS','BeringS','DenmarkS','IclShtl',
          'ShtlScot','LaManch','NAtl39','BaffNAFram','AlaskaIcld','GoMCarib',
          'Yucatan2','FlorCabl']

pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
DX, DY = mhycom.dx_dy(LON,LAT)
idm  = HH.shape[1]
jdm  = HH.shape[0]

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

nscts = len(SCTNMS)
for isct in range(nscts):
# Either strait or ocean section:
  sctnm = SCTNMS[isct]
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

  if isct==0:
    m = Basemap(projection='ortho', lon_0=-50, lat_0=45, resolution='l')

    xh, yh = m(LON,LAT)  # modl grid coordinates on the projections

    m.drawcoastlines()
    m.fillcontinents(color=(0.2,0.2,0.2))
    m.drawparallels(np.arange(-90.,120.,10.))
    m.drawmeridians(np.arange(-180.,180.,10.))
    m.contour(xh, yh, HH, [-4000,-3000,-2000,-1000],
              colors=[(0.8,0.8,0.8)],
              linestyles='solid')


    sttl='Sections'
    ax1.set_title(sttl)

    f_lon = True
    if f_lon:
      ax1.contour(LON,list(range(-180,180,5)), 
                  colors=[(0.9,0.9,0.9)], 
                  linestyles='solid')
      ax1.contour(LAT,list(range(-80,89,2)), 
                  colors=[(0.9,0.9,0.9)], 
                  linestyles='solid')
    
  m.plot(xh[JJ,II],yh[JJ,II],'.')

btx='plot_allsections_orthomap.py'
bottom_text(btx)




