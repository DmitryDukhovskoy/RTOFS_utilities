"""
  Plot T/S sections across straits and ocean sections
  from GDEM climatology

  Note: The original GDEM data on 78 z-levels is in-situ temperature. 
 
  HYCOM - T is potential (wrt surface)
  MOM6 - T is potential (?)
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
#from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/mom6_utils')

import mod_time as mtime
from mod_utils_fig import bottom_text

import mod_mom6_valid as mom6vld
importlib.reload(mom6vld)

sctnm = 'Fram79s2'
#sctnm = 'DavisS2'
#sctnm = 'Yucatan2'  # slanted section
#sctnm = 'BarentsS'
#sctnm = 'BeringS'
#sctnm = 'DenmarkS'
#sctnm = 'IclShtl'
#sctnm = 'ShtlScot'
#sctnm = 'LaManch'
#sctnm = 'NAtl39'
#======= Ocean Sections =====
#sctnm = 'BaffNAFram'
#sctnm = 'AlaskaIcld' 

#fld2d = 'salt'
fld2d = 'temp'
f_pot = True    # true - convert to potential wrt to P=0 from tn situ
moS   = 1
moE   = 12

if fld2d == 'temp' and f_pot:
  f_conv = True
else:
  f_conv = False

import mod_misc1 as mmisc
import mod_mom6 as mom6util
import mod_colormaps as mcmp
import mod_read_hycom as mhycom
import mod_gdem as mgdem
importlib.reload(mom6util)
importlib.reload(mmisc)
importlib.reload(mcmp)
importlib.reload(mgdem)

pthgdem = '/scratch2/NCEPDEV/marine/Zulema.Garraffo/ncoda/fix/gdem/'
pthout  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/GDEM/'

LON, LAT, ZM = mgdem.gdem_grid()
NI = LON.shape[0]
NJ = LAT.shape[0]
NK = ZM.shape[0]

LON = np.where(LON > 180., LON-360., LON)


# Load GDEM land mask:
Lmsk = mgdem.read_gdem_Lmask()

# Section indices and coordinates:
# if do not exist - extract using MOM6 section info
fsect  = f"indx_coord_{sctnm}.py"
ffsect = os.path.join(pthout, fsect)

try: 
  with open(ffsect, 'rb') as fid:
    [II, JJ, XX, YY, Lsgm, Hbtm] = pickle.load(fid)
except:
# Load section info from MOM6:
  dnmb1  = mtime.datenum([2021,1,1])
  dnmb2  = mtime.datenum([2021,12,31])
  dv1    = mtime.datevec(dnmb1)
  dv2    = mtime.datevec(dnmb2)
  pthmom = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/' + \
            'MOM6_CICE6/expt003/'
  flmom  = f"mom6-003_saltVFlx_{dv1[0]}" + \
            f"{dv1[1]:02d}-{dv2[0]}{dv2[1]:02d}_{sctnm}.pkl"
  ffmom = pthmom + flmom
  print('Loading ' + ffmom)

  with open(ffmom, 'rb') as fid:
    F2D = pickle.load(fid)
  #TM   = F2D.TM
  #A2dT = F2D.Fld2D  # 2D fields: Time x depth x Width
#  IM   = F2D.Iindx
#  JM   = F2D.Jindx
  XM   = F2D.LON
  YM   = F2D.LAT
  #Lsgm = F2D.Lsgm
  HbtmM = F2D.Hbtm
  #ZZi  = F2D.ZZi

  # Find section indices on GDEM grid:
  II, JJ, XX, YY, Lsgm, Hbtm = mgdem.find_gdem_indx(XM, YM, LON, LAT, Hb0=HbtmM)
  with open(ffsect, 'wb') as fid:
    pickle.dump([II, JJ, XX, YY, Lsgm, Hbtm], fid)
  
# Read GDEM:
ncc = 0
for imo in range(moS, moE+1):
  A3d = mgdem.read_gdem3d(pthgdem, imo, fld2d, NI, NJ, NK)

  if ncc == 0:
    A2d = np.squeeze(A3d[:,JJ,II])
  else:
    A2d = A2d + np.squeeze(A3d[:,JJ,II])

  if f_conv:
    S3d = mgdem.read_gdem3d(pthgdem, imo, 'salt', NI, NJ, NK)
    if ncc == 0:
      S2d = np.squeeze(S3d[:,JJ,II])
    else:
      S2d = S2d + np.squeeze(S3d[:,JJ,II])
 
  ncc += 1
 
A2d  = A2d/ncc 
if f_conv:
  import mod_swstate as msw
  
  S2d  = S2d/ncc
# for depth loop
# convert depths to pressure dbar
# compute Tpot
  TPOT = A2d.copy()
  nk = ZM.shape[0]
  for ik in range(nk):
    zm0 = np.zeros((len(YY))) + ZM[ik]
    pr0_db, pr0_pa = msw.sw_press(zm0, YY)
    pr_ref = 0.    # reference pressure
    tpot   = msw.sw_ptmp(S2d[ik,:], A2d[ik,:], pr0_db, pr_ref)
    TPOT[ik,:] = tpot

  A2d = TPOT.copy()

A2di = mom6vld.box_fltr(A2d, npnts=3)

STR = mom6vld.ocean_straits()
if sctnm in STR:
  oc_strait = True
else:
  STR = mom6vld.ocean_sections()
  if sctnm in STR:
    oc_strait = False
  else:
    raise Exception(f"Name {sctnm} is not defined as a strait or section")

IJ = np.zeros((len(II),2))
IJ[:,0] = II
IJ[:,1] = JJ

btx = 'plot_TSxsect_gdem.py'

import mod_utils as mutil
import plot_sect as psct
#importlib.reload(mcmp)

if fld2d == 'salt':
#  cmpr = mcmp.colormap_salin(clr_ramp=[0.94,0.95,1])
  cmpr = mcmp.colormap_haline()
  rmin = STR[sctnm]["smin"]
  rmax = STR[sctnm]["smax"]
  cntr1 = STR[sctnm]["scntr"]
elif fld2d == 'temp':
#  cmpr = mcmp.colormap_temp()
  cmpr = mcmp.colormap_temp2()
  rmin = STR[sctnm]["tmin"]
  rmax = STR[sctnm]["tmax"]
  cntr1 = STR[sctnm]["tcntr"]

if fld2d == 'temp':
  if f_pot:
    addf = 'potent'
  else:
    addf = 'in situ'
else:
  addf = ' '

stl = f"GDEMv3 {sctnm} {addf} {fld2d} months: {moS}-{moE}"

HH = Lmsk.copy()
HH = np.where(HH<1., -0.1, HH)

plt.ion()
ni = len(II)
XI = np.arange(0, ni, 1, dtype=int)
XI = (np.cumsum(Lsgm)-Lsgm[0])*1.e-3 # distance along section, km
mom6vld.plot_xsect(XI, Hbtm, ZM, A2di, HH, stl=stl,\
                   rmin = rmin, rmax = rmax, clrmp=cmpr,\
                   IJs=IJ, btx=btx, btm_midpnt=True, cntr2=cntr1)



f_plt = False
if f_plt:
  plt.ion()
  fig1 = plt.figure(1,figsize=(9,9))

  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  ax1.plot(XM,YM,'.-')
  #ax1.contour(HH, [0.0], colors=[(0,0,0)], linewidths=1)
  ax1.axis('scaled')
  #ax1.set_xlim([xlim1,xlim2])
  #ax1.set_ylim([ylim1,ylim2])











