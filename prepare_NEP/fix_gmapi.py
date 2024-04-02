# A quick fix in some of gmapi bouding box indices
# due to a bug in the previous version of find_gridpnts_box module
#  
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
#import struct
import datetime
import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
import yaml
from netCDF4 import Dataset as ncFile

#PPTHN = '/home/Dmitry.Dukhovskoy/python'
PPTHN = []
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mom6util
import mod_misc1 as mmisc1
#import mod_valid_utils as mvutil
importlib.reload(mcmp)


# Change: model run name/expt that is used to generate MOM restart
#         model restart date (input)
# Specify MOM grid: symmetric / non-symmetric (or "standard")
# u,v,h, q -points 
# in MOM6 - variables oriented in North-East order wrt p.ij
nrun       = "GOFS3.0"
expt       = "19.0"
YR         = 1993
jday       = 1
f_cont     = True        # True - continue unfinished file
grid_var   = 'vgrid'     # hgrid, ugrid, vgrid, qgrid
grid_shape = 'symmetr'   # MOM grid: symmetr/nonsymmetr

dnmb  = mtime.jday2dnmb(YR, jday)
dv    = mtime.datevec(dnmb)
ds    = mtime.datestr(dnmb)
MM    = dv[1]
DM    = dv[2]
rdate = f"{YR}{MM:02d}{DM:02d}"

pthout  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_restart/'
flgmap  = f"gofs2mom_nep_gmapi-{grid_var}_{grid_shape}.pkl"
dflgmap = os.path.join(pthout,flgmap)

print(f"FIXING gmapi NEP MOM6 --> GOFS: {grid_shape} {grid_var}\n")
print(f"Output will be saved --> {dflgmap}")

with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthrun  = dct[nrun][expt]["pthrun"]
pthgrid = dct[nrun][expt]["pthgrid"]
ftopo   = dct[nrun][expt]["ftopo"]
fgrid   = dct[nrun][expt]["fgrid"]

LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
HH           = np.where(HH>=0, np.nan, HH)

#
# Read MOM6 NEP domain nx=342, ny=816
pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom) 
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom) 
hlon, hlat  = mom6util.read_mom6grid(dfgrid_mom, grid=grid_shape, grdpnt=grid_var)
#qlon, qlat  = mom6util.read_mom6grid(dfgrid_mom, grid='symmetr', grdpnt='qgrid')
HHM         = mom6util.read_mom6depth(dftopo_mom) 
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]

Iall = np.where(HHM.flatten()<=0.)[0]
nall = Iall.shape[0]

# Convert to -180/180:
hlon   = np.where(hlon > 180.0, hlon-360., hlon)

import timeit
import mod_regmom as mrmom
import mod_bilinear as mblnr
importlib.reload(mblnr)

Irsb = []
print('Start from last saved record')
print('Loading ' + dflgmap)
try:
  with open(dflgmap,'rb') as fid:
    INDX, JNDX = pickle.load(fid)
  aa = np.squeeze(INDX[:,:,0])
  Irsb  = np.where(aa.flatten()>=0)[0]
  Ndone = len(Irsb)/nall
  print(f"In saved output: {Ndone*100:3.1f}% finished")
except:
  print(f"Cannot open output {dflgmap}")
  print("Starting from 0")


kcc = 0
tic = timeit.default_timer()
ticR = timeit.default_timer()
print('Searching HYCOM indices gmapi for MOM6 hpoints')

for ikk in range(nall):
  I1 = Iall[ikk]
  jj, ii = np.unravel_index(I1,HHM.shape)
  kcc += 1

  if (kcc % 2000) == 0 or I1 == Irsb[-1]:
    print(f" ---> Processed {kcc/nall*100.:3.1f}%...")

  II = np.squeeze(INDX[jj,ii,:])
  JJ = np.squeeze(JNDX[jj,ii,:])
  if np.min(II) < 0:
    print(f"End of processed indices ikk={ikk} II={II}")
    break
  II, JJ = mblnr.sort_gridcell_indx(II,JJ)

  x0  = hlon[jj,ii]
  y0  = hlat[jj,ii]

  ii1, jj1 = II[0], JJ[0]
  ii2, jj2 = II[1], JJ[1]
  ii3, jj3 = II[2], JJ[2]
  ii4, jj4 = II[3], JJ[3]

  xx  = np.array([LON[jj1,ii1], LON[jj2,ii2], LON[jj3,ii3], LON[jj4,ii4]])
  yy  = np.array([LAT[jj1,ii1], LAT[jj2,ii2], LAT[jj3,ii3], LAT[jj4,ii4]])

  XV, YV, x0c, y0c = mblnr.lonlat2xy_wrtX0(xx,yy,x0,y0)

# Map X,Y ---> Xhat, Yhat on reference quadrialteral 
# i.e. map GOFS grid coordinate to a reference quadrilateral 
# to do bilinear interpolation 
#  xht, yht = mblnr.map_x2xhat(xx, yy, x0, y0)
  xht, yht = mblnr.map_x2xhat(XV, YV, x0c, y0c)
  if abs(xht) > 1.001 or abs(yht) > 1.001:
    print(f"Fixing Error in bounding grid indices xht={xht} yht={yht}:")
    print(f" ==== x0={x0} y0={y0} ikk={ikk} jj={jj} ii={ii}")
    x0  = hlon[jj,ii]
    y0  = hlat[jj,ii]
    ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=0.12)
 
# Check:
    II = ixx.copy()
    JJ = jxx.copy()
    ii1, jj1 = II[0], JJ[0]
    ii2, jj2 = II[1], JJ[1]
    ii3, jj3 = II[2], JJ[2]
    ii4, jj4 = II[3], JJ[3]

    xx  = np.array([LON[jj1,ii1], LON[jj2,ii2], LON[jj3,ii3], LON[jj4,ii4]])
    yy  = np.array([LAT[jj1,ii1], LAT[jj2,ii2], LAT[jj3,ii3], LAT[jj4,ii4]])
    XV, YV, x0c, y0c = mblnr.lonlat2xy_wrtX0(xx,yy,x0,y0)

#    INp = mmisc1.inpolygon_1pnt(x0c, y0c, XV, YV)   
#    INp = mmisc1.inpolygon_1pnt(x0, y0, xx, yy)   
#    xht, yht = mblnr.map_x2xhat(xx, yy, x0, y0)
    xht, yht = mblnr.map_x2xhat(XV, YV, x0c, y0c)
#    xchk, ychk = mblnr.map_xhat2x(XV, YV, xht, yht)  # reverse mapping not working?
    if abs(xht) > 1.001 or abs(yht) > 1.001:
      print(f"Fixing Error failed:  xht={xht} yht={yht}:")
#      raise Exception(f"Fixing Error failed:  xht={xht} yht={yht}:")

    INDX[jj,ii,:] = ixx
    JNDX[jj,ii,:] = jxx

print('END FIXING: Saving to '+dflgmap)
with open(dflgmap,'wb') as fid:
  pickle.dump([INDX,JNDX],fid)

f_plt = False
if f_plt:
  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.3, 0.3, 0.3])
  ax1.plot(x0,y0,'r*')
#  ax1.plot(LON[jxx,ixx],LAT[jxx,ixx],'o-')
  ax1.plot(xx,yy,'o-')
  ax1.set_title(f'Bounding grid box pnt{x0:5.2f} {y0:5.2f}')

  ax2 = plt.axes([0.6, 0.3, 0.3, 0.3])
  ax2.plot(0,0,'r*')
  ax2.plot(XV,YV,'o-')
  ax2.set_title(f'Bounding grid box in cartes coord')




