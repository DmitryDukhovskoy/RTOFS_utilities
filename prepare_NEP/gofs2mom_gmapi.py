# To prepare MOM restart:
# 1st step: interpoalte T,S,U, fields horizontally 
# from HYCOM archv file onto MOM grid
# 
# Keep fields on HYCOM vertical grid
# 
# Plot NEP domain on GOFS grid
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
grid_var   = 'ugrid'     # hgrid, ugrid, vgrid, qgrid
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

print(f"Preparing gmapi NEP MOM6 --> GOFS: {grid_shape} {grid_var}\n")
print(f"Output will be saved --> {dflgmap}")

with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthrun  = dct[nrun][expt]["pthrun"]
pthgrid = dct[nrun][expt]["pthgrid"]
ftopo   = dct[nrun][expt]["ftopo"]
fgrid   = dct[nrun][expt]["fgrid"]

LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
HH           = np.where(HH>=0, np.nan, HH)
DX, DY       = mhycom.dx_dy(LON, LAT)
RR           = np.sqrt(DX**2 + DY**2)

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


# To speed up, subset global domain:
# Longitudes - does not work on tri-polar grid
# Do manual subsetting for now
#lat_min = np.min(hlat)
#lat_max = np.max(hlat)
#if np.min(hlon)<0:
#  aa = hlon.copy()
#  aa = np.where(aa<0, aa+360., aa)
#  lon_min = np.min(aa) 
#  if lon_min > 180.:
#    lon_min = lon_min - 360.
#else:
#  lon_min = np.min(hlon)
#
#if np.min(hlon)<0:
#  lon_max = np.max(aa)
#  if lon_max > 180.:
#    lon_max = lon_max - 360.
#else:
#  lon_max = np.max(hlon)
#JJ,II = np.where(LAT < lat_min)
#Jsub1 = np.max(JJ)
#JJ,II = np.where(LAT > lat_max)
#Jsub2 = np.min(JJ)

#  aa = LON[Jsub1:Jsub2,:]
#  bb = LAT[Jsub1:Jsub2,:]
#  # Only for domains not across the boundaries:
#  JJ,II = np.where(aa < lon_min)
#  Isub1 = np.max(II)
#  JJ,II = np.where(aa > lon_max)
#  Isub2 = np.min(II)

Jsub1 = 1600
Jsub2 = 3030
Isub1 = 1010 
Isub2 = 2290

#LONsub = LON[Jsub1:Jsub2+1, Isub1:Isub2+1]
#LATsub = LAT[Jsub1:Jsub2+1, Isub1:Isub2+1]
#RRsub  = RR[Jsub1:Jsub2+1, Isub1:Isub2+1]

INDX = np.zeros((jdm,idm,4)).astype(int)-999
JNDX = np.zeros((jdm,idm,4)).astype(int)-999

import timeit
import mod_regmom as mrmom

Irsb = []
if f_cont:
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

# If continue from saved:
  if len(Irsb) > 0:
    if np.min(abs(Irsb-I1)) == 0:
      if (kcc % 1000) == 0 or I1 == Irsb[-1]:
        print(f" ---> Skipping {kcc/nall*100.:3.1f}%...")
      continue

  if (kcc % 2000) == 0:
    toc    = timeit.default_timer()
    pprc   = kcc/nall*100.
    dtic   = (toc-ticR)/60.
    tictot = (toc-tic)/60.
    print(f"  {pprc:4.1f}% done dt={dtic:6.2f} min, ttot={tictot:7.1f} min")
    ticR = timeit.default_timer()

#  t1 = timeit.default_timer() 
  x0  = hlon[jj,ii]
  y0  = hlat[jj,ii]
  ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=0.12)
#  t2 = timeit.default_timer()
#  print(f"Process time: {t2-t1}")

  INDX[jj,ii,:] = ixx
  JNDX[jj,ii,:] = jxx

  if kcc%10000 == 0:
    print('Saving to ' + dflgmap)
    with open(dflgmap,'wb') as fid:
      pickle.dump([INDX,JNDX],fid)


print('END: Saving to '+dflgmap)
with open(dflgmap,'wb') as fid:
  pickle.dump([INDX,JNDX],fid)

f_plt = False
if f_plt:
  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])

  ax1.plot(x0,y0,'r*')
  ax1.plot(LON[jxx,ixx],LAT[jxx,ixx],'o')






