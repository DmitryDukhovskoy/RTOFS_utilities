"""
  T/S Profile Observations - downloaded from WOD18 website
  https://www.ncei.noaa.gov/access/world-ocean-database/bin/getwodyearlydata.pl
  Derive WOD UID for selected regions and time

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
from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text
import mod_mom6_valid as mom6vld
import mod_read_hycom as mhycom
import mod_misc1 as mmisc
import mod_time as mtime
import mod_WODdata as mwod
importlib.reload(mwod)


# WOD data
pthwod = '/scratch1/NCEPDEV/stmp4/Dmitry.Dukhovskoy/WOD_profiles/'

pthoutp= '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/MOM6_CICE6/ts_prof/'


# Select lon, lat to search for WOD profiles
f_save = True
YR1    = 2018
YR2    = 2022
mo1    = 1
mo2    = 12
#regn   = 'AmundsAO'
#regn   = 'NansenAO'
#regn   = 'MakarovAO'
regn  = 'CanadaAO'

REGNS = mom6vld.ts_prof_regions()
x0 = REGNS[regn]["lon0"]
y0 = REGNS[regn]["lat0"] 
dx = REGNS[regn]["dlon"]
dy = REGNS[regn]["dlat"]

print(f'Finding WOD UID for {regn} {YR1} {YR2} mo={mo1}-{mo2}\n')

obtype = 'CTD'
furl  = '{0}{1}/info_{1}.nc'.format(pthwod,obtype)
Yctd, Xctd, TMctd, UIDctd = mwod.search_UID(furl, x0, y0, dx, dy, \
                            YR1=YR1, YR2=YR2, mnth1=mo1, mnth2=mo2)
DVctd = mtime.datevec1D(TMctd)

obtype = 'DRB'
furl  = '{0}{1}/info_{1}.nc'.format(pthwod,obtype)
Ydrb, Xdrb, TMdrb, UIDdrb = mwod.search_UID(furl, x0, y0, dx, dy, \
                            YR1=YR1, YR2=YR2, mnth1=mo1, mnth2=mo2)
DVdrb = mtime.datevec1D(TMdrb)

obtype = 'PFL'
furl  = '{0}{1}/info_{1}.nc'.format(pthwod,obtype)
Ypfl, Xpfl, TMpfl, UIDpfl = mwod.search_UID(furl, x0, y0, dx, dy, \
                            YR1=YR1, YR2=YR2, mnth1=mo1, mnth2=mo2)
DVpfl = mtime.datevec1D(TMpfl)

fuid_out = f'uidWOD_{regn}_{YR1}-{YR2}.pkl'
dflout   = os.path.join(pthoutp,fuid_out)

#
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

# Subset by depths:
Hmin = -1000. 

obtype = 'CTD'
UID = UIDctd  
pthdata = os.path.join(pthwod, obtype)
Jctd    = []
SIDctd  = []
if len(UID) > 0:
  print('Depth: Subsetting CTD')
  Jctd, SIDctd = mwod.select_WODdepth(pthdata, UID, Hmin, LON, LAT, HH, qflag = True)

obtype = 'DRB'
UID = UIDdrb  
pthdata = os.path.join(pthwod, obtype)
Jdrb    = []
SIDdrb  = []
if len(UID) > 0:
  print('Depth: Subsetting DRB')
  Jdrb, SIDdrb = mwod.select_WODdepth(pthdata, UID, Hmin, LON, LAT, HH, qflag = True)

obtype = 'PFL'
UID = UIDpfl  
pthdata = os.path.join(pthwod, obtype)
Jpfl    = []
SIDpfl  = []
if len(UID) > 0:
  print('Depth: Subsetting PFL')
  Jpfl, SIDpfl = mwod.select_WODdepth(pthdata, UID, Hmin, LON, LAT, HH, qflag = True)

# Subset coordinates of obs:
Xctd0 = Xctd.copy()
Yctd0 = Yctd.copy()
Xdrb0 = Xdrb.copy()
Ydrb0 = Ydrb.copy()
Xpfl0 = Xpfl.copy()
Ypfl0 = Ypfl.copy()
Xctd = Xctd[Jctd]
Yctd = Yctd[Jctd]
Xdrb = Xdrb[Jdrb]
Ydrb = Ydrb[Jdrb]
Xpfl = Xpfl[Jpfl]
Ypfl = Ypfl[Jpfl]

print('\n ======================')
print(f'Selected: CTD={len(SIDctd)} DRB={len(SIDdrb)} PFL={len(SIDpfl)}\n')

if f_save:
  print(f'Saving --> {dflout}')
  with open(dflout, 'wb') as fid:
    pickle.dump([SIDctd, SIDdrb, SIDpfl],fid)
#pickle.dump([II, JJ, XX, YY, Lsgm, Hbtm], fid)

btx = 'derive_WODuid.py'
f_pltobs = True
if f_pltobs:
  from mpl_toolkits.basemap import Basemap, cm
  import matplotlib.colors as colors
  import matplotlib.mlab as mlab
  plt.ion()
  clr_drb = [0.8, 0.4, 0]
  clr_pfl = [0, 0.4, 0.8]
  clr_ctd = [1, 0, 0.8]

  lon0=-10
  lat0=70
  obtype = 'DRB'
  sttl = f'{regn} WOD obs {YR1}-{YR2} {mo1}-{mo2}'
  ax1 = mom6vld.plot_points_orthomap(Xdrb,Ydrb,LON,LAT,HH, clr=clr_drb,\
                                     lon0=lon0, lat0=lat0, btx=btx, sttl=sttl)

  m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution='l')
  xpp, ypp = m(Xpfl,Ypfl)  
  m.plot(xpp,ypp,'.', color=clr_pfl)

# Plot specified region:
  f_regn = True
  if f_regn:
    ddx = 5
    nstps = int(2*dx/ddx) + 1
    XR  = [x0-dx]
    YR  = [y0-dy]
    for icc in range(nstps):
      xnxt = (x0-dx) + icc*ddx
      if xnxt > x0+dx:
        break
      XR.append(xnxt)
      YR.append(y0+dy)

    for icc in range(nstps):
      xnxt = (x0+dx) - icc*ddx
      if xnxt < x0-dx:
        break
      XR.append(xnxt)
      YR.append(y0-dy)

    XRp, YRp = m(XR,YR)
    m.plot(XRp,YRp,'-', color=(1,0,0))

  if len(Xctd) > 0:
    xcc, ycc = m(Xctd, Yctd)
    m.plot(xcc, ycc, '.', color=clr_ctd)

  ax2 = plt.axes([0.78,0.82,0.18,0.15])
  ax2.plot(0, 0.9, '.', ms=14, color=clr_ctd)
  ax2.text(0.07, 0.9, 'CTD')
  ax2.plot(0, 0.5, '.', ms=14, color=clr_drb)
  ax2.text(0.07, 0.5, 'DRB')
  ax2.plot(0, 0.1, '.', ms=14, color=clr_pfl)
  ax2.text(0.07, 0.1, 'PFL')
  ax2.set_xlim([-0.1, 0.4])
  ax2.set_ylim([0, 2])
  ax2.axis('off')




