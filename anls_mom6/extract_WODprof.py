"""
  T/S Profile Observations - downloaded from WOD18 website
  https://www.ncei.noaa.gov/access/world-ocean-database/bin/getwodyearlydata.pl

  First, run derive_WODuid.py to select profiles for specified regions/time

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
importlib.reload(mom6vld)


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

print(f'Extracting WOD T/S for selected UID {regn} {YR1} {YR2} mo={mo1}-{mo2}\n')


fuid_out = f'uidWOD_{regn}_{YR1}-{YR2}.pkl'
dflout   = os.path.join(pthoutp,fuid_out)
flts_out = f'WODTS_{regn}_{YR1}-{YR2}.pkl'
dfltsout = os.path.join(pthoutp, flts_out)

# Load saved UID:
print(f'Loading {dflout}')
with open(dflout, 'rb') as fid:
  [SIDctd, SIDdrb, SIDpfl] = pickle.load(fid)

lctd = len(SIDctd)
ldrb = len(SIDdrb)
lpfl = len(SIDpfl)

print(f'# profiles CTD: {lctd}, DRB: {ldrb}, PFL: {lpfl}')

ZZi = mom6vld.zlevels()

import mod_swstate as msw
lat_ref = 85.
pr_ref  = 0.
f_new   = True

# CTD
for obtype  in ['CTD', 'DRB', 'PFL']:
  if obtype == 'CTD':
    ll  = lctd
    SID = SIDctd
  elif obtype == 'DRB':
    ll  = ldrb
    SID = SIDdrb
  elif obtype == 'PFL':
    ll  = lpfl
    SID = SIDpfl

  if ll > 0:
    print(f'Interpolating T/S prof onto Z {obtype}')
    for ii in range(ll):
      if ii%100 == 0:
        print(f' {float(ii/ll)*100.:.1f}% done ...')

      uid     = SID[ii]
      flnm    = f'wod_0{uid}O.nc'
      pthdata = os.path.join(pthwod, obtype)
      dflnm   = os.path.join(pthdata, flnm)

      ZZ, T, S, pr_db = mwod.readTS_WODuid(dflnm, f_info=False)  # in situ T !

  # COnvert in situ -> potential T:
      if len(pr_db) == 0:
        lat = mwod.read_ncfld(dflnm, 'lat', finfo=False)
        lat_obs = np.zeros((len(ZZ)))*0.0 + lat
        pr_db, pr_pa = msw.sw_press(ZZ, lat_obs)

      Tpot = msw.sw_ptmp(S, T, pr_db, pr_ref)
  # Interpolate
      Ti = mom6vld.interp_prof2zlev(Tpot, ZZ, ZZi)
      Si = mom6vld.interp_prof2zlev(S, ZZ, ZZi)

      if f_new:
        SPROF = mom6vld.PROF1D(ZZi,Si)
        TPROF = mom6vld.PROF1D(ZZi,Ti)
        f_new = False
      else:
        SPROF.add_array(Si)
        TPROF.add_array(Ti)

if f_save:
  print(f'Saving --> {dfltsout}')
  with open(dfltsout, 'wb') as fid:
    pickle.dump([SPROF, TPROF],fid)  

print('All Done\n')



btx = 'extract_WODprof.py'
f_pltobs = False
if f_pltobs:
  from mpl_toolkits.basemap import Basemap, cm
  import matplotlib.colors as colors
  import matplotlib.mlab as mlab
  plt.ion()
  clr_drb = [0.8, 0.4, 0]
  clr_pfl = [0, 0.4, 0.8]
  clr_ctd = [1, 0, 0.8]


  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()

  sttl = f'UID={uid} {obtype} Tpot'
  zlim = np.floor(ZZ[-1])
  ax1 = plt.axes([0.08, 0.1, 0.4, 0.8])
  ln1, = ax1.plot(Tpot,ZZ,'.-', label="obs")
  ln2, = ax1.plot(Ti,ZZi,'.-', label="intrp")
  ax1.set_ylim([zlim, 0])
  ax1.grid('on')
  ax1.set_title(sttl)

  sttl = f'UID={uid} {obtype} S'
  zlim = np.floor(ZZ[-1])
  ax2 = plt.axes([0.58, 0.1, 0.4, 0.8])
  ax2.plot(S,ZZ,'.-', label="obs")
  ax2.plot(Si,ZZi,'.-', label="intrp")
  ax2.set_ylim([zlim, 0])
  ax2.grid('on')
  ax2.set_title(sttl)






