"""
  Find T/S profiles from MOM6 using WOD observed T/S profiles
  Subsampling for the same locations / months

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

expt  = '003'
hg    = 1.e15

# WOD data
pthwod = '/scratch1/NCEPDEV/stmp4/Dmitry.Dukhovskoy/WOD_profiles/'
pthoutp= '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/MOM6_CICE6/ts_prof/'
pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_' + expt + '/'


# Select lon, lat to search for WOD profiles
f_save = True
YR1    = 2018
YR2    = 2022
mo1    = 1
mo2    = 12
regn   = 'AmundsAO'
#regn   = 'NansenAO'
#regn   = 'MakarovAO'
#regn  = 'CanadaAO'

REGNS = mom6vld.ts_prof_regions()
x0 = REGNS[regn]["lon0"]
y0 = REGNS[regn]["lat0"]
dx = REGNS[regn]["dlon"]
dy = REGNS[regn]["dlat"]

print(f'Extracting MOM6 T/S for WOD UID {regn} {YR1} {YR2} mo={mo1}-{mo2}\n')


fuid_out = f'uidWOD_{regn}_{YR1}-{YR2}.pkl'
dflout   = os.path.join(pthoutp,fuid_out)
flts_out = f'mom6TS_{regn}_{YR1}-{YR2}.pkl'
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

pthgrid   = pthrun + 'INPUT/'
fgrd_mom  = pthgrid + 'regional.mom6.nc'
ftopo_mom = pthgrid + 'ocean_topog.nc'
LON, LAT  = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
HH        = mom6util.read_mom6depth(ftopo_mom)

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
      lat = mwod.read_ncfld(dflnm, 'lat', finfo=False)
      lon = mwod.read_ncfld(dflnm, 'lon', finfo=False)
      TMM = mwod.read_ncfld(dflnm,'time')
    # time:units = "days since 1770-01-01 00:00:00" ;
      dnmb_ref = mmisc.datenum([1770,1,1])
      TM  = TMM + dnmb_ref
      DV  = mmisc.datevec(TM)
      YY  = DV[0]
      MM  = DV[1]  # months
      DM  = DV[2]  # month days
      jday= int(mtime.date2jday([YR,MM,DD]))
      HR  = 12

      pthbin = pthrun + 'tarmom_{0}{1:02d}/'.format(YR,MM)
      flmom  = 'ocnm_{0}_{1:03d}_{2}.nc'.format(YR,jday,HR)
      flin   = pthbin + flmom

      A3d    = mom6util.read_mom6(flin, 'potT', finfo=False)
# Find indices in MOM6 grid:

# Get T profile

#
      A3d    = mom6util.read_mom6(flin, 'salt', finfo=False)
