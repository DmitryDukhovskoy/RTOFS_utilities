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
sys.path.append(PPTHN + '/MyPython/ncoda_utils')

from mod_utils_fig import bottom_text
import mod_mom6_valid as mom6vld
import mod_read_hycom as mhycom
import mod_misc1 as mmisc
import mod_time as mtime
import mod_WODdata as mwod
import mod_mom6 as mom6util
import mod_utils as mutil
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
#regn   = 'AmundsAO'
#regn   = 'NansenAO'
#regn   = 'MakarovAO'
regn  = 'CanadaAO'

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

# For model use only 1 year - the only output available
YRmom = 2021
YR    = YRmom
MMin  = 0
# Keep truck of model indices to avoid duplicates:
Indx  = np.array(([]))
Jndx  = np.array(([]))
LI    = []
LJ    = []
for ii in range(12):
  LI.append(Indx)
  LJ.append(Jndx)

f_new = True
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
      lat0 = mwod.read_ncfld(dflnm, 'lat', finfo=False)
      lon0 = mwod.read_ncfld(dflnm, 'lon', finfo=False)
      TMM = mwod.read_ncfld(dflnm,'time')
    # time:units = "days since 1770-01-01 00:00:00" ;
      dnmb_ref = mmisc.datenum([1770,1,1])
      TM  = TMM + dnmb_ref
      DV  = mmisc.datevec(TM)
      YY  = DV[0]
      MM  = DV[1]  # months
      DM  = DV[2]  # month days
      jday= int(mtime.date2jday([YY,MM,DM]))
      HR  = 12

# Find indices in MOM6 grid:
      ii0, jj0 = mutil.find_indx_lonlat(lon0, lat0, LON, LAT)
# Check duplicates for this month:
      Indx = LI[MM-1]
      Jndx = LJ[MM-1]
      if len(Indx) > 0:
        dmm = np.sqrt((Indx - ii0)**2 + (Jndx - jj0)**2)
        if np.min(dmm) < 1.e-10:
          print(f'Skipping Duplicate i/j: {ii0}/{jj0}')
          continue

      Indx = np.append(Indx,[ii0])
      Jndx = np.append(Jndx,[jj0])
      LI[MM-1] = Indx
      LJ[MM-1] = Jndx

# If loaded data from the same month - use what is loaded
# Day is not important for this analysis
      if MM != MMin:
        pthbin = pthrun + 'tarmom_{0}{1:02d}/'.format(YR,MM)
        flmom  = 'ocnm_{0}_{1:03d}_{2}.nc'.format(YR,jday,HR)
        flin   = pthbin + flmom
        MMin   = MM

        print(f"Reading thknss {flin}")
# Read layer thicknesses:
        dH     = mom6util.read_mom6(flin, 'h', finfo=False)
        ssh    = mom6util.read_mom6(flin, 'SSH', finfo=False)
        ZZ, ZM = mom6util.zz_zm_fromDP(dH, ssh, f_intrp=True, finfo=False)

# Get T profile
        print(f"Reading pot T")
        T3d   = mom6util.read_mom6(flin, 'potT', finfo=False)
# Get S profile:
        S3d   = mom6util.read_mom6(flin, 'salt', finfo=False)

      zm    = np.squeeze(ZM[:,jj0,ii0])
      Tprf  = np.squeeze(T3d[:,jj0,ii0])
      Ti    = mom6vld.interp_prof2zlev(Tprf, zm, ZZi)
      Sprf  = np.squeeze(S3d[:,jj0,ii0])
      Si    = mom6vld.interp_prof2zlev(Sprf, zm, ZZi)

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

btx = 'extract_TSprofWOD_mom6.py'
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

  sttl = f'MOM6-003 T \n{flin}'
  zlim = np.floor(zm[-1])
  ax1 = plt.axes([0.08, 0.1, 0.4, 0.8])
  ln1, = ax1.plot(Tprf,zm,'.-', label="obs")
  ln2, = ax1.plot(Ti,ZZi,'.-', label="intrp")
  ax1.set_ylim([zlim, 0])
  ax1.grid('on')
  ax1.set_title(sttl)

  sttl = f'MOM6-003 S '
  ax2 = plt.axes([0.58, 0.1, 0.4, 0.8])
  ax2.plot(Sprf,zm,'.-', label="obs")
  ax2.plot(Si,ZZi,'.-', label="intrp")
  ax2.set_ylim([zlim, 0])
  ax2.grid('on')
  ax2.set_title(sttl)



