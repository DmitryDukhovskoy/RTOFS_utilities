"""
  Find T/S profiles from GOFS3.1 using WOD observed T/S profiles
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

hg    = 1.e15
hg    = 1.e15
huge  = 1.e15
rg    = 9806.

# 
pthwod  = '/scratch1/NCEPDEV/stmp4/Dmitry.Dukhovskoy/WOD_profiles/'
pthrun  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/GOFS3.1/restart/'
#pthoutp = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/GOFS3.1/ts_prof/'
pthoutp= '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/MOM6_CICE6/ts_prof/'
pthbin  = pthrun


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

print(f'Extracting GOFS3.1 T/S for WOD UID {regn} {YR1} {YR2} mo={mo1}-{mo2}\n')

fuid_out = f'uidWOD_{regn}_{YR1}-{YR2}.pkl'
dflout   = os.path.join(pthoutp,fuid_out)
flts_out = f'gofsTS_{regn}_{YR1}-{YR2}.pkl'
dfltsout = os.path.join(pthoutp, flts_out)

print(f'Extracted T/S saved to ---> {dfltsout}')

# Load saved UID:
print(f'Loading {dflout}')
with open(dflout, 'rb') as fid:
  [SIDctd, SIDdrb, SIDpfl] = pickle.load(fid)

lctd = len(SIDctd)
ldrb = len(SIDdrb)
lpfl = len(SIDpfl)

print(f'# profiles CTD: {lctd}, DRB: {ldrb}, PFL: {lpfl}')

ZZi = mom6vld.zlevels()

pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
IDM = HH.shape[1]
JDM = HH.shape[0]

# For model use only 1 year - the only output available
YR    = 2021
MMin  = 0

# List of existing model files for processing, dates:
fall = os.listdir(pthrun)
LL   = [fnm for fnm in fall if not os.path.isdir(os.path.join(pthrun,fnm))]
LL   = [fnm for fnm in LL if (fnm[-2:] == '.a')]
nfls = len(LL)
TMG  = np.zeros((nfls))
JDG  = np.zeros((nfls))
for ifl in range(nfls):
  rdate    = LL[ifl][9:19]
  dnmb     = mtime.rdate2datenum(rdate)
  _, jday = mtime.dnmb2jday(dnmb)
  TMG[ifl] = dnmb
  JDG[ifl] = jday

TMG = np.sort(TMG)
JDG = np.sort(JDG)

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
      DD  = DM
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
# Find closest date
      if MM != MMin:
        dmm   = abs(JDG - jday)
        imm   = np.argmin(dmm)
        dnmbG = TMG[imm]
        DVG   = mmisc.datevec(dnmbG)
        YYG   = DVG[0]
        MMG   = DVG[1]  # months
        DMG   = DVG[2]  # month days
        
        flhycom = f"restart_r{YYG}{MMG:02d}{DMG:02d}00_930"
        fina    = pthbin + flhycom + '.a'
        finb    = pthbin + flhycom + '.b'
        MMin    = MM

        print(f"Reading thknss {fina}")
# Read layer thicknesses:
        dH, kdm = mhycom.read_hycom_restart(fina, finb, 'dp', IDM, JDM, finfo=True)
        dH = dH/rg
        dH = np.where(dH>huge, np.nan, dH)
        dH = np.where(dH<0.001, 0., dH)
        ZZ, ZM = mhycom.zz_zm_fromDP(dH, f_btm=False, finfo=False)

# Get T profile
        print(f"Reading temp {fina}")
        T3d, _ = mhycom.read_hycom_restart(fina, finb, 'temp', IDM, JDM, finfo=False)
        T3d = np.where(T3d >= huge, np.nan, T3d)
# Get S profile:
        print(f"Reading saln {fina}")
        S3d, _ = mhycom.read_hycom_restart(fina, finb, 'saln', IDM, JDM, finfo=False)
        S3d = np.where(S3d >= huge, np.nan, S3d)

      zm    = np.squeeze(ZM[:,jj0,ii0])
      Tprf  = np.squeeze(T3d[:,jj0,ii0])
      Ti    = mom6vld.interp_prof2zlev(Tprf, zm, ZZi, fill_surf=True)
      Sprf  = np.squeeze(S3d[:,jj0,ii0])
      Si    = mom6vld.interp_prof2zlev(Sprf, zm, ZZi, fill_surf=True)

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

  sttl = f'GOFS3.1 T \n{fina}'
  zlim = np.floor(zm[-1])
  ax1 = plt.axes([0.08, 0.1, 0.4, 0.8])
  ln1, = ax1.plot(Tprf,zm,'.-', label="obs")
  ln2, = ax1.plot(Ti,ZZi,'.-', label="intrp")
  ax1.set_ylim([zlim, 0])
  ax1.grid('on')
  ax1.set_title(sttl)

  sttl = f'GOFS3.1 S '
  ax2 = plt.axes([0.58, 0.1, 0.4, 0.8])
  ax2.plot(Sprf,zm,'.-', label="obs")
  ax2.plot(Si,ZZi,'.-', label="intrp")
  ax2.set_ylim([zlim, 0])
  ax2.grid('on')
  ax2.set_title(sttl)



