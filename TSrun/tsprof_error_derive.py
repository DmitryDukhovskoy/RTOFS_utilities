"""
  Compute error Argo and HYCOM T/S profile
  Argo profiles downloaded from HPSS QC log files
  see /scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/scripts/rtofs
  get_qctar.ssh

 use script to untar from HPSS HYCOM fields:
 scripts/rtofs/get_rtofs_archv.sh
 ./get_rtofs_archv.sh 20230128 f24
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

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text

import mod_utils as mutil
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_read_ncoda as rncoda
#importlib.reload(rncoda)
import mod_extrargo as exargo
importlib.reload(exargo)
from mod_utils import tserr


# What RTOFS output field to plot:
# bkgrd  - background, f/cast from the previous day  n00 day before
# incup  - incr. updated NCODA increments from incup files (6hr run) -n24 
# fcst0  - state after 24hr hindcast ready for forecast n00 fcast day
# fcst12 - 12hr forecast 

fcont      = True   # if true - continue from last saved, otherwise rewrite
rtofs_outp = 'incup'
#rtofs_outp = 'bkgrd'
rdate0     = '20230308'

# Specify Argo to plot
Srjct   = 4.  # QC STD rejection in NCODA for S profile
Trjct   = 4.  # QC STD rejection for T profile

yr, mo, mday, hr = rncoda.parse_rdate(rdate0)
time_n00    = datetime.datetime(yr,mo,mday,hr,0,0)  # forecast time
# incup:
time_n24    = time_n00 - datetime.timedelta(days=1) # time of hindcast incup incorpor
hdate_n24   = time_n24.strftime('%Y%m%d%H')
rdate_n24   = rdate0
# background:
time_bkgrd  = time_n00 - datetime.timedelta(days=1)
hdate_bkgrd = time_bkgrd.strftime('%Y%m%d%H') 
rdate_bkgrd = time_bkgrd.strftime('%Y%m%d')
# 12hr f/cast:
time_f12    = datetime.datetime(yr,mo,mday,hr,12,0)
hdate_f12   = time_f12.strftime('%Y%m%d%H')
rdate_f12   = rdate0

pthbin  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+\
          rdate0+'/ocnqc_logs/profile_qc/'
pthout  = pthbin
pthhcm0 = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'
#pthhcmB = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+\
#          rdate_bkgrd+'/'

# Read Argo profiles
flnm = 'prof_argo_rpt.'+hdate_n24+'.txt'
TS = exargo.extract_argo(pthbin,flnm)
# dir(TS)  to see attributes
if len(TS.Tprof) == 0:
  raise Exception('Argo float was not found')

import mod_read_hycom as mhycom
import mod_extrTS_hycom as exhcm
importlib.reload(exhcm)

IDM = 4500
JDM = 3298
KDM = 41

pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo   = 'regional.depth'
fgrid   = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid) 


if rtofs_outp == "incup":
  time_rtofs = time_n24     # actual time of the output field
  flhcm      = 'rtofs_glo.t00z.n-24.archv'
  pthhcm     = pthhcm0 + rdate_n24 + "/"
elif rtofs_outp == "bkgrd":
  time_rtofs = time_bkgrd     # actual time of the output field
  flhcm      = 'rtofs_glo.t00z.n00.archv'
  pthhcm     = pthhcm0 + rdate_bkgrd + '/'
 
SL2   = []
TL2   = []
SRMSE = []
TRMSE = []
SINF  = []
TINF  = []



TSERR  = tserr() 
nprof  = len(TS.Sprof)
nprc  = 40
kk1   = 440
kk2    = kk1+(nprc-1)
# kk2   = nprof
if kk2 >= nprof:
  kk2 = nprof-1

#fout   = pthout + 'TSprof_stat_{2}_{0:04d}-{1:04d}.pkl'.format(kk1,kk2,rtofs_outp)
fout      = pthout + 'TSprof_stat_{0:04d}-{1:04d}.pkl'.format(kk1,kk2)
frqout    = 10 
ncc       = 0
nlast     = 0
rcrd_last = kk1-1
print('Total profiles = {0}, saving {1}-{2}'.format(nprof,kk1,kk2))
print('Output file: '+fout)

if fcont:
  print('Loading saved output '+fout)
  try:
    with open(fout,'rb') as fid:
      TSERR   = pickle.load(fid)
    nlast     = len(TSERR.numb)
    rcrd_last = TSERR.recn[-1] 
    print('Found {0} saved records, last #={1}'.format(nlast,rcrd_last))

  except:
    print('Saved output not found, starting from {0} rec={1}'.\
           format(ncc,kk1))


for kk in range(kk1,kk2+1):
  tkk   = datetime.datetime.now()
  Anmb  = TS.numb[kk]
  latar = TS.lat[kk]
  lonar = TS.lon[kk]
  
# Check time: ARGO profiles should be ~1 day back from the forecast rdate
# Saved is observation time, received time is Rcpt in txt ~2 hrs later
  dnmb_argo = TS.ptime[kk]
  dstr_argo = mmisc.datestr(dnmb_argo)
  rnmb_argo = TS.rtime[kk]  # time of data recieved
  rstr_argo = mmisc.datestr(rnmb_argo)
  print('{0} Argo# {1}, Lon={4} Lat={5} TimeObs: {2} TimeReceived: {3}'.\
        format(kk+1,Anmb,mmisc.datestr(dnmb_argo), mmisc.datestr(rnmb_argo), \
               lonar,latar))

  ncc += 1   # record counter
  if kk <= rcrd_last:
    print('Record {0} #{1} exist, skipping ...'.format(ncc,kk+1))
    continue

#print(TS.__dir__()R.numb) see attributes of TS
  Th, ZZh, ZMh, dHh = exhcm.extract_1prof(pthhcm,flhcm,'temp',lonar,latar,\
                                        LON,LAT,KDM=KDM)
  Sh, dm1, dm2, dm3 = exhcm.extract_1prof(pthhcm,flhcm,'salin',lonar,latar,\
                                        LON,LAT,fthknss=False,KDM=KDM)
#
# Argo profiles:
  ZSar = np.array(TS.ZS[kk]).squeeze()
  Sar  = np.array(TS.Sprof[kk]).squeeze()
  ZTar = np.array(TS.ZT[kk]).squeeze()
  Tar  = np.array(TS.Tprof[kk]).squeeze()
  ZSar_min = np.min(ZSar)
  ZTar_min = np.min(ZTar)


# Interpolate ARGO profiles onto HYCOM layers:
  iSar = mutil.interp_1Dhycom(ZSar,Sar,ZZh,ZMh)
  iTar = mutil.interp_1Dhycom(ZSar,Tar,ZZh,ZMh)

# HYCOM indices for Argo obs:
# nearest point:
  Ih, Jh = mutil.find_indx_lonlat(lonar,latar,LON,LAT)
# exact index of Argo float:
  Iex, Jex = mutil.interp_indx_lonlat(lonar,latar,LON,LAT) 

# Check, Iex and Ih should be close:
  if abs(Iex-Ih)>1:
    raise Exception('Check indices: Ih={0} Iex={1}'.format(Ih,Iex))
  elif abs(Jex-Jh)>1:
    raise Exception('Check indices: Jh={0} Jex={1}'.format(Jh,Jex))
# Error:
  Sl2, Srmse, Sinf = mutil.err_stat(Sh, iSar)
  Tl2, Trmse, Tinf = mutil.err_stat(Th, iTar)

# Pool data together:
# for scatter plot
  Zbins = np.array([0,-200,-6000])
  Tbin = mutil.pool_data_zbins(Th,iTar,ZMh,Zbins)
  Sbin = mutil.pool_data_zbins(Sh,iSar,ZMh,Zbins)
  SS01 = Sbin.T1b       # depth-binned S data profile 1 (hycom)
  SS02 = Sbin.T2b       # depth-binned S data profile 2 (argo)
  TT01 = Tbin.T1b       # depth-binned T data prof 1
  TT02 = Tbin.T2b       # depth-binned T data prof 2
 
# Rejected?
  Tstd = TS.Tstd[kk]
  Sstd = TS.Sstd[kk]
  Sqc_rjct = -999
  Tqc_rjct = -999
  if Sstd > Srjct:
    Sqc_rjct = kk
  if Tstd > Trjct:
    Tqc_rjct = kk

  TSERR.add_data(kk, Anmb, lonar, latar, dnmb_argo, rnmb_argo, Tstd, Sstd, \
                 Sl2, Tl2, Srmse, Trmse, Sinf, Tinf, SS01, TT01,\
                 SS02, TT02, Sqc_rjct, Tqc_rjct)    

  if kk%frqout == 0 and kk > kk1:
    print('Saving TS stat '+fout)
    with open(fout,'wb') as fid:
      pickle.dump(TSERR,fid)

  tend = datetime.datetime.now()
  print('Processing 1 rec: {0}'.format(tend-tkk))


print(' End of cycle ')
print('Saving TS stat '+fout)
with open(fout,'wb') as fid:
  pickle.dump(TSERR,fid)
  
print(' All done ' )




