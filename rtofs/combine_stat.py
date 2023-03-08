"""
  Combine pkl error stat files into a single file
  after running serial tsprof_error_anls.py 
  for a set of Argo profiles
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

sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

import mod_utils as mutil
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_read_ncoda as rncoda
#importlib.reload(rncoda)
import mod_extrargo as exargo
importlib.reload(exargo)

# What RTOFS output field to plot:
# bkgrd  - background, f/cast from the previous day  n00 day before
# incup  - incr. updated NCODA increments from incup files (6hr run) -n24 
# fcst0  - state after 24hr hindcast ready for forecast n00 fcast day
# fcst12 - 12hr forecast 

rtofs_outp = 'incup'
#rtofs_outp = 'bkgrd'
rdate0     = '20230129'  # forecast date where Argo is assimilated

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


class tserr():
  kind = 'TS Argo error stat'

  def __init__(self):
    self.recn  = []
    self.numb  = []
    self.lon   = []
    self.lat   = []
    self.ptime = []
    self.rtime = []
    self.Tstd  = []
    self.Sstd  = []
    self.SL2   = []
    self.TL2   = []
    self.SRMSE = []
    self.TRMSE = []
    self.SINF  = []
    self.TINF  = []
    self.SS01  = []  # pooled S values in the top level 01
    self.TT01  = []
    self.SS02  = []
    self.TT02  = []
    self.Sqcrj = []  # indices of rejected S profiles
    self.Tqcrj = []  # omdoces pf rejected T profiles

  def add_data(self, recn, numb, lon, lat, dnmb, rnmb, Tstd, Sstd, \
               SL2, TL2, SRMSE, TRMSE, SINF, TINF, SS01, TT01, \
               SS02, TT02, Sqcrj, Tqcrj):
    self.recn.append(recn)    # record #
    self.numb.append(numb)    # Argo #
    self.lon.append(lon)
    self.lat.append(lat)
    self.ptime.append(dnmb)
    self.rtime.append(rnmb)
    self.Tstd.append(Tstd)
    self.Sstd.append(Sstd)
    self.SL2.append(SL2)
    self.TL2.append(TL2)
    self.SRMSE.append(SRMSE)
    self.TRMSE.append(TRMSE)
    self.SINF.append(SINF)
    self.TINF.append(TINF)
    self.SS01.append(SS01)     #values in the top level 01
    self.TT01.append(TT01)
    self.SS02.append(SS02)
    self.TT02.append(TT02)
    self.Sqcrj.append(Sqcrj) #S rejected profiles
    self.Tqcrj.append(Tqcrj)

class tsarray():
  def __init__(self,recn, numb, lon, lat, dnmb, rnmb, Tstd, Sstd, \
               SL2, TL2, SRMSE, TRMSE, SINF, TINF, SS01, TT01, \
               SS02, TT02, Sqcrj, Tqcrj):
    self.recn  = recn
    self.numb  = numb
    self.lon   = lon
    self.lat   = lat
    self.ptime = dnmb
    self.rtime = rnmb
    self.Tstd  = Tstd
    self.Sstd  = Sstd
    self.SL2   = SL2
    self.TL2   = TL2
    self.SRMSE = SRMSE
    self.TRMSE = TRMSE
    self.SINF  = SINF
    self.TINF  = TINF
    self.SS01  = SS01
    self.TT01  = TT01
    self.SS02  = SS02
    self.TT02  = TT02
    self.Sqcrj = Sqcrj
    self.Tqcrj = Tqcrj

  def add_array(self, recn, numb, lon, lat, dnmb, rnmb, Tstd, Sstd, \
               SL2, TL2, SRMSE, TRMSE, SINF, TINF, SS01, TT01, \
               SS02, TT02, Sqcrj, Tqcrj):
    self.recn  = np.append(self.recn,recn, axis=0)
    self.numb  = np.append(self.numb,numb, axis=0)
    self.lon   = np.append(self.lon,lon, axis=0) 
    self.lat   = np.append(self.lat,lat, axis=0) 
    self.ptime = np.append(self.ptime,dnmb, axis=0)
    self.rtime = np.append(self.rtime,rnmb, axis=0)
    self.Tstd  = np.append(self.Tstd,Tstd, axis=0)
    self.SL2   = np.append(self.SL2,SL2, axis=0) 
    self.TL2   = np.append(self.TL2,TL2, axis=0) 
    self.SRMSE = np.append(self.SRMSE,SRMSE, axis=0)
    self.TRMSE = np.append(self.TRMSE,TRMSE, axis=0)
    self.SINF  = np.append(self.SINF,SINF, axis=0)
    self.TINF  = np.append(self.TINF,TINF, axis=0)
    self.SS01  = np.append(self.SS01,SS01, axis=0)
    self.TT01  = np.append(self.TT01,TT01, axis=0)
    self.SS02  = np.append(self.SS02,SS02, axis=0)
    self.TT02  = np.append(self.TT02,TT02, axis=0)
    self.Sqcrj = np.append(self.Sqcrj,Sqcrj, axis=0)
    self.Tqcrj = np.append(self.Tqcrj,Tqcrj, axis=0)

# Save as np array
#TSERR   = tsarray()
crc     = 0
pthbin  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+\
          rdate0+'/ocnqc_logs/profile_qc/'
pthout  = pthbin
flbase  = 'TSprof_stat_'
lnnm    = len(flbase)
for fll in os.listdir(pthbin):
  ffpth = os.path.join(pthbin,fll)
  if os.path.isfile(ffpth) and fll[:lnnm]==flbase:
    print(fll)

    with open(ffpth, 'rb') as fid:
      DMM = pickle.load(fid)

    recn = np.array(DMM.recn)
    numb = np.array(DMM.numb)
    lon  = np.array(DMM.lon)
    lat  = np.array(DMM.lat)
    dnmb = np.array(DMM.ptime)
    rnmb = np.array(DMM.rtime)
    Tstd = np.array(DMM.Tstd)
    Sstd = np.array(DMM.Sstd)
    SL2  = np.array(DMM.SL2)
    TL2  = np.array(DMM.TL2)
    SRMS = np.array(DMM.SRMSE)
    TRMS = np.array(DMM.TRMSE)
    SINF = np.array(DMM.SINF)
    TINF = np.array(DMM.TINF)
    S01  = np.array(DMM.SS01, dtype=object)
    S02  = np.array(DMM.SS02, dtype=object)
    T01  = np.array(DMM.TT01, dtype=object)
    T02  = np.array(DMM.TT02, dtype=object)
    Srj  = np.array(DMM.Sqcrj)
    Trj  = np.array(DMM.Tqcrj)

    crc += 1
    if crc == 1:
      TSERR = tsarray(recn, numb, lon, lat, dnmb, rnmb, Tstd, Sstd, \
                   SL2, TL2, SRMS, TRMS, SINF, TINF, S01, T01, \
                   S02, T02, Srj, Trj) 
    else:
      TSERR.add_array(recn, numb, lon, lat, dnmb, rnmb, Tstd, Sstd, \
                   SL2, TL2, SRMS, TRMS, SINF, TINF, S01, T01, \
                   S02, T02, Srj, Trj) 


nlast = len(TSERR.numb)
print('Found {0} saved records'.format(nlast))
#
# Save as np arrays:
ffinal = pthbin+'TSprofALL_stat_'+rtofs_outp+'.pkl' 
print('Saving to '+ffinal)
with open(ffinal,'wb') as fid:
  pickle.dump(TSERR,fid)

print(' All done ' )

#aa = TSERR.SS01


