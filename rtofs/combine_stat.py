"""
  Combine pkl error stat files into a single file
  after running serial tsprof_error_anls.py 
  for a set of Argo profiles

  for parallel runs on N nodes 
  use 
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

import mod_utils as mutil
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_read_ncoda as rncoda
#importlib.reload(rncoda)
import mod_extrargo as exargo
importlib.reload(exargo)
from mod_utils import tsarray 
from mod_utils import tserr

# What RTOFS output field to plot:
# bkgrd  - background, f/cast from the previous day  n00 day before
# incup  - incr. updated NCODA increments from incup files (6hr run) -n24 
# fcst0  - state after 24hr hindcast ready for forecast n00 fcast day
# fcst12 - 12hr forecast 

rtofs_outp = 'incup'
#rtofs_outp = 'bkgrd'
rdate0     = '20230130'  # forecast date where Argo is assimilated

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
rdate_anls  = rdate_bkgrd    # analysis time
# 12hr f/cast:
time_f12    = datetime.datetime(yr,mo,mday,hr,12,0)
hdate_f12   = time_f12.strftime('%Y%m%d%H')
rdate_f12   = rdate0


#class tserr():
# kind = 'TS Argo error stat'

# def __init__(self):
#   self.recn  = []
#   self.numb  = []
#   self.lon   = []
#   self.lat   = []
#   self.ptime = []
#   self.rtime = []
#   self.Tstd  = []
#   self.Sstd  = []
#   self.SL2   = []
#   self.TL2   = []
#   self.SRMSE = []
#   self.TRMSE = []
#   self.SINF  = []
#   self.TINF  = []
#   self.SS01  = []  # pooled S values in the top level 01
#   self.TT01  = []
#   self.SS02  = []
#   self.TT02  = []
#   self.Sqcrj = []  # indices of rejected S profiles
#   self.Tqcrj = []  # omdoces pf rejected T profiles

# def add_data(self, recn, numb, lon, lat, dnmb, rnmb, Tstd, Sstd, \
#              SL2, TL2, SRMSE, TRMSE, SINF, TINF, SS01, TT01, \
#              SS02, TT02, Sqcrj, Tqcrj):
#   self.recn.append(recn)    # record #
#   self.numb.append(numb)    # Argo #
#   self.lon.append(lon)
#   self.lat.append(lat)
#   self.ptime.append(dnmb)
#   self.rtime.append(rnmb)
#   self.Tstd.append(Tstd)
#   self.Sstd.append(Sstd)
#   self.SL2.append(SL2)
#   self.TL2.append(TL2)
#   self.SRMSE.append(SRMSE)
#   self.TRMSE.append(TRMSE)
#   self.SINF.append(SINF)
#   self.TINF.append(TINF)
#   self.SS01.append(SS01)     #values in the top level 01
#   self.TT01.append(TT01)
#   self.SS02.append(SS02)
#   self.TT02.append(TT02)
#   self.Sqcrj.append(Sqcrj) #S rejected profiles
#   self.Tqcrj.append(Tqcrj)

#
# Save as np array
#TSERR   = tsarray()
crc     = 0
pthbin  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+\
          rdate0+'/ocnqc_logs/profile_qc/'
pthout  = pthbin
flbase  = 'TSprof_stat_'
lnnm    = len(flbase)
# Read error stat files in order 
iFmin = 100
iFmax = 0
for fll in os.listdir(pthbin):
  ffpth = os.path.join(pthbin,fll)
  if os.path.isfile(ffpth) and fll[:lnnm]==flbase:
    print(fll)
    smm   = fll.split('_')[2]
    ii1   = int(smm[:4])
    ii2   = int(smm[5:9])
    dFnmb = ii2-ii1+1
    iFmin = min([ii1,iFmin])
    iFmax = max([ii2,iFmax])

print('Reading files from # {0}:{2}:{1}'.format(iFmin,iFmax,dFnmb))
for inmb in range(iFmin,iFmax,dFnmb):
  ii1 = inmb
  ii2 = inmb+dFnmb-1
  if ii2 > iFmax:
    ii2 = iFmax
  fnmb = '{0:04d}-{1:04d}'.format(ii1,ii2)
  fll = flbase+fnmb+'.pkl'
 
  ffpth = os.path.join(pthbin,fll)
  if os.path.isfile(ffpth) and fll[:lnnm]==flbase:
    print('Reading '+fll)

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
ffinal = pthbin+'TSargo_stat_'+rdate_anls+rtofs_outp+'.pkl' 
print('Saving to '+ffinal)
with open(ffinal,'wb') as fid:
  pickle.dump(TSERR,fid)

print(' All done ' )

#aa = TSERR.SS01


