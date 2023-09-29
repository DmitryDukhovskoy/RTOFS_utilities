"""
  Analyze error stat from combined pkl file
  i.e. after running combine_stat.py and 
  after running serial tsprof_error_derive.py 
  for a set of Argo profiles

  Histograms of errors
  by lat / depth bins

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
from mod_utils import tsarray

np.set_printoptions(precision=3)

# What RTOFS output field to plot:
# bkgrd  - background, f/cast from the previous day  n00 day before
# incup  - incr. updated NCODA increments from incup files (6hr run) -n24 
# fcst0  - state after 24hr hindcast ready for forecast n00 fcast day
# fcst12 - 12hr forecast 

rtofs_outp = 'incup'
#rtofs_outp = 'bkgrd'
rdate0     = '20230130'  # forecast date with Argo is assimilated

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

pthbin = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+\
          rdate0+'/ocnqc_logs/profile_qc/'
ffinal = pthbin+'TSargo_stat_'+rdate_anls+rtofs_outp+'.pkl' 

print('Loading '+ffinal)
with open(ffinal,'rb') as fid:
  TSERR = pickle.load(fid)

Zbins = np.array([0,-200,-6000])
sbin1 = '{0}-{1} m'.format(abs(Zbins[0]),abs(Zbins[1]))
sbin2 = '{0}-{1} m'.format(abs(Zbins[1]),abs(Zbins[2]))

# Find NCODA consecutive #s of rejected/not rejected profiles
Srj   = TSERR.Sqcrj.astype(int)
Trj   = TSERR.Tqcrj.astype(int)
rnmb  = TSERR.recn.astype(int)
# Not rejected indices:
Snrj  = Srj[np.where(Srj < 0)]
Tnrj  = Trj[np.where(Trj < 0)]
TSnrj = np.where( (Srj < 0) & (Trj < 0) )[0]
# Rejected indices:
TSrj  = np.where( (Srj >= 0) | (Trj >= 0) )[0]
Srj   = Srj[np.where(Srj >= 0)]
Trj   = Trj[np.where(Trj >= 0)]
ns    = Srj.shape[0]
nt    = Trj.shape[0]

# Unpack T, S binned by depth intervals
# HYCOM
dmm = TSERR.SS01
SHb1 = mutil.unpack_TSbins_Iprf(dmm,TSnrj,0)
SHb2 = mutil.unpack_TSbins_Iprf(dmm,TSnrj,1)
dmm = TSERR.TT01
THb1 = mutil.unpack_TSbins_Iprf(dmm,TSnrj,0)
THb2 = mutil.unpack_TSbins_Iprf(dmm,TSnrj,1)
# Argo
dmm = TSERR.SS02
SAb1 = mutil.unpack_TSbins_Iprf(dmm,TSnrj,0)
SAb2 = mutil.unpack_TSbins_Iprf(dmm,TSnrj,1)
dmm = TSERR.TT02
TAb1 = mutil.unpack_TSbins_Iprf(dmm,TSnrj,0)
TAb2 = mutil.unpack_TSbins_Iprf(dmm,TSnrj,1)


#
# Argo numbers:
Anmb = TSERR.numb
# Rejected Argo numbers:
Anmb_rj = Anmb[TSrj]
# Histogram of 
# Rejected by lat:
Lat      = TSERR.lat
Lbin     = np.arange(-90,110,20)
Nbins    = Lbin.shape[0]
Lat_rjct = Lat[TSrj]
Lat_hist, bins = np.histogram(Lat_rjct,Lbin)
  
# Find profiles rec# for specified values:
# not rejected profiles
vmin  = 27.9
vmax  = 28.2
ibin  = 0
dmm = TSERR.SS02  # Argo S binned values
vindx = mutil.find_val_TSbinned(dmm,vmin,vmax,ibin,TSrj)

# Time receipt of ob data:
Trcpt = TSERR.rtime

#
# Calc difference for depth binsi, global:
dTb1 = THb1-TAb1
dTb2 = THb2-TAb2
dSb1 = SHb1-SAb1
dSb2 = SHb2-SAb2

# Calc difference by lat for reject & notrejected profiles:
for kk in range(Nbins-1):
  lt1  = Lbin[kk]
  lt2  = Lbin[kk+1]
  Ilat = np.where((Lat > lt1) & (Lat <= lt2))[0]
  if len(Ilat) == 0:
    continue

# HYCOM S
  dmm = TSERR.SS01
  SHb1 = mutil.unpack_TSbins_Iprf(dmm,Ilat,0)
  SHb2 = mutil.unpack_TSbins_Iprf(dmm,Ilat,1)
  dmm = TSERR.TT01
  THb1 = mutil.unpack_TSbins_Iprf(dmm,Ilat,0)
  THb2 = mutil.unpack_TSbins_Iprf(dmm,Ilat,1)
  # Argo
  dmm = TSERR.SS02
  SAb1 = mutil.unpack_TSbins_Iprf(dmm,Ilat,0)
  SAb2 = mutil.unpack_TSbins_Iprf(dmm,Ilat,1)
  dmm = TSERR.TT02
  TAb1 = mutil.unpack_TSbins_Iprf(dmm,Ilat,0)
  TAb2 = mutil.unpack_TSbins_Iprf(dmm,Ilat,1)

  dS1 = SHb1-SAb1
  dS2 = SHb2-SAb2
  dT1 = THb1-TAb1
  dT2 = THb2-TAb2

# S depth bin 1
  medS1   = np.nanmedian(dS1)
  iqrbS1  = np.nanpercentile(dS1,25)
  iqruS1  = np.nanpercentile(dS1,75)
  p10S1   = np.nanpercentile(dS1,10)
  p90S1   = np.nanpercentile(dS1,90)
  vmin    = np.nanmin(dS1)
  vmax    = np.nanmax(dS1)
  if kk == 0:
    STS1 = mutil.med_iqr(medS1, iqrbS1, iqruS1, p10S1, p90S1, vmin, vmax)
  else:
    STS1.add_data(medS1, iqrbS1, iqruS1, p10S1, p90S1, vmin, vmax)
# S depth bin 2
  medS2   = np.nanmedian(dS2)
  iqrbS2  = np.nanpercentile(dS2,25)
  iqruS2  = np.nanpercentile(dS2,75)
  p10S1   = np.nanpercentile(dS2,10)
  p90S1   = np.nanpercentile(dS2,90)
  vmin    = np.nanmin(dS2)
  vmax    = np.nanmax(dS2)
  if kk == 0:
    STS2 = mutil.med_iqr(medS2, iqrbS2, iqruS2, p10S1, p90S1, vmin, vmax)
  else:
    STS2.add_data(medS2, iqrbS2, iqruS2, p10S1, p90S1, vmin, vmax)

# T depth bin 1
  medT1   = np.nanmedian(dT1)
  iqrbT1  = np.nanpercentile(dT1,25)
  iqruT1  = np.nanpercentile(dT1,75)
  p10T1   = np.nanpercentile(dT1,10)
  p90T1   = np.nanpercentile(dT1,90)
  vmin    = np.nanmin(dT1)
  vmax    = np.nanmax(dT1)
  if kk == 0:
    STT1 = mutil.med_iqr(medT1, iqrbT1, iqruT1, p10T1, p90T1, vmin, vmax)
  else:
    STT1.add_data(medT1, iqrbT1, iqruT1, p10T1, p90T1, vmin, vmax)

# T depth bin 2:
  medT2   = np.nanmedian(dT2)
  iqrbT2  = np.nanpercentile(dT2,25)
  iqruT2  = np.nanpercentile(dT2,75)
  p10T2   = np.nanpercentile(dT2,10)
  p90T2   = np.nanpercentile(dT2,90)
  vmin    = np.nanmin(dT2)
  vmax    = np.nanmax(dT2)
  if kk == 0:
    STT2 = mutil.med_iqr(medT2, iqrbT2, iqruT2, p10T2, p90T2, vmin, vmax)
  else:
    STT2.add_data(medT2, iqrbT2, iqruT2, p10T2, p90T2, vmin, vmax)



# ==================
#   Plotting
# ==================
print('Plotting ...')
plt.ion()

#t2l1, t2l2 = mutil.axis_limits2(THb2,TAb2,cff=0.1)
#s2l1, s2l2 = mutil.axis_limits2(SHb2,SAb2,cff=0.1)

# List of xaxis labels:
Xlbls = []
for kk in range(Nbins-1):
  lt1  = Lbin[kk]
  lt2  = Lbin[kk+1]
  if lt1 < 0:
    ns1 = 'S'
  else:
    ns1 = 'N'

  if lt2 < 0:
    ns2 = 'S'
  else:
    ns2 = 'N'

  sstr = '{0}{1}-\n{2}{3}'.format(abs(lt1),ns1,abs(lt2),ns2)
  Xlbls.append(sstr)

fig1 = plt.figure(1,figsize=(9,9), constrained_layout=False)
fig1.clf()
# S profiles, shallow bin = 1
ax1  = plt.axes([0.07,0.51,0.4,0.3])
ctl = 'dS (RTOFS - Argo), {0}\n anls={1}'.\
           format(sbin1,rdate_anls)
ax1 = mutil.plot_boxplot(ax1, STS1, ctl=ctl, Xlbls=Xlbls)

# dlt T upper bin
ax2  = plt.axes([0.55,0.51,0.4,0.3])
ctl = 'dT (RTOFS - Argo), {0}\n anls={1}'.\
           format(sbin1,rdate_anls)
ax2 = mutil.plot_boxplot(ax2, STT1, ctl=ctl, Xlbls=Xlbls)

# dlt S deep bin
ax3  = plt.axes([0.07,0.1,0.4,0.3])
ctl = 'dS (RTOFS - Argo), {0}\n anls={1}'.\
           format(sbin2,rdate_anls)
ax3 = mutil.plot_boxplot(ax3, STS2, ctl=ctl, Xlbls=Xlbls)

# dlt T deep bin
ax4  = plt.axes([0.55,0.1,0.4,0.3])
ctl = 'dT (RTOFS - Argo), {0}\n anls={1}'.\
           format(sbin2,rdate_anls)
ax4 = mutil.plot_boxplot(ax4, STT2, ctl=ctl, Xlbls=Xlbls)


# Information:
ss1 = 'RTOFS f/cast date: {0}/{1}/{2}\n'.\
       format(rdate0[0:4],rdate0[4:6],rdate0[6:8])
ss1 = ss1 + 'RTOFS {0}/Argo T/S \n'.format(rtofs_outp)

ax5 = plt.axes([0.1,0.86,0.8,0.12])
ax5.text(0.0,0.1,ss1)
ax5.axis('off')



try:
  rFile = __file__
  print('rFile = '+rFile)
  dmm = rFile.split("/")
  btx = dmm[-1]
except:
  btx = 'anls_TSargo_errbxplt.py'

bottom_text(btx,pos=[0.02, 0.02])

plt.show()










