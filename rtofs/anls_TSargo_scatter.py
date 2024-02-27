"""
  Analyze error stat from combined pkl file
  i.e. after running combine_stat.py and 
  after running serial tsprof_error_derive.py 
  for a set of Argo profiles

  Scatter plots of T/S

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
# Unpack T, S binned by depth intervals
# HYCOM
dmm = TSERR.SS01
SHb1 = mutil.unpack_TSbins(dmm,0)
SHb2 = mutil.unpack_TSbins(dmm,1)
dmm = TSERR.TT01.copy()
THb1 = mutil.unpack_TSbins(dmm,0)
THb2 = mutil.unpack_TSbins(dmm,1)
# Argo
dmm = TSERR.SS02
SAb1 = mutil.unpack_TSbins(dmm,0)
SAb2 = mutil.unpack_TSbins(dmm,1)
dmm = TSERR.TT02
TAb1 = mutil.unpack_TSbins(dmm,0)
TAb2 = mutil.unpack_TSbins(dmm,1)

# Find NCODA consecutive #s of rejected profiles
Srj  = TSERR.Sqcrj.astype(int)
Trj  = TSERR.Tqcrj.astype(int)
rnmb = TSERR.recn.astype(int)
Srj = Srj[np.where(Srj >= 0)]
Trj = Trj[np.where(Trj >= 0)]
ns  = Srj.shape[0]
nt  = Trj.shape[0]

# Indices of T or S rejected profile - total
TSrj = Srj.copy()
for ii in range(nt):
  smm = abs(Srj-Trj[ii])
  if min(smm) > 0:
    TSrj = np.append(TSrj,Trj[ii])

# Find indices of NCODA prof nmbers (1, 2, ....) for rejected T/S profiles
# in the current order of profiles (unordered) <-- obsolete
# should be ordered now
#nts = TSrj.shape[0]
#Irj = np.array([], dtype=int)
#for ii in range(nts):
#  irj = np.where(rnmb == TSrj[ii])
#  Irj = np.append(Irj,irj)
#Irj = TSrj

#Add algorithm for identifying binned values 
#of rejected profiles
# HYCOM
dmm = TSERR.SS01
SHb1_trj = mutil.find_rejected_TSbinned(Trj,dmm,0)
SHb2_trj = mutil.find_rejected_TSbinned(Trj,dmm,1)
SHb1_srj = mutil.find_rejected_TSbinned(Srj,dmm,0)
SHb2_srj = mutil.find_rejected_TSbinned(Srj,dmm,1)
dmm = TSERR.TT01
THb1_trj = mutil.find_rejected_TSbinned(Trj,dmm,0)
THb2_trj = mutil.find_rejected_TSbinned(Trj,dmm,1)
THb1_srj = mutil.find_rejected_TSbinned(Srj,dmm,0)
THb2_srj = mutil.find_rejected_TSbinned(Srj,dmm,1)
## Argo
dmm = TSERR.SS02
SAb1_trj = mutil.find_rejected_TSbinned(Trj,dmm,0)
SAb2_trj = mutil.find_rejected_TSbinned(Trj,dmm,1)
SAb1_srj = mutil.find_rejected_TSbinned(Srj,dmm,0)
SAb2_srj = mutil.find_rejected_TSbinned(Srj,dmm,1)
dmm = TSERR.TT02
TAb1_trj = mutil.find_rejected_TSbinned(Trj,dmm,0)
TAb2_trj = mutil.find_rejected_TSbinned(Trj,dmm,1)
TAb1_srj = mutil.find_rejected_TSbinned(Srj,dmm,0)
TAb2_srj = mutil.find_rejected_TSbinned(Srj,dmm,1)

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
  
# Find profiles rec# for specific quest. values:
# not rejected profiles
vmin  = 27.9
vmax  = 28.2
ibin  = 0
dmm = TSERR.SS02  # Argo S binned values
vindx = mutil.find_val_TSbinned(dmm,vmin,vmax,ibin,TSrj)

# Time receipt of ob data:
Trcpt = TSERR.rtime

# ==================
#   Plotting
# ==================
print('Plotting ...')
plt.ion()

clrtrj = [0.9,0.6,0.0]  # Rejected T data
clrsrj = [0.,0.9,0]    # rejected S data
t1l1, t1l2 = mutil.axis_limits2(THb1,TAb1,cff=0.1)
s1l1, s1l2 = mutil.axis_limits2(SHb1,SAb1,cff=0.1)
dltT = 5
dltS = 2

t2l1, t2l2 = mutil.axis_limits2(THb2,TAb2,cff=0.1)
s2l1, s2l2 = mutil.axis_limits2(SHb2,SAb2,cff=0.1)

fig1 = plt.figure(1,figsize=(9,9), constrained_layout=False)
fig1.clf()
# T profiles, shallow bin = 1
ax1 = plt.axes([0.07,0.5,0.3,0.3])
ax1.plot(THb1,TAb1,'.')
# Rejected
if len(Trj)>0:
  for ii in range(len(Trj)):
    th = THb1_trj[ii]
    ta = TAb1_trj[ii]
    ax1.plot(th,ta,'.',color=clrtrj)
if len(Srj)>0:
  for ii in range(len(Srj)):
    th = THb1_srj[ii]
    ta = TAb1_srj[ii]
    ax1.plot(th,ta,'.',color=clrsrj)

ax1.plot([t1l1,t1l2],[t1l1,t1l2],'--',color=[1,0,0])
ax1.set_xticks(np.arange(-5,40,dltT))
ax1.set_yticks(np.arange(-5,40,dltT))
ax1.set_xlabel('HYCOM')
ax1.set_ylabel('Observ')
ax1.set_xlim(t1l1,t1l2)
ax1.set_ylim(t1l1,t1l2)
ax1.grid(True)
ctl = 'T RTOFS {0} vs Argo {1}\n anls={2}'.\
           format(rtofs_outp,sbin1,rdate_anls)
ax1.set_title(ctl)

# S profiles, shallow bin
ax2 = plt.axes([0.47,0.5,0.3,0.3])
ax2.plot(SHb1,SAb1,'.')
# Rejected:
if len(Trj)>0:
  for ii in range(len(Trj)):
    th = SHb1_trj[ii]
    ta = SAb1_trj[ii]
    ax2.plot(th,ta,'.',color=clrtrj)
if len(Srj)>0:
  for ii in range(len(Srj)):
    th = SHb1_srj[ii]
    ta = SAb1_srj[ii]
    ax2.plot(th,ta,'.',color=clrsrj)

#ax2.plot(SHb1[Irj],SAb1[Irj],'.',color=clrj)
ax2.plot([s1l1,s1l2],[s1l1,s1l2],'--',color=[1,0,0])
ax2.set_xticks(np.arange(4,48,dltS))
ax2.set_yticks(np.arange(4,48,dltS))
ax2.set_xlabel('HYCOM')
ax2.set_ylabel('Observ')
ax2.set_xlim(s1l1,s1l2)
ax2.set_ylim(s1l1,s1l2)
ax2.grid(True)
ctl = 'S RTOFS {0} vs Argo {1}\n anls={2}'.\
           format(rtofs_outp,sbin1,rdate_anls)
ax2.set_title(ctl)


# T profiles, deep bin
ax3 = plt.axes([0.07,0.1,0.3,0.3])
ax3.plot(THb2,TAb2,'.')
# Rejected:
if len(Trj)>0:
  for ii in range(len(Trj)):
    th = THb2_trj[ii]
    ta = TAb2_trj[ii]
    ax3.plot(th,ta,'.',color=clrtrj)
if len(Srj)>0:
  for ii in range(len(Srj)):
    th = THb2_srj[ii]
    ta = TAb2_srj[ii]
    ax3.plot(th,ta,'.',color=clrsrj)

ax3.plot([t2l1,t2l2],[t2l1,t2l2],'--',color=[1,0,0])
ax3.set_xticks(np.arange(-5,40,dltT))
ax3.set_yticks(np.arange(-5,40,dltT))
ax3.set_xlabel('HYCOM')
ax3.set_ylabel('Observ')
ax3.set_xlim(t2l1,t2l2)
ax3.set_ylim(t2l1,t2l2)
ax3.grid(True)
ctl = 'T RTOFS {0} vs Argo {1}\n anls={2}'.\
           format(rtofs_outp,sbin2,rdate_anls)
ax3.set_title(ctl)

# S profiles, deep bin
ax4 = plt.axes([0.47,0.1,0.3,0.3])
ax4.plot(SHb2,SAb2,'.')
# Rejected:
if len(Trj)>0:
  for ii in range(len(Trj)):
    th = SHb2_trj[ii]
    ta = SAb2_trj[ii]
    ax3.plot(th,ta,'.',color=clrtrj)
if len(Srj)>0:
  for ii in range(len(Srj)):
    th = SHb2_srj[ii]
    ta = SAb2_srj[ii]
    ax4.plot(th,ta,'.',color=clrsrj)
ax4.plot([s2l1,s2l2],[s2l1,s2l2],'--',color=[1,0,0])
ax4.set_xticks(np.arange(4,48,dltS))
ax4.set_yticks(np.arange(4,48,dltS))
ax4.set_xlabel('HYCOM')
ax4.set_ylabel('Observ')
ax4.set_xlim(s1l1,s1l2)
ax4.set_ylim(s1l1,s1l2)
ax4.grid(True)
ctl = 'S RTOFS {0} vs Argo {1}\n anls={2}'.\
           format(rtofs_outp,sbin2,rdate_anls)
ax4.set_title(ctl)

# Information - rejected profiles by lat:
ax5 = plt.axes([0.82,0.7,0.15,0.1])
ax5.hist(Lat_rjct,Lbin, histtype='bar', rwidth=0.95, color=([0.4,0.4,0.4]))
ax5.set_xticks(Lbin)
ax5.set_xlim(-80,90)
ax5.set_xticklabels(['','','50S','','10S','','30N','','70N',''])
ax5.set_title('Rejected profiles \nby lat')

#
# Information:
ss1 = 'RTOFS f/cast date: \n{0}/{1}/{2}\n'.\
       format(rdate0[0:4],rdate0[4:6],rdate0[6:8])
ss1 = ss1 + 'RTOFS {0}/Argo T/S \n'.format(rtofs_outp)
ss1 = ss1 + 'T-rejected prof = {0}\n'.format(Trj.shape[0])
ss1 = ss1 + 'S-rejected prof = {0}\n'.format(Srj.shape[0])
ss1 = ss1 + 'Total rejected  = {0}\n'.format(TSrj.shape[0])

ax6 = plt.axes([0.8,0.5,0.19,0.3])
ax6.text(0.0,0.1,ss1)
ax6.axis('off')

ax7 = plt.axes([0.82,0.4,0.12,0.05])
ax7.plot(0.,0.1, '.', color=clrtrj)
ax7.plot(0.,0.7,'.', color=clrsrj)
ax7.text(0.2,-0.2,'T-rejected')
ax7.text(0.2,0.3,'S-rejected')
ax7.set_xlim(-0.5,2)
ax7.set_ylim(-0.5,1.3)
ax7.set_xticks([])
ax7.set_yticks([])


try:
  rFile = __file__
  print('rFile = '+rFile)
  dmm = rFile.split("/")
  btx = dmm[-1]
except:
  btx = 'tsprof_argo_hycom.py'

bottom_text(btx,pos=[0.02, 0.03])

plt.show()










