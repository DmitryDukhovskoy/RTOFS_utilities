"""
  Plot Argo and HYCOM T/S profile
  Specify Argo by # or lon/lat
  Argo profiles downloaded from HPSS QC log files
  see /scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/scripts/rtofs
  get_qctar.ssh

  ARGO time UTC, NCODA is UTC
  read below how to select rdate of RTOFS given Argo date/time

  For rtofs.20230122, the analysis time is 20230121, 
  and if the observation is new it should be included
  up to 2023012112.  
  So 20230121 10:31 should be included if it arrived before 20230122 
  For Argo profiles there is an observation time and an observation arrival time.

1) Check how well the profiles are being assimilated into RTOFS
Example for (1):
Compare Argo on 2023/1/28 6:41  vs RTOFS n-24 
with rdate ("forecast 00 date") = 20230129, 
i.e. ocean fields after the 6hr incremental update of the incup fields.
rtofs.20230129/rtofs_glo.t00z.n-24
background field: 
rtofs.20230128/rtofs_glo.t00z.f24

2) Check predictive skills of RTOFS in forward runs, i.e. (true) forecasts
plot TS profiles form 
rtofs.20230128/rtofs_glo.t00z.f06 (6-hr f/cast)
rtofs.20230127/rtofs_glo.t00z.f30  (30-hr f/cast)
rtofs.20230126/rtofs_glo.t00z.f54  (54-hr f/cast)
...
96-hr f/cast


# Derive HYCOM profiles:
# in rtofs.20230129/rtofs.ab.tar
# rtofs_glo.t00z.n-24.archv.a.tgz <- incrementally updated NCODA increments
# rtofs_glo.t00z.n00.archs.a.tgz  <- 24hr hindcast with atm analysis
# rtofs_glo.t00z.f12.archv.a.tgz  <- 12hr f/cast
# 
# in rtofs.20230129/rtofs_archv_1_inc.tar
# rtofs_glo.archv_1_inc.2023_028_00.a  <- NCODA increments + background
# rtofs_glo.incupd.2023_027_18.a       <- difference incup (inc fields - previous 
#                                                             18hr f/cast)
# rtofs.20230128/rtofs.ab.tar
# rtofs_glo.t00z.n00.archs.a.tgz  <--- background fields for 20230129 f/cast
#
# use script to untar from HPSS HYCOM fields:
# scripts/rtofs/get_rtofs_archv.sh
# ./get_rtofs_archv.sh 20230128 f24
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
#import pickle
#from netCDF4 import Dataset as ncFile
import matplotlib.colors as colors
import matplotlib.mlab as mlab
#from matplotlib.patches import Polygon
#from matplotlib.colors import ListedColormap
#from mpl_toolkits.basemap import Basemap, shiftgrid

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

# What RTOFS output field to plot:
# bkgrd  - background, f/cast from the previous day  n00 day before
# incup  - incr. updated NCODA increments from incup files (6hr run) -n24 
# fcst0  - state after 24hr hindcast ready for forecast n00 fcast day
# fcst12 - 12hr forecast 

rtofs_outp = 'incup'
#rtofs_outp = 'bkgrd'
rdate0     = '20230129'  # forecast date where Argo is assimilated

# Specify Argo to plot
# Note there could be duplicates and 
# also several profiles of the same
# Argo # 
# Use argoLT, argoLN to specify exact profile to plot 
# and make argoNmb < 0
#argoNmb = 6903111
#argoNmb = 4903493
#argoNmb = 6903239
argoNmb = 2902265
argoLT  = 16.69
argoLN  = 68.42 
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
pthhcm0 = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'
#pthhcmB = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+\
#          rdate_bkgrd+'/'

# Formatted text file
flnm = 'prof_argo_rpt.'+hdate_n24+'.txt'

if argoNmb > 0:
  TS = exargo.extract_argo(pthbin,flnm,Anmb=argoNmb)
else:
  TS = exargo.extract_argo(pthbin,flnm,Alat=argoLT,Alon=argoLN)

# dir(TS)  to see attributes
if len(TS.Tprof) == 0:
  raise Exception('Argo float was not found')


# Check time: ARGO profiles should be ~1 day back from the forecast rdate
# Saved is observation time, received time is Rcpt in txt ~2 hrs later
dnmb_argo = TS.ptime[0]
dstr_argo = mmisc.datestr(dnmb_argo)
print('Argo observed: '+mmisc.datestr(dnmb_argo))
rnmb_argo = TS.rtime[0]  # time of data recieved
rstr_argo = mmisc.datestr(rnmb_argo)
print('Argo Received: '+mmisc.datestr(rnmb_argo))
#print(TS.__dir__())  # to see attributes of TS

import mod_read_hycom as mhycom

IDM = 4500
JDM = 3298
KDM = 41

pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'
ftopo   = 'regional.depth'
fgrid   = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)


import mod_extrTS_hycom as exhcm
importlib.reload(exhcm)

if rtofs_outp == "incup":
  time_rtofs = time_n24     # actual time of the output field
  flhcm      = 'rtofs_glo.t00z.n-24.archv'
  pthhcm     = pthhcm0 + rdate_n24 + "/"
elif rtofs_outp == "bkgrd":
  time_rtofs = time_bkgrd     # actual time of the output field
  flhcm      = 'rtofs_glo.t00z.n00.archv'
  pthhcm     = pthhcm0 + rdate_bkgrd + '/'
 
fld   = 'temp'
latar = TS.lat[0]
lonar = TS.lon[0]

Th, ZZh, ZMh, dHh = exhcm.extract_1prof(pthhcm,flhcm,'temp',lonar,latar,LON,LAT)
Sh, dm1, dm2, dm3 = exhcm.extract_1prof(pthhcm,flhcm,'salin',lonar,latar,\
                                        LON,LAT,fthknss=False)

# Selected profiles  
import plot_sect as psct
importlib.reload(psct)
#
# Argo profiles:
ZSar = np.array(TS.ZS).squeeze()
Sar  = np.array(TS.Sprof).squeeze()
ZTar = np.array(TS.ZT).squeeze()
Tar  = np.array(TS.Tprof).squeeze()
#psct.plot_2prof(ZSar,Sar,ZTar,Tar) 

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

# Define T, S limits:
zz0  = -400.
nInt = 5
th1, th2, dTh = mutil.prof_limits(Th, ZMh, nT=nInt, z0=zz0)
ta1, ta2, dTa = mutil.prof_limits(iTar, ZMh, nT=nInt, z0=zz0)
tlim1 = min([th1,ta1])
tlim2 = max([th2,ta2])
#dltT  = min([dTh,dTa])
dltT  = 0.1*(int((tlim2-tlim1)*10/nInt))

sh1, sh2, dSh = mutil.prof_limits(Sh,ZMh,z0=zz0)
sa1, sa2, dSa = mutil.prof_limits(iSar,ZMh,z0=zz0)
slim1 = min([sh1,sa1])
slim2 = max([sh2,sa2])
#dltS  = min([dSh,dSa])
dltS  = 0.1*(int((slim2-slim1)*10/nInt))


# Error:
Sl2, Srmse, Sinf = mutil.err_stat(Sh, iSar)
Tl2, Trmse, Tinf = mutil.err_stat(Th, iTar)


clra  = [0.,0.4,0.8]  # argo interpolated profile
clrh  = [0.9,0.3,0]
clra0 = [0.,0.9,0.4]  # original profile

ctl1 = '$T^oC$ '+rtofs_outp
ctl2 = 'S '+rtofs_outp

# =====================
#   Plotting:
# =====================
f_argorg = True  # plot original Argo profiles

plt.ion()
fig1 = plt.figure(1,figsize=(9,8), constrained_layout=False)
fig1.clf()
#
# T profile
ax1 = plt.axes([0.08, 0.1, 0.25, 0.8]) 
plt.plot(iTar,ZMh, '.-', color=clra, linewidth=2, label="Argo")
plt.plot(Th, ZMh, '.-', color=clrh, linewidth=2, label="RTOFS")
if f_argorg:
  plt.plot(Tar, ZTar, '.-', color=clra0, linewidth=2, label="Argo orig")
plt.yticks(np.arange(-2000,0,50))
plt.xticks(np.arange(tlim1,tlim2,dltT))
plt.ylim(zz0,0)
plt.xlim(tlim1,tlim2)
ax1.grid(True)
ax1.set_title(ctl1)

# S profile
ax2 = plt.axes([0.4, 0.1, 0.25, 0.8]) 
#plt.plot(Sar, ZSar, '.-', color=clra, linewidth=2, label="Argo")
ln1, = plt.plot(iSar, ZMh, '.-', color=clra, linewidth=2, label="Argo")
ln2, = plt.plot(Sh, ZMh, '.-', color=clrh, linewidth=2, label="RTOFS")
if f_argorg:
  ln3, = plt.plot(Sar, ZSar, '.-', color=clra0, linewidth=2, label="Argo orig")
plt.yticks(np.arange(-2000,0,50))
plt.xticks(np.arange(slim1,slim2,dltS))
plt.ylim(zz0,0)
plt.xlim(slim1,slim2)
ax2.grid(True)
ax2.set_title(ctl2)

# Plot map:
LMSK = HH.copy()
LMSK = np.where(LMSK<0.,0.,1.)
lcmp = mutil.clrmp_lmask(2, clr_land=[0.5,0.5,0.5])
ax3 = plt.axes([0.7, 0.4, 0.25, 0.25])
ax3.pcolormesh(LMSK, shading='flat',\
                cmap=lcmp, \
                vmin=0, \
                vmax=1)
ax3.contour(HH, [-5000,-4000,-3000,-2000,-1000], \
            colors=[(0.9,0.9,0.9)], linestyles='solid',linewidths=1)
ax3.plot(Iex, Jex, marker='.', markersize=12, color=[1.,0.2,0])
if lonar<0:
  wfx = 'W'
else:
  wfx = 'E'
if latar<0:
  nfx = 'S'
else:
  nfx = 'N'

dii = 250   # +/- dii map range
sst = '{0:5.2f}{1}, {2:5.2f}{3}'.format(lonar,wfx,latar,nfx)
ax3.text(Iex-int(dii/7),Jex+int(dii/10),sst)

ilim1 = max([Iex-dii,0])
ilim2 = min([Iex+dii,IDM])
jlim1 = max([Jex-dii,0])
jlim2 = min([Jex+dii,JDM])

ax3.set_xlim([ilim1,ilim2])
ax3.set_ylim([jlim1,jlim2])

# Set y, x axis not visible
#ahd = plt.gca()
xax = ax3.axes.get_xaxis()
xax = xax.set_visible(False)
yax = ax3.axes.get_yaxis()
yax = yax.set_visible(False)

# Info:
Sstd    = TS.Sstd[0]
Tstd    = TS.Tstd[0]
frjct_s = 'accpt'
frjct_t = 'accpt'
if Sstd > Srjct:
  frjct_s = 'reject'
if Tstd > Trjct:
  frjct_t = 'reject'

dv  = mmisc.datevec(dnmb_argo)
dvr = mmisc.datevec(rnmb_argo)

ax4 = plt.axes([0.66, 0.7, 0.25, 0.25])
ss1 = 'ARGO #{0}\n'.format(TS.numb[0])
ss1 = ss1 + 'ARGO   date: {0}/{1}/{2} {3:02d}:{4:02d} UTC\n'.\
            format(dv[0],dv[1],dv[2],dv[3],dv[4])
ss1 = ss1 + 'ARGO receiv: {0}/{1}/{2} {3:02d}:{4:02d} UTC\n'.\
            format(dvr[0],dvr[1],dvr[2],dvr[3],dvr[4])
astr = flhcm.split(".")
YY = time_rtofs.year
MM = time_rtofs.month
DM = time_rtofs.day
HR = time_rtofs.hour

ss1 = ss1 + 'RTOFS date: {0}/{1}/{2} {3:02d}:00 UTC\n'.format(YY,MM,DM,HR)
ss1 = ss1 + 'RTOFS f/cast date : {0}/{1}/{2}\n'.\
          format(rdate0[0:4],rdate0[4:6],rdate0[6:8])
ss1 = ss1 + 'RTOFS h/cast date : {0}/{1}/{2}\n'.\
          format(rdate_bkgrd[0:4],rdate_bkgrd[4:6],rdate_bkgrd[6:8])
amm = pthhcm.split("/")
ss1 = ss1 + 'RTOFS output:\n {0}/{1}\n'.format(amm[-2],flhcm)
ss1 = ss1 + 'T QC STD = {0} {1}\n'.format(Tstd,frjct_t)
ss1 = ss1 + 'S QC STD = {0} {1}\n'.format(Sstd,frjct_s)

ax4.text(0.0,0.1,ss1)
ax4.axis('off')


# Error:
ss2 = 'Error Stat    T             S:\n'
ss2 = ss2+'L2 norm : {0:9.7f}  {1:9.7f}\n'.format(Tl2,Sl2)
ss2 = ss2+'RMSE    : {0:9.7f}  {1:9.7f}\n'.format(Trmse,Srmse)
ss2 = ss2+'Inf norm: {0:9.7f}  {1:9.7f}\n'.format(Tinf,Sinf)


ax5 = plt.axes([0.7, 0.1, 0.25, 0.25])
if f_argorg:
  lgd = plt.legend(handles=[ln1,ln2,ln3], loc='upper left')
else:
  lgd = plt.legend(handles=[ln1,ln2], loc='upper left')
ax5.text(0,0.02,ss2)
ax5.axis('off')


try:
  rFile = __file__
  print('rFile = '+rFile)
  dmm = rFile.split("/")
  btx = dmm[-1]
except:
  btx = 'tsprof_argo_hycom.py'

bottom_text(btx,pos=[0.02, 0.03])

plt.show()



