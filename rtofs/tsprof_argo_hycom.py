"""
  Plot Argo and HYCOM T/S profile
  Specify Argo by # or lon/lat
  Argo profiles downloaded from HPSS QC log files
  see /scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/scripts/rtofs
  get_qctar.ssh
  For Argo profile on date MM/DD 
  need rtofs.YYMM(DD+1)

  ARGO time:

  Generally Argo arrives soon after observation. Other data arrive later.
  I guess yes all times are UTC, as is ncoda's analysis time.

  For rtofs.20230122, the analysis time is 20230121, 
  and if the observation is new it should be included
  up to 2023012112.  
  So yes 202301211031 should be included if it arrived before 20230122 
  (plus maybe 1-2 hours, production
  launching the daily run), 
  This is if the observation has not have been received before, which for that 
  day/time is not possible,
  and if you are using the Argo txt files those duplicate checks are already done. 
  There is an observation time and an observation arrival time.

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




rdate   = '20230129'
pthbin  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+\
          rdate+'/ocnqc_logs/profile_qc/'
pthhcm  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+rdate+'/'
# 
# Specify Argo to plot
argoNmb = 6903111
argoLT  = 14.41
argoLN  = -61.42 

yr, mo, mday, hr = rncoda.parse_rdate(rdate)
time_n00  = datetime.datetime(yr,mo,mday,hr,0,0)  # forecast time
time_n24  = time_n00 - datetime.timedelta(days=1) # time of hindcast incup incorpor
rdate_n24 = time_n24.strftime('%Y%m%d%H')

# Formatted text file
flnm = 'prof_argo_rpt.'+rdate_n24+'.txt'

TS = exargo.extract_argo(pthbin,flnm,Anmb=argoNmb)
# dir(TS)  to see attributes
if len(TS.Tprof) == 0:
  raise Exception('Argo float was not found')


# Check time: ARGO profiles should be ~1 day back from the forecast rdate
# Saved is observation time, received time is Rcpt in txt ~2 hrs later
dnmb_argo = TS.ptime[0]
dstr_argo = mmisc.datestr(dnmb_argo)
print(mmisc.datestr(dnmb_argo))


#
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
# Previous f/cast - used as a background:
# rtofs/
#
# use script to untar from HPSS HYCOM fields:
# scripts/rtofs/get_rtofs_archv.sh
#

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

flhcm = 'rtofs_glo.t00z.n-24.archv'
fld   = 'temp'
latar = TS.lat[0]
lonar = TS.lon[0]

Th, ZZh, ZMh, dHh = exhcm.extract_1prof(pthhcm,flhcm,'temp',lonar,latar,LON,LAT)
Sh, dm1, dm2, dm3 = exhcm.extract_1prof(pthhcm,flhcm,'salin',lonar,latar,\
                                        LON,LAT,fthknss=False)

# Selected profiles  
import plot_sect as psct
importlib.reload(psct)

ZSar = np.array(TS.ZS).squeeze()
Sar  = np.array(TS.Sprof).squeeze()
ZTar = np.array(TS.ZT).squeeze()
Tar  = np.array(TS.Tprof).squeeze()
#psct.plot_2prof(ZSar,Sar,ZTar,Tar) 

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

btx = 'tsprof_argo_hycom.py'

clra = [0.,0.4,0.8]
clrh = [0.9,0.3,0]

ctl1 = '$T^oC$, '+rdate
ctl2 = 'S '+rdate

# Plotting:
plt.ion()
fig1 = plt.figure(1,figsize=(9,8), constrained_layout=False)
fig1.clf()
#
# T profile
ax1 = plt.axes([0.08, 0.1, 0.25, 0.8]) 
plt.plot(Tar, ZTar, '.-', color=clra, linewidth=2, label="Argo")
plt.plot(Th, ZMh, '.-', color=clrh, linewidth=2, label="RTOFS")
plt.yticks(np.arange(-2000,0,50))
plt.xticks(np.arange(10,30,5))
plt.ylim(-400,0)
plt.xlim(11,29)
ax1.grid(True)
ax1.set_title(ctl1)

# S profile
ax2 = plt.axes([0.4, 0.1, 0.25, 0.8]) 
plt.plot(Sar, ZSar, '.-', color=clra, linewidth=2, label="Argo")
plt.plot(Sh, ZMh, '.-', color=clrh, linewidth=2, label="RTOFS")
plt.yticks(np.arange(-2000,0,50))
plt.xticks(np.arange(35,38,0.5))
plt.ylim(-400,0)
plt.xlim(35,37.5)
ax2.grid(True)
ax2.set_title(ctl2)


# Plot map:
LMSK = HH.copy()
LMSK = np.where(LMSK<0.,0.,1.)
lcmp = mutil.clrmp_lmask(2)

ax3 = plt.axes([0.7, 0.5, 0.25, 0.25])
ax3.pcolormesh(LMSK, shading='flat',\
                cmap=lcmp, \
                vmin=0, \
                vmax=1)
#ax3.contour(HH, [0.0], colors=[(0,0,0)], linewidth=1)
ax3.plot(Iex, Jex, marker='.', markersize=8, color=[1.,0.4,0])






