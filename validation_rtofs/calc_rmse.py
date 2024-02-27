# Calc RMSE for specified regions
#
# Validation: RTOFS ssh vs altimetry-derived SSH global gridded 
# CMEMS near-real time copernicus
# https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046/description
#
# This product is processed by the DUACS multimission altimeter 
# data processing system.
# DUACS  is the operationnal multimission production system of 
# altimeter data developed by CNES/CLS.
#
# downloaded from
# https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046/download?dataset=dataset-duacs-nrt-global-merged-allsat-phy-l4
# 

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
import time
from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text

import mod_utils as mutil
import mod_misc1 as mmisc
#importlib.reload(mmisc)
import mod_time as mtime

#expt   = 'paraD5'  # paraXX - parallel runs or product - production
expt   = 'product'
rdateS  = '20230303'
rdateE  = '20230508'
#rdate0 = '20230314'
sfx     = 'n-24'
regn_nm = 'Agulhas' # GOM, Carib, GulfStr, SOcean, Kurosh, Agulhas

f_save      = True
f_save_mean = True   # saved global time-mean RMSE - need for diferent regions
                     # because of different spatial SSH means used for 
                     # demeaning regional SSH
pthsave = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/data_outp/ssh_intrp2cmems/'

#
# Create land-ocean mask from SSH cmems
pthssh = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/OBS/ssh_altimetry/'
flnm   = 'ssh_duacs_nrt_global_20230301_20230430.nc'
flssh  = pthssh + flnm

nc    = ncFile(flssh,'r')
XA    = nc.variables['longitude'][:].squeeze().data
YA    = nc.variables['latitude'][:].squeeze().data
TMA   = nc.variables['time'][:].squeeze().data
idma  = XA.shape[0]
jdma  = YA.shape[0]
ijdma = idma*jdma
SSHA  = nc['adt'][0,:,:]

Lmsk = np.where(np.isnan(SSHA),0,1)
Iall = np.where(Lmsk.flatten()>0)[0]
nall = Iall.shape[0]

# Get time array of near-real time ssh altimetry
dnmbR = mtime.rdate2datenum('19500101')  # ref day 
TMA   = TMA + dnmbR
TMA   = TMA.astype(np.int32)

# Define region:
import mod_valid_utils as mvutil
importlib.reload(mvutil)
REGNS = mvutil.ssh_regions()
REGNMS = list(REGNS.keys())
xl1    = REGNS[regn_nm]["xl1"]
xl2    = REGNS[regn_nm]["xl2"]
yl1    = REGNS[regn_nm]["yl1"]
yl2    = REGNS[regn_nm]["yl2"]
Ip     = np.array(REGNS[regn_nm]["Ip"])
Jp     = np.array(REGNS[regn_nm]["Jp"])
# Points inside the region
X, Y = np.meshgrid(np.arange(idma), np.arange(jdma))
Rmsk, IRg, JRg = mmisc.inpolygon_v2(X,Y,Ip,Jp)  # region

# Interpolated SSH dates:
dnmbS = int(mtime.rdate2datenum(rdateS))
dnmbE = int(mtime.rdate2datenum(rdateE))
TM    = np.arange(int(dnmbS),int(dnmbE)+1)

RMSE_TS = TM*0.

icnt  = 0
icc   = 0
for dnmb0 in range(dnmbS,dnmbE+1):
  rdate0 = mtime.dnumb2rdate(dnmb0, ihours=False)
  print('Processing ' + rdate0)

  fsshout = 'sshintrp_' + expt + '_' + rdate0 + sfx + '.pkl'
  flintrp  = pthsave + fsshout
  if not os.path.isfile(flintrp):
    print('File not found ' + flintrp)
    print('skipping ...')
    continue

  print('Reading  ' + flintrp)
  with open(flintrp,'rb') as fid:
    SSHi = pickle.load(fid)

# Find time in altimetry data
# for n-24 - previous date:
  if sfx == 'n-24':
    dd0 = dnmb0-1
  else:
    dd0 = dnmb0

  SSHA = mvutil.read_ssh_cmems(dd0)
#  rdate_cmems = mtime.dnumb2rdate(dd0,ihours=False)
#  iTime = np.where(TMA == dd0)[0][0]
#  dnmbA = TMA[iTime]
#  SSHA = nc['adt'][iTime,:,:]
 
# Demean ssh using regional mean:
  sshi_sub = SSHi[JRg,IRg]
  sshi_mn  = np.nanmean(sshi_sub)
  dSSHi    = SSHi - sshi_mn
 
  ssha_sub = SSHA[JRg,IRg]
  ssha_mn  = np.nanmean(ssha_sub)
  dSSHa    = SSHA - ssha_mn

# Compute squared difference (!!! this is not RMSE yet - need to sq.root)
  RMSE2        = (dSSHi-dSSHa)**2
  npnts        = IRg.shape[0]
  rmse_mn      = np.sqrt(np.sum(RMSE2[JRg,IRg])/npnts) 
  RMSE_TS[icc] = rmse_mn

  icc += 1
  if icc==1:
    RMSE_mn = RMSE2
  else:
    RMSE_mn = RMSE_mn + RMSE2


RMSE_mn = np.sqrt(RMSE_mn/float(icc))

if f_save:
#  frmseout = 'ssh_rmse2_' + expt + sfx + '_' + regn_nm + '_' + rdate0 +  '.pkl'
  frmseout = 'ssh_rmse_' + expt + sfx + '_' + regn_nm + '.pkl'
  flrmse  = pthsave + frmseout
  print('Saving  ' + flrmse + ' \n')
  with open(flrmse,'wb') as fid:
    pickle.dump([RMSE_TS,TM],fid)
   
# Saves global spatial mean RMSE 2D field
if f_save_mean:
  frmsemn = 'ssh_rmsemn2D_' + expt + sfx + '_'  + regn_nm + '.pkl'
  flrmsemn  = pthsave + frmsemn
  print('Saving  ' + flrmsemn + ' \n')
  with open(flrmsemn,'wb') as fid:
    pickle.dump([RMSE_mn,Lmsk],fid)

 
print('All done')


# Plot interpolated SSH
# Function to print mouse click event coordinates
def onclick(event):
   print([event.xdata, event.ydata])

f_plotssh = False
if f_plotssh:
  plt.ion()

# Plot region and select points if needed:
  fig1 = plt.figure(1,figsize=(8,7), constrained_layout=False)
  fig1.clf()

  ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8]) 
  im1 = ax1.contour(Lmsk,[0.1],colors=[(0,0,0)])
  f_setrmu = True
  if f_setrmu:
  # Bind the button_press_event with the onclick() method
    fig1.canvas.mpl_connect('button_press_event', onclick)

  ax1.axis('scaled')
  ax1.set_xlim([317,456])
  ax1.set_ylim([384,492])

#
# Plot demeaned SSH in the region
  fig2 = plt.figure(2,figsize=(8,7), constrained_layout=False)
  fig2.clf()

  ax21 = fig2.add_axes([0.1, 0.55, 0.8, 0.4])
  im21 = ax21.pcolormesh(dSSHi, vmin=-0.5, vmax=0.5)
  ax21.axis('scaled')
  ax21.set_xlim([xl1,xl2])
  ax21.set_ylim([yl1,yl2])
  ax21.set_title(expt + ' ssh-<ssh> ' + rdate0 + sfx + \
                ' <ssh>={0:5.3f}m'.format(sshi_mn))

  cax = fig2.add_axes([0.9, 0.55, 0.015, 0.4])
  clb = fig2.colorbar(im21, cax=cax, extend='both')
  cax.set_yticklabels(cax.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=5)

  ax22 = fig2.add_axes([0.1, 0.07, 0.8, 0.4])
  im22 = ax22.pcolormesh(dSSHA, vmin=-0.5, vmax=0.5)
  ax22.axis('scaled')
  ax22.set_xlim([xl1,xl2])
  ax22.set_ylim([yl1,yl2])
  ax22.set_title('CMEMS ssh-<ssh> ' + rdate_cmems + \
              ' <ssh>={0:5.3f}m'.format(ssha_mn))

  cax2 = fig2.add_axes([0.9, 0.07, 0.015, 0.4])
  clb2 = fig2.colorbar(im22, cax=cax2, extend='both')
  cax2.set_yticklabels(cax2.get_yticks())
  ticklabs = clb2.ax.get_yticklabels()
#  clb2.ax.set_yticklabels(ticklabs,fontsize=10)
  clb2.ax.set_yticklabels(["{:.2f}".format(i) \
                           for i in clb2.get_ticks()], fontsize=10)
  clb2.ax.tick_params(direction='in', length=5)

  btx = 'calc_rmse.py'
  bottom_text(btx, pos=[0.02,0.01])





