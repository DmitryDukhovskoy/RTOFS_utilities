"""
  Plot variables averaged over some area 
  for ensemble runs
  Show spectra
  SSH fields

  SPEAR seasonal f/casts
  for NEP domain

"""

NOT FINISHED - 
see check_ssh_SPEAR.py 
for OBs spectra  

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
from copy import copy
import matplotlib.colors as colors
from yaml import safe_load
#import scipy.fftpack as fftpack
import scipy.fft as sfft

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

from mod_utils_fig import bottom_text
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
importlib.reload(mutob)

#
varnm  = 'zos'  # temp (potential) / salin
# expt_nmb - for runs with dailyOB, climatOB - not needed
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1    
nens   = 1  # ens #
regn   = 'poly_south'  # region to do the averaging over
dnmbS   = mtime.datenum([YRS,MOS,DDS])
dv_start = mtime.datevec(dnmbS)

# Get grid for NEP domain
# seasonal_daily :  seas f/casts with dailyOB from SPEAR
# seasonal_fcst : seas f/casts with climatological OB from SPEAR 
expt     = "seasonal_daily"  # 
dnmbS    = mtime.datenum([YRS,MOS,DDS]) 
dv_start = mtime.datevec(dnmbS)

match expt:
  case "seasonal_fcst":
    expt_ob = "climatOB"
  case "seasonal_daily":
    expt_ob = "dailyOB"

print(f'Expt: {expt} {expt_ob} init date: {YRS}/{MOS}/{DDS}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

# Select some region:
# Get indices of the polygon:
II = pthseas['ANLS_NEP'][regn]['II']
JJ = pthseas['ANLS_NEP'][regn]['JJ']
Xmom = hlon[JJ,II]
Ymom = hlat[JJ,II]


pthglorys = '/archive/e1n/datasets/GLORYS/'
pthdata   = os.path.join(pthglorys, f'{YRS}/nep_10')
foutp     = f'GLORYS_REANALYSIS_NEP_{YRS}-{MOS:02d}-{DDS:02d}.nc' 
ds_glorys = xarray.open_dataset(os.path.join(pthdata,foutp)) 

lon = ds_glorys['longitude'].data
lat = ds_glorys['latitude'].data
ssh  = ds_glorys['zos'].data.squeeze()
Glon, Glat = np.meshgrid(lon ,lat)
jdm, idm = Glon.shape[:2]

IIG = []
JJG = []
for ii in range(len(Xmom)):
  x0, y0 = Xmom[ii], Ymom[ii]
  i0, j0 = mutil.find_indx_lonlat(x0, y0, Glon, Glat)
  IIG.append(i0)
  JJG.append(j0)


MSKBS  = np.where(np.isnan(ssh),0,1)
DX, DY = mmom6.dx_dy(Glon, Glat)
Acell  = DX*DY
X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MS, _, _ = mmisc.inpolygon_v2(X, Y, IIG, JJG)  # 
JBS, IBS = np.where( (MS == 1) & (MSKBS == 1) ) 

print(f"Processing GLORYS start: {YRS}/{MOS}")
lr = -1
ndays = 365
Fts, TM = manseas.timeser_spatavrg_GLORYS(pthglorys, YRS, MOS, varnm, lr, MSKBS, Acell, ndays=ndays)


DAYS = TM-TM[0]
dT  = DAYS[1] - DAYS[0]

#Qfr = fftpack.fft(SSH)  # Fourier transfer
#Frw = fftpack.fftfreq(len(DAYS), dT)

# Detrend:
Pcoef = np.polyfit(DAYS,Fts,4)
Plnm  = np.poly1d(Pcoef)
Pfit = Plnm(DAYS)

Fts_dtr = Fts-Pfit
N = len(Fts_dtr)
Qfr = sfft.rfft(Fts_dtr) / N
Frw = sfft.rfftfreq(n=N, d=dT/ndays) # unit = 1/12 of sampling period
Frw[0] = 1e-20
Fcday = Frw/ndays  # cyc / day
Qfr = np.abs(Qfr)
Qfr[0] = np.nan

# Bin averaging
# Skip the 1st bin which is nan - low-freq removed
Nav = 2
Iav = [x for x in range(1,len(Frw),Nav)]
Iav.append(N)
Frq_avg = []
Qfr_avg = []
for ii in range(len(Iav)-1):
  i1 = Iav[ii]
  i2 = Iav[ii+1]-1
  if i2 < i1: i2=i1
  Frq_avg.append(np.mean(Frw[i1:i2+1]))
  Qfr_avg.append(np.mean(Qfr[i1:i2+1]))


match regn:
  case "poly_south":
    yl1, yl2 = 0.22, 0.4

  case "poly_central":
    yl1, yl2 = 0., 0.35

# ===================
# Plotting
# ===================

plt.ion()

fgnmb = 1
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()
ax1  = plt.axes([0.1, 0.54, 0.85, 0.4])
ax1.plot(DAYS,Fts,'-')
ax1.plot(DAYS,Pfit)

ax1.set_ylim([yl1,yl2])
ax1.grid('on')
ax1.set_ylabel(f'{varnm}')

#archv_fl = 'oceanm_YYYY_DAY.nc'
dstart = f'{dv_start[0]}/{dv_start[1]}/{dv_start[2]}'
sttl = f'GLORYS_nep Reanalysis\n Seas f/cast init: {dstart}, {varnm} {regn}'
ax1.set_title(sttl)
ax1.set_xlabel('F/cast days')

ax2  = plt.axes([0.1, 0.12, 0.85, 0.34])
ax2.plot(Fcday, np.abs(Qfr))
#ax2.plot(Frq_avg, Qfr_avg)
ax2.set_yscale('log')
ax2.set_xscale('log')

sttl2 = 'Spectrum, m2/day2'
ticks = ax2.get_xticks()
#ax2.set_xticklabels([f'{tick/N:6.2f}' if tick!=0 else '$\infty$' for tick in ticks])
ax2.set_xlabel('cyc/day')

btx = 'timeser_2Davrg_spectra_GLORYS.py'
bottom_text(btx, pos=[0.04, 0.02])


