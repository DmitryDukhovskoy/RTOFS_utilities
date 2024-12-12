"""
  Plot variables averaged over some area 
  for ensemble runs
  Show spectra
  SSH fields

  Use N-day av output 3D fields
  oceanm_XXX.nc

  In standar output fields:
# variables in ocean_daily.nc:
# sos - Sea Surface Salinity
# ssh 
# tos - Sea Surface Temperature
# tob - Sea Water Potential Temperature at Sea Floor
# sob - Sea Water Salinity at Sea Floor


"""
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

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# variables in ocean_daily.nc:
# sos - Sea Surface Salinity
# ssh 
# tos - Sea Surface Temperature
# tob - Sea Water Potential Temperature at Sea Floor
# sob - Sea Water Salinity at Sea Floor
#
varnm  = 'ssh'  # temp (potential) / salin
# expt_nmb - for runs with dailyOB, climatOB - not needed
expt_nmb = 3    # =1 - OBs from fixed SPEAR ens #, =2 - OBs from multi-ens. SPEAR 
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1    
nens   = 1  # ens #
regn   = 'poly_south'  # region to do the averaging over
dnmbS   = mtime.datenum([YRS,MOS,DDS])
dv_start = mtime.datevec(dnmbS)
plt_stdoutp = False # also plot std output 

# Standard output in ocean_monthly.nc, ocean_daily.nc
archv_fl = 'ocean_daily.nc' # output file name 

# For 3D N-daily av. output:
lr     = -1             # vertical layer to analyze, =-1 for 2D output fields
match varnm:
  case "tos":
    lr = 1
    varnm_nc = 'potT'
  case "sos":
    lr = 1
    varnm_nc = 'salt'
  case "tob":
    lr = 100
    varnm_nc = 'potT'
  case "sob":
    lr = 100
    varnm_nc = 'salt'
  case "ssh":
    lr = -1
    varnm_nc = 'ssh'

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
ndav       = pthseas['MOM6_NEP'][expt]['ndav']  # # of days output averaged

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

# Select some region:
# Get indices of the polygon:
II = pthseas['ANLS_NEP'][regn]['II']
JJ = pthseas['ANLS_NEP'][regn]['JJ']
jdm, idm = HH.shape

DX, DY = mmom6.dx_dy(hlon, hlat)
Acell  = DX*DY
X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MS, _, _ = mmisc.inpolygon_v2(X, Y, II, JJ)  # 
JBS, IBS = np.where( (MS == 1) & (HH < 0) ) #exclude deeep regions
MSKBS  = np.zeros((jdm,idm))
MSKBS[JBS,IBS] = 1


match expt:
  case "seasonal_fcst":
    runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
    pthwoutp   = pthseas['MOM6_NEP'][expt]['pthwoutp'].format(runname=runname)
    pthfcst  = pthseas['MOM6_NEP'][expt]['pthwoutp'].format(runname=runname)
  case "seasonal_daily":
    runname  = f'NEPphys_frcst_dailyOB-expt{expt_nmb:02d}'
    pthoutp    = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=expt_nmb)
    pthfcst   = os.path.join(pthoutp,f"{YRS}-{MOS:02d}-e{nens:02d}","history")

print(f"Processing ens={nens:02d} {pthfcst}")
# N-day average output:
Fts, TM = manseas.timeser_spatavrg_dayoutp(pthfcst, YRS, MOS, varnm_nc, lr, MSKBS, Acell)
# Daily output from standard output files:
Fts_std, TM_std = manseas.timeser_spatavrg_stdoutp(pthfcst, archv_fl, varnm, lr, MSKBS, Acell)


DAYS = TM-TM[0]
dT  = DAYS[1] - DAYS[0]

DAYS_std = TM_std-TM_std[0]
dT_std = DAYS_std[1] - DAYS_std[0]

#Qfr = fftpack.fft(SSH)  # Fourier transfer
#Frw = fftpack.fftfreq(len(DAYS), dT)

# Detrend:
Pcoef = np.polyfit(DAYS,Fts,4)
Plnm  = np.poly1d(Pcoef)
Pfit = Plnm(DAYS)

Fts_dtr = Fts-Pfit

N = len(Fts_dtr)
Qfr = sfft.rfft(Fts_dtr) / N
Frw = sfft.rfftfreq(n=N, d=dT/365.) # unit = 1/12 of sampling period
Frw[0] = 1e-20
Fcday = Frw/365.  # cyc / day
Qfr = np.abs(Qfr)
Qfr[0] = np.nan

Pcoef = np.polyfit(DAYS_std,Fts_std,4)
Plnm  = np.poly1d(Pcoef)
Pfit_std = Plnm(DAYS_std)
Fts_std_dtr = Fts_std-Pfit_std

N2 = len(Fts_std_dtr)
Qfr_std = sfft.rfft(Fts_std_dtr) / N2
Frw_std = sfft.rfftfreq(n=N2, d=dT_std/365.) # unit = 1/12 of sampling period
Frw_std[0] = 1.e-20
Fcday_std = Frw_std/365.  # cyc/day
Qfr_std = np.abs(Qfr_std)
Qfr_std[0] = np.nan


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
    yl1, yl2 = 0.2, 0.5

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

ax1.plot(DAYS_std,Fts_std,'-', color=[0.,0.9, 0.2])
ax1.plot(DAYS_std,Pfit_std, '-', color=[1.,0., 0.9])

ax1.set_ylim([yl1,yl2])
ax1.grid('on')
ax1.set_ylabel(f'{varnm}')

#archv_fl = 'oceanm_YYYY_DAY.nc'
dstart = f'{dv_start[0]}/{dv_start[1]}/{dv_start[2]}'
sttl = f'{runname}\n Seas f/cast init: {dstart}, {varnm} {archv_fl} {regn}'
ax1.set_title(sttl)
ax1.set_xlabel('F/cast days')

ax2  = plt.axes([0.1, 0.12, 0.85, 0.34])
ax2.plot(Fcday, np.abs(Qfr))
ax2.plot(Fcday_std,Qfr_std, '-', color=[0.,0.9, 0.2])
#ax2.plot(Frq_avg, Qfr_avg)
ax2.set_yscale('log')
ax2.set_xscale('log')

sttl2 = 'Spectrum, m2/day2'
ticks = ax2.get_xticks()
#ax2.set_xticklabels([f'{tick/N:6.2f}' if tick!=0 else '$\infty$' for tick in ticks])
ax2.set_xlabel('cyc/day')

btx = 'timeser_2Davrg_spectra_dayoutp.py'
bottom_text(btx, pos=[0.04, 0.02])


