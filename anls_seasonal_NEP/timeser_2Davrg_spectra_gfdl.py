"""
  Plot SSH  averaged over some area and its power spectral density
  for ensemble runs
  GFDL simulations (Liz's paper, original hindcast

  Open lateral boundary and initial conditions for temperature, salinity, 
  sea surface height and momentum were prescribed as daily means from the 
  1/12° Global Ocean Physics Reanalysis 

  Amt: ERA-5

  Use standard output daily files
  post-processed into separate files
  e.g.: 
  ocean_daily.19930101-19971231.ssh.nc 

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

expt_nmb = 1    # =1 - OBs from fixed SPEAR ens #, =2 - OBs from multi-ens. SPEAR 
YRS    = 1993 # year start of the forecast
nyrav  = 5
MOS    = 4
DDS    = 1    
nensR  = 1  #ens # for reference ensemble run
regn   = 'poly_north'  # region to do the averaging over
outp_span = f'{YRS}0101-{YRS+nyrav-1}1231'
archv_fl = f'ocean_daily.{outp_span}.{varnm}.nc' # output file name 

run_nm = 'MOM6_NEP_GFDL'
expt   = 'NEP_BGC_seas'
runname= 'GFDL NEP hindcast 1993-2017'
dnmbS    = mtime.datenum([YRS,MOS,DDS]) 
dv_start = mtime.datevec(dnmbS)

print(f'Expt: {expt} Run: {runname} init date: {YRS}/{MOS}/{DDS}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

ptharch    = '/archive/e1n/fre/cefi/NEP/2024_07/NEP_cefi_bgc_072024/gfdl.ncrc5-intel22-prod/'+\
             'pp/ocean_daily/ts/daily/5yr/'

#pthoutp    = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=expt_nmb)
#pthwoutp   = pthseas['MOM6_NEP'][expt]['pthwoutp'].format(runname=runname)
#pthfcst    = pthseas[run_nm][expt]['pthoutp']
pthtopo    = pthseas[run_nm][expt]['pthgrid']
fgrid      = pthseas[run_nm][expt]['fgrid']
ftopo_mom  = pthseas[run_nm][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

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

#ENSR = [x for x in range(1,11)]
#NensR = len(ENSR)
pthfcst = ptharch

print(f"Processing {pthfcst}")
lr = -1
nens = 1
Fts, TM = manseas.timeser_spatavrg_stdoutp(pthfcst, archv_fl, varnm, lr, MSKBS, Acell)

#if iens == 0:
dim1 = "Time"
darr_var = xarray.DataArray(Fts, dims=(dim1), coords={dim1: TM})
dset1D = xarray.Dataset({f"{varnm}_e{nens:02d}": darr_var}).copy()
#else:
#  darr_var = xarray.DataArray(Fts, dims=(dim1), coords={dim1: TM})
#  dset_var = xarray.Dataset({f"{varnm}_e{nens:02d}": darr_var})
#  dset1D = xarray.merge([dset1D, dset_var])

# 1 year of data:
ndays = 365
Yts = Fts[:ndays]
YTM = TM[:ndays]
DAYS = TM-TM[0]
dT  = DAYS[1] - DAYS[0]
YDAYS = YTM-YTM[0]

#Qfr = fftpack.fft(SSH)  # Fourier transfer
#Frw = fftpack.fftfreq(len(DAYS), dT)

# Detrend:
Pcoef = np.polyfit(YDAYS,Yts,4)
Plnm  = np.poly1d(Pcoef)
Pfit = Plnm(YDAYS)

Yts_dtr = Yts-Pfit
N = len(Yts_dtr)
Qfr = sfft.rfft(Yts_dtr) / N
Frw = sfft.rfftfreq(n=N, d=dT/ndays) # unit = 1/12 of sampling period
Frw[0] = 1e-20
Fcday = Frw/ndays  # cyc / day
Qfr = np.abs(Qfr)
Qfr[0] = np.nan


# 
match regn:
  case "poly_south":
    if varnm == 'tos':
      yl1, yl2 = 17.5, 27.0
    elif varnm == 'sos':
      yl1, yl2 = 33.95, 34.3
    elif varnm == 'ssh':
      yl1, yl2 = 0.15, 0.4
    elif varnm == 'tob':
      yl1, yl2 = 3.0, 4.5
    elif varnm == 'sob':
      yl1, yl2 = 34.59, 34.63

  case "poly_central":
    if varnm == 'tos':
      yl1, yl2 = 17.5, 27.0
    elif varnm == 'sos':
      yl1, yl2 = 32.0, 32.8
    elif varnm == 'ssh':
      yl1, yl2 = 0.05, 0.25
    elif varnm == 'tob':
      yl1, yl2 = 1.9, 2.3
    elif varnm == 'sob':
      yl1, yl2 = 34.48, 34.6

  case "poly_north":
    if varnm == 'tos':
      yl1, yl2 = 10.5, 27.0
    elif varnm == 'sos':
      yl1, yl2 = 32.0, 32.8
    elif varnm == 'ssh':
      yl1, yl2 = -0.15, 0.15
    elif varnm == 'tob':
      yl1, yl2 = 1.9, 2.3
    elif varnm == 'sob':
      yl1, yl2 = 34.48, 34.6


# ===================
# Plotting
# ===================
#btx = 'timeser_2Davrg_ensmbl.py' 
btx = 'timeser_2Davrg_spectra_gfdl.py'

plt.ion()

fgnmb = 1
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()

ax1  = plt.axes([0.1, 0.54, 0.85, 0.4])
LNS  = []

vards = f"{varnm}_e{nens:02d}"
Ts = dset1D[vards].data.squeeze()
Ts = Ts[:ndays]
TM = dset1D['Time'].data[:ndays]
Tday = TM-TM[0]

ax1.plot(Tday, Ts, '-', label=f'ens{nens:02d}')
ax1.plot(YDAYS,Pfit)


ax1.set_ylim([yl1,yl2])
ax1.grid('on')
ax1.set_xlabel('Forecast days')
ax1.set_ylabel(f'{varnm}')
ax1.set_xlim([0, ndays])

dstart = f'{dv_start[0]}/{dv_start[1]}/{dv_start[2]}'
sttl = f'{runname}\n GFDL NEP hindcast, 1993 {varnm} {archv_fl} {regn}'
ax1.set_title(sttl)

ax2  = plt.axes([0.1, 0.12, 0.85, 0.34])
ax2.plot(Fcday, np.abs(Qfr))
#ax2.plot(Frq_avg, Qfr_avg)
ax2.set_yscale('log')
ax2.set_xscale('log')

sttl2 = 'Spectrum, m2/day2'
ticks = ax2.get_xticks()
#ax2.set_xticklabels([f'{tick/N:6.2f}' if tick!=0 else '$\infty$' for tick in ticks])
ax2.set_xlabel('cyc/day')

bottom_text(btx)


