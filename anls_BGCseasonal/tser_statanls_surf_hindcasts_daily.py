"""
  Stat analysis of 2 time series from 
  2 hindcasts

  use daily fileds (high auto-correlation)

  from BGC and phys hindcasts
  ice_daily.nc:

  ssh, tos, sos, ssu, ssv

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
import argparse

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
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

# Select location
jj0 = 113
ii0 = 205

parser = argparse.ArgumentParser()
parser.add_argument("--yrs", help="year start to calc RMSE: 1993, ..., 2020", type=int, required=True)
parser.add_argument("--yre", help="year end RMSE", type=int)
parser.add_argument("--varnm", help="field to analyze: ssh, tos, sos, ssu, ssv", type=str, required=True)
parser.add_argument("--ij", help=f"Test i and j indices, default = {ii0}, {jj0}", 
                    nargs=2, type=int)

args = parser.parse_args()

YRS   = args.yrs if args.yrs else None
YRE   = args.yre if args.yre else YRS
varnm = args.varnm if args.varnm else None
if args.ij is not None:
  ii0, jj0 = args.ij
else:
  ii0, jj0 = ii0, jj0

use_era_crpt = True 
hcast_bgc = 'NEPbgc_nudged_hindcast02'
hcast_phys = 'NEPphys_nudged_hindcast' 


hcst_time = 3 # f/csat time interval, months
hcst_interv = np.array([x for x in range(1,12+hcst_time,hcst_time)], dtype=int)


fyaml = 'bgc_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP']['hindcast']['pthgrid']
fgrid      = pthseas['MOM6_NEP']['hindcast']['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"]['hindcast']["ftopo"]
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
jdm, idm = HH.shape
LMsk = np.where(HH<0, 1, 0)

assert HH[jj0,ii0] < 0., f"Test points j={jj0} i={ii0} are on land {HH[jj0,ii0]:.2f}"

MCAL = []
YCAL = []
for yr in range(YRS,YRE+1):
  for mm in range(1,13,3):
    if yr == 1993 and mm == 1:
      continue
    MCAL.append(mm)
    YCAL.append(yr)

len_cal = len(MCAL)  
  
icc = 0
TM  = [] 
tser_crpt = []  # hindcast with corrupted ERA5 BGC
tser_crct = []  # hindcast with corrected ERA5 BGC
tser_phys  = []  # h/cast with physics

for ii in range(len(MCAL)):
  MINIT = MCAL[ii]
  YY = YCAL[ii] 

  pthhcst_bgc = (
       pthseas['MOM6_NEP']['hindcast']['pthoutp'].format(
         hindcast_name=hcast_bgc, YR=YY, MM=MINIT, DD=1
    )
  )
  # Hindcast with corrputed ERA5:
  # Segment 1993-1999 of hindcast02 was replaced with hindcast03 
  # that was run with corrected ERA5
  # the rest of hindcast02 not changed
  if YY <= 1999:
    pthhcst_crpt = '/archive/Dmitry.Dukhovskoy/fre/NEP/hindcast_bgc/NEPbgc_nudged_hindcast02/history/' +\
                  f'{YY}{MINIT:02d}01crpt_forc'

  pthhcst_phys = '/archive/Dmitry.Dukhovskoy/fre/NEP/2024/NEP_physics_202404_nudging-15d/' + \
                 f'gfdl.ncrc5-intel22-repro/history/{YY}-{MINIT:02d}'


  print(f'Extracting {varnm} for {YY}/{MINIT}')
  dnmb_ref = mtime.datenum([1993,1,1])

  # Read hindcast:
  docn_bgc = os.path.join(pthhcst_bgc,f'ocean_daily.nc')
  print(f'Reading {varnm} from {docn_bgc}')
  with xarray.open_dataset(docn_bgc, decode_times=False) as ds_bgc:
    A3d_bgc = ds_bgc[varnm].data.squeeze()
    TIME = ds_bgc['time'].data + dnmb_ref

  docn_crpt = os.path.join(pthhcst_crpt,f'ocean_daily.nc')
  print(f'Reading {varnm} from {docn_crpt}')
  with xarray.open_dataset(docn_crpt, decode_times=False) as ds_bgc:
    A3d_crpt = ds_bgc[varnm].data.squeeze()

  docn_phys = os.path.join(pthhcst_phys,f'ocean_daily.nc')
  print(f'Reading {varnm} from {docn_phys}')
  with xarray.open_dataset(docn_phys) as ds_phys:
    A3d_phys = ds_phys[varnm].data.squeeze()

  icc = -1
  rmse_max = 0.
  for dnmb0 in TIME:
    icc += 1
    A2d_bgc  = A3d_bgc[icc,:,:].squeeze()
    A2d_crpt = A3d_crpt[icc,:,:].squeeze()
    A2d_phys = A3d_phys[icc,:,:].squeeze()
 
    tser_crct.append(A3d_bgc[icc,jj0,ii0])
    tser_crpt.append(A3d_crpt[icc,jj0,ii0])
    tser_phys.append(A3d_phys[icc,jj0,ii0])
    TM.append(dnmb0)  

tser_crct = np.array(tser_crct)
tser_crpt = np.array(tser_crpt)
tser_phys = np.array(tser_phys)

TM = np.array(TM)
DV = mtime.datevec1D(TM, fHR=False)
DV = np.array(DV).transpose()
TMyr = TM*0.
for ik in range(len(TM)):
  dnmb = TM[ik]
  yr, jday = mtime.dnmb2jday(dnmb)
  if yr%4 == 0:
    ndays = 366
  else:
    ndays = 365
  TMyr[ik] = yr + np.floor(jday-1)/ndays


# Remove seasonality:
from scipy.linalg import lstsq
from scipy.signal import correlate

nrec = tser_crct.size
npar = 11
xtime = np.arange(nrec).astype('float')
freq = 2*np.pi/np.max(xtime)

A = np.ones((nrec,npar))
A[:,1] = np.cos(0.5*freq*xtime)
A[:,2] = np.sin(0.5*freq*xtime)
A[:,3] = np.cos(freq*xtime)
A[:,4] = np.sin(freq*xtime)
A[:,5] = np.cos(2*freq*xtime)
A[:,6] = np.sin(2*freq*xtime)
A[:,7] = np.cos(3*freq*xtime)
A[:,8] = np.sin(3*freq*xtime)
A[:,9] = np.cos(4*freq*xtime)
A[:,10] = np.sin(4*freq*xtime)

Xcrct, residuals, rank, s = lstsq(A, tser_crct)
Xcrpt, residuals, rank, s = lstsq(A, tser_crpt)
Xphys, residuals, rank, s = lstsq(A, tser_phys)

#tser_crct1D = tser_crct[:, np.newaxis]
anom_crct = tser_crct - (A @ Xcrct)
anom_crpt = tser_crpt - (A @ Xcrpt)
anom_phys = tser_phys - (A @ Xphys)


def xcorr(x,y):
  xcorr = correlate(x - np.mean(x), y - np.mean(y), mode='full')
  lags = np.arange(-len(x)+1, len(x))
  xcorr_norm = xcorr / (np.std(x) * np.std(y) * len(x))

  return xcorr_norm, lags

xr_crct_crpt, lags = xcorr(anom_crct, anom_crpt)
xr_crpt_phys, lags = xcorr(anom_crpt, anom_phys)


# Differencing:
dfcrct = np.diff(anom_crct)
dfcrpt = np.diff(anom_crpt)
dfphys = np.diff(anom_phys)

xr_dcrct_dphys, lags = xcorr(dfcrct, dfphys) 




plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
# Time series RMSE
ax1 = plt.axes([0.1, 0.5, 0.85, 0.4])
ax1.plot(TMyr,RMSE, '-', linewidth=2, color=[0.0,0.4,0.9])
sttl = f'RMSE {varnm} NEPbgc vs NEPphys hindcasts {YRS}-{YRE}'
if f_crct:
  ax1.plot(TMyr_crct,RMSE_crct, '-', linewidth=2, color=[0.0,0.9,0.4])
  sttl = f'RMSE {varnm} NEPbgc02 /NEPbgc03 (green) vs NEPphys hindcasts {YRS}-{YRE}'

ax1.set_title(sttl)
ax1.grid('on')

btx = 'calc_RMSE_surf_hindcasts_daily.py'
bottom_text(btx, pos=[0.05,0.35])

