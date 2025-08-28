"""
  Calc and plot RMSE for daily fields
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


parser = argparse.ArgumentParser()
parser.add_argument("--yrs", help="year start to calc RMSE: 1993, ..., 2020", type=int, required=True)
parser.add_argument("--yre", help="year end RMSE", type=int)
parser.add_argument("--varnm", help="field to analyze: ssh, tos, sos, ssu, ssv", type=str, required=True)
args = parser.parse_args()

YRS   = args.yrs if args.yrs else None
YRE   = args.yre if args.yre else YRS
varnm = args.varnm if args.varnm else None

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
RMSE = []
for ii in range(len(MCAL)):
  MINIT = MCAL[ii]
  YY = YCAL[ii] 

  # Find init date for given month, assuming hcst_time (n months) f/cast interval
  #kint = np.searchsorted(hcst_interv, MM, side='right') - 1
  #assert(hcst_interv[kint] <= MM < hcst_interv[kint+1]), f'Wrong time bin {kint} for {MMA}'
  #MINIT = hcst_interv[kint]
  #imo = MM-MINIT      # current month in the archive output

  pthhcst_bgc = (
       pthseas['MOM6_NEP']['hindcast']['pthoutp'].format(
         hindcast_name=hcast_bgc, YR=YY, MM=MINIT, DD=1
    )
  )

  pthhcst_phys = '/archive/Dmitry.Dukhovskoy/fre/NEP/2024/NEP_physics_202404_nudging-15d/' + \
                 f'gfdl.ncrc5-intel22-repro/history/{YY}-{MINIT:02d}'


  print(f'Calculating RMSE for {YY}/{MINIT}')
  dnmb_ref = mtime.datenum([1993,1,1])

  # Read hindcast:
  docn_bgc = os.path.join(pthhcst_bgc,f'ocean_daily.nc')
  print(f'Reading {varnm} from {docn_bgc}')
  with xarray.open_dataset(docn_bgc, decode_times=False) as ds_bgc:
    A3d_bgc = ds_bgc[varnm].data.squeeze()
    TIME = ds_bgc['time'].data + dnmb_ref

  docn_phys = os.path.join(pthhcst_phys,f'ocean_daily.nc')
  print(f'Reading {varnm} from {docn_phys}')
  with xarray.open_dataset(docn_phys) as ds_phys:
    A3d_phys = ds_phys[varnm].data.squeeze()

  icc = -1
  rmse_max = 0.
  for dnmb0 in TIME:
    icc += 1
    A2d_bgc  = A3d_bgc[icc,:,:].squeeze()
    A2d_phys = A3d_phys[icc,:,:].squeeze()
 
    sqerr = (A2d_bgc-A2d_phys)**2
    nB = np.count_nonzero(~np.isnan(sqerr))
    rmseBP = np.sqrt(np.nansum(sqerr)/nB)

    RMSE.append(rmseBP)
    TM.append(dnmb0)  
    rmse_max = np.max([rmse_max,rmseBP])

  print(f"max RMSE: {rmse_max:.4f}")


RMSE = np.array(RMSE)

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


f_crct = True
if f_crct:
  # Check corrected hindcasts with corrected ERA5
  # to replace impacted hindcasts 1993-1995 that were run
  # with corrupted ERA5 forcing
  print('Computing RMSE for corrected hindcast 03')
  RMSE_crct = []
  TM_crct = []
  for YY in range(1993,1997): 
    for MINIT in range(1,13,3):
      if YY == 1993 and MINIT == 1:
        continue
      # hindcast with corrected ERA5
      pthcrct = '/archive/Dmitry.Dukhovskoy/fre/NEP/hindcast_bgc/NEPbgc_nudged_hindcast03/history/'+\
                f'{YY}{MINIT:02d}01'
      docn_bgc = os.path.join(pthcrct,f'ocean_daily.nc')
      print(f'Reading {varnm} from {docn_bgc}')
      with xarray.open_dataset(docn_bgc, decode_times=False) as ds_bgc:
        A3d_bgc = ds_bgc[varnm].data.squeeze()
        TIME = ds_bgc['time'].data + dnmb_ref

      pthhcst_phys = '/archive/Dmitry.Dukhovskoy/fre/NEP/2024/NEP_physics_202404_nudging-15d/' + \
                     f'gfdl.ncrc5-intel22-repro/history/{YY}-{MINIT:02d}'
      docn_phys = os.path.join(pthhcst_phys,f'ocean_daily.nc')
      print(f'Reading {varnm} from {docn_phys}')
      with xarray.open_dataset(docn_phys) as ds_phys:
        A3d_phys = ds_phys[varnm].data.squeeze()

      icc = -1
      rmse_max = 0.
      for dnmb0 in TIME:
        icc += 1
        A2d_bgc  = A3d_bgc[icc,:,:].squeeze()
        A2d_phys = A3d_phys[icc,:,:].squeeze()

        sqerr = (A2d_bgc-A2d_phys)**2
        nB = np.count_nonzero(~np.isnan(sqerr))
        rmseBP = np.sqrt(np.nansum(sqerr)/nB)

        RMSE_crct.append(rmseBP)
        TM_crct.append(dnmb0)
        rmse_max = np.max([rmse_max,rmseBP])

  RMSE_crct = np.array(RMSE_crct)
  TM_crct = np.array(TM_crct) 

  TMyr_crct = TM_crct*0.
  for ik in range(len(TM_crct)):
    dnmb = TM_crct[ik]
    yr, jday = mtime.dnmb2jday(dnmb)
    if yr%4 == 0:
      ndays = 366
    else:
      ndays = 365
    TMyr_crct[ik] = yr + np.floor(jday-1)/ndays


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

