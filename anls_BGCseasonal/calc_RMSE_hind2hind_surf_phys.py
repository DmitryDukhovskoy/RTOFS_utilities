"""
  Calc and plot RMSE for daily surface fields (ocean_daily.nc)
  from BGC hindcasts 02 (corrupted ERA5) and 03 (corrected ERA5)
  for physical variables

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
parser.add_argument("--varnm", help="field to analyze: ssh, tos, sos, ssu, ssv", type=str, required=True)
args = parser.parse_args()

varnm = args.varnm if args.varnm else None

hcast2_bgc = 'NEPbgc_nudged_hindcast02'
hcast3_bgc = 'NEPbgc_nudged_hindcast03'

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


YRS = 1993
YRE = 1996
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

  pthhcst2 = (
       pthseas['MOM6_NEP']['hindcast']['pthoutp'].format(
         hindcast_name=hcast2_bgc, YR=YY, MM=MINIT, DD=1
    )
  )

  pthhcst3 = '/archive/Dmitry.Dukhovskoy/fre/NEP/hindcast_bgc/NEPbgc_nudged_hindcast03/history/' +\
            f'{YY}{MINIT:02d}01/'


  print(f'Calculating RMSE for {YY}/{MINIT}')
  dnmb_ref = mtime.datenum([1993,1,1])

  # Read hindcasts:
  docn2 = os.path.join(pthhcst2,f'ocean_daily.nc')
  print(f'Reading {varnm} from {docn2}')
  with xarray.open_dataset(docn2, decode_times=False) as ds_bgc2:
    A3d_bgc2 = ds_bgc2[varnm].data.squeeze()
    TIME2 = ds_bgc2['time'].data + dnmb_ref

  docn3 = os.path.join(pthhcst3,f'ocean_daily.nc')
  print(f'Reading {varnm} from {docn3}')
  with xarray.open_dataset(docn3, decode_times=False) as ds_bgc3:
    A3d_bgc3 = ds_bgc3[varnm].data.squeeze()
    TIME3 = ds_bgc3['time'].data + dnmb_ref

  # Time series should be over the same time intervals:
  D = np.abs(TIME3-TIME2)
  assert np.max(D) < 1.e-19, f'ERROR: Time arrays in hindcast2 and hindcast3 are not the same' 
  TIME = TIME2.copy()
   
  icc = -1
  rmse_max = 0.
  for dnmb0 in TIME:
    icc += 1
    A2d_bgc2 = A3d_bgc2[icc,:,:].squeeze()
    A2d_bgc3 = A3d_bgc3[icc,:,:].squeeze()
 
    sqerr = (A2d_bgc2-A2d_bgc3)**2
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


plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
# Time series RMSE
ax1 = plt.axes([0.1, 0.5, 0.85, 0.4])
ax1.plot(TMyr,RMSE, '-', linewidth=2, color=[0.0,0.4,0.9])

sttl = f'RMSE {varnm} NEPbgc vs NEPphys hindcasts {YRS}-{YRE}'
ax1.set_title(sttl)
ax1.grid('on')

btx = 'calc_RMSE_surf_hindcasts_daily.py'
bottom_text(btx, pos=[0.05,0.35])

