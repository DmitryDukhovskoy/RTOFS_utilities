"""
  NEPphys hindcast GLORYS nudging (no ice relaxation)
  Calc and plot time series of RMSE for monthly sea ice conc/thickness 
  from test simulations against target relaxation fields from PIOMAS

  Statistics are computed for 2 regions: Ber Sea and Arctic Region (Chukchi)

  from test experiments with relaxation and / or target fields (monthly mean)
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
parser.add_argument("--varnm", help="field to use: ithkn or iconc", type=str)
args = parser.parse_args()

varnm = args.varnm if args.varnm else 'iconc'
YRS   = args.yrs if args.yrs else None
YRE   = args.yre if args.yre else YRS
hcast_name = f'NEPphys_nudged_hindcast'

pthrlx = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax' 

outfld  = 'icem'
prfx = ''  # 19930401 - time stamp used in SIS2 output in file names, note that find
           # closest archive file does not work for 19930401.icem*.nc file names
           # rename files using ./rename_archive_v0.sh 0 in the output dir
hcst_time = 3 # f/csat time interval, months
hcst_interv = np.array([x for x in range(1,12+hcst_time,hcst_time)], dtype=int)

if varnm == 'iconc':
  ifld = 'siconc'
elif varnm == 'ithkn':
  ifld = 'sithick'

def read_test_run_monthly(dnmb0, outfld, pthtest, prfx, varnm):
  """
    Some test have monthly fields saved, if not - derive from N-day average
  """
  YYR, MMR = mtime.datevec(dnmb0)[:2]

  # Check if monthly file exists:
  dfmnth = os.path.join(pthtest,'ice_month.nc')
  if os.path.isfile(dfmnth):
    import pandas as pd
    # Read monthly mean from saved ice_month.nc if exists
    dset = xarray.open_dataset(dfmnth)
    tm_nep = dset['time'].data
    tmP = pd.to_datetime(tm_nep)
    nrec = len(tmP)
    TNEP = np.zeros((nrec,3), dtype=int)
    for irec in range(nrec):
      yr0 = int(tmP.year[irec])
      mo0 = int(tmP.month[irec])
      dd0 = int(tmP.day[irec])
      TNEP[irec,:] = [yr0,mo0,dd0]

    yrfcst = TNEP[:,0]
    mmfcst = TNEP[:,1]
    if YYR < np.min(yrfcst) or YYR > np.max(yrfcst):
      raise Exception(f"{YYR} is out of range for saved ice_month: {np.min(yrfcst)}/{np.max(yrfcst)}")
    ifcst  = np.where((mmfcst==MMR) & (yrfcst==YYR))[0][0]

    HIce = dset['sithick'].isel(time=ifcst).data
    CIce = dset['siconc'].isel(time=ifcst).data
    if varnm == 'iconc':
      A2d = CIce
    elif varnm == 'ithkn':
      A2d = CIce*HIce
  else:
    print(f'ice_mmonth.nc not found in {pthtest}, deriving mean ...')

    A2d = msisrlx.calc_iconc_ithkn_mnthmean(dnmb0, pthtest, varnm)

  return A2d

def read_relax_piomas(dnmb0, pthsis, varnm):
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  flthck  = f'piomas_heff{YR0}_v21.nc'
  varthck = 'heff'
  flconc  = f'piomas_area{YR0}_v21.nc'
  varconc = 'area'

  # Read saved relax. fields:
  YR1 = YR0
  YR2 = YR0+1
  flout = f'PIOMASv21_ithkn_iconc_{YR1}_{YR2}_monthly.nc'
  diclim = os.path.join(pthsis, flout)
  print(f'Reading relax fields from {diclim}')
  ds_rlx = xarray.open_dataset(diclim)
  Time = ds_rlx['time'].data
  TM = mmisc.convert_nptime_to_datenum(Time)
  dnmb0 = mtime.datenum([YR0,MM0,15,12])
  D = abs(TM-dnmb0)
  itime = np.argmin(D)
  dv0 = mtime.datevec(TM[itime])
  assert dv0[0]==YR0, f'Requested YR={YR0}, year in rlx file={dv0[0]}'
  assert dv0[1]==MM0, f'Requested month={MM0}, month in rlx file={dv0[1]}'

  match varnm:
    case('ithkn'):
      ifld = 'ithkn'
    case('iconc'):
      ifld = 'iarea'

  A2dS = ds_rlx[ifld].isel(time=itime).data

  return A2dS


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

expt = 'seasonal_daily'
pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
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

# Bering Sea and Arctic region masks:
BMsk = manseas.mask_BeringSea_NEPgrid(HH,hlon,hlat)
AMsk = manseas.mask_Arctic_NEPgrid(HH,hlon,hlat)
JB,IB = np.where(BMsk==1)
JA,IA = np.where(AMsk==1)


MCAL = []
YCAL = []
for yr in range(YRS,YRE+1):
  for mm in range(1,13):
    MCAL.append(mm)
    YCAL.append(yr)

len_cal = len(MCAL)  
  
icc = 0
TM  = [] 
RMSE_Ber = []
RMSE_Arc = []
BIAS_Ber = []
BIAS_Arc = []
for ii in range(len(MCAL)):
  MM = MCAL[ii]
  YY = YCAL[ii] 

  # Find init date for given month, assuming hcst_time (n months) f/cast interval
  kint = np.searchsorted(hcst_interv, MM, side='right') - 1
  assert(hcst_interv[kint] <= MM < hcst_interv[kint+1]), f'Wrong time bin {kint} for {MMA}'
  MINIT = hcst_interv[kint]
  imo = MM-MINIT      # current month in the archive output

  pthroot = '/archive/Dmitry.Dukhovskoy/fre/NEP/2024/NEP_physics_202404_nudging-15d/'+\
            'gfdl.ncrc5-intel22-repro/history'
  pthhcst = os.path.join(pthroot,f'{YY}-{MINIT:02d}')

  print(f'Calculating RMSE for {YY}/{MM} {pthhcst}')
  dnmbR = mtime.datenum([YY,MM,15])

  # Read target rlx field:
  Arlx = read_relax_piomas(dnmbR, pthrlx, varnm)

  # Read hindcast:
  dcice = os.path.join(pthhcst,f'ice_month.nc')
  print(f'Reading {dcice}')
  ds = xarray.open_dataset(dcice)
  match varnm:
    case('ithkn'):
      HIce = ds[ifld].isel(time=imo).data.squeeze()
      CIce = ds['siconc'].isel(time=imo).data.squeeze()
      Anep = HIce*CIce   
    case('iconc'):
      Anep = ds[ifld].isel(time=imo).data.squeeze()

  sqerr = (Anep - Arlx)**2
  nB = np.count_nonzero(~np.isnan(sqerr[(JB, IB)]))
  nA = np.count_nonzero(~np.isnan(sqerr[(JA, IA)]))
  rmseB = np.sqrt(np.nansum(sqerr[(JB, IB)]) / nB)
  rmseA = np.sqrt(np.nansum(sqerr[(JA, IA)]) / nA)
  biasB = np.nanmean(Anep[(JB,IB)] - Arlx[(JB,IB)])
  biasA = np.nanmean(Anep[(JA,IA)] - Arlx[(JA,IA)])

  RMSE_Ber.append(rmseB)
  RMSE_Arc.append(rmseA)
  BIAS_Ber.append(biasB)
  BIAS_Arc.append(biasA)
  TM.append(dnmbR)  

  print(f"RMSE: BerSea={rmseB:.2f} ArcReg={rmseA:.2f}")


nyr = int(len(RMSE_Ber)/12)
assert(nyr*12==len(RMSE_Ber)), f'Check record length RMSE_Ber={len(RMSE_Ber)}'

RMSE_Ber = np.array(RMSE_Ber)
RMSE_Arc = np.array(RMSE_Arc)
BIAS_Ber = np.array(BIAS_Ber)
BIAS_Arc = np.array(BIAS_Arc)
R2dB  = RMSE_Ber.reshape(nyr,12)
R2dA  = RMSE_Arc.reshape(nyr,12)
B2dB  = BIAS_Ber.reshape(nyr,12)
B2dA  = BIAS_Arc.reshape(nyr,12)

# Bering Sea
R2dB_md = np.median(R2dB, axis=0)
R2dB_pu = np.percentile(R2dB, 90, axis=0)
R2dB_pl = np.percentile(R2dB, 10, axis=0)

B2dB_md = np.median(B2dB, axis=0)
B2dB_pu = np.percentile(B2dB, 90, axis=0)
B2dB_pl = np.percentile(B2dB, 10, axis=0)

# Arctic
R2dA_md = np.median(R2dA, axis=0)
R2dA_pu = np.percentile(R2dA, 90, axis=0)
R2dA_pl = np.percentile(R2dA, 10, axis=0)

B2dA_md = np.median(B2dA, axis=0)
B2dA_pu = np.percentile(B2dA, 90, axis=0)
B2dA_pl = np.percentile(B2dA, 10, axis=0)

TM = np.array(TM)
DV = mtime.datevec1D(TM, fHR=False)
DV = np.array(DV).transpose()


plt.ion()
def plot_rmse_bias(fgnmb,time_yrs, rmse, rmse_md,rmse_pu,rmse_pl,\
                   bias, bias_md, bias_pu, bias_pl,clr1,clr2,regn,varnm):
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  # Time series RMSE
  time_yrs = DV[:,0]+(DV[:,1]-1)/12
  ax1 = plt.axes([0.06, 0.56, 0.4, 0.38])
  ax1.plot(time_yrs,rmse, '-', linewidth=2, color=clr1)
  sttl = f'NEPphys hcast GLORYS vs PIOMAS RMSE {varnm}\n {regn} {YRS}-{YRE}'
  ax1.set_title(sttl)

  # Monthly RMSE
  time_mnth = np.arange(1,13)
  ax2 = plt.axes([0.06,0.08,0.4,0.38])
  ax2.plot(time_mnth,rmse_md,'-',linewidth=2, color=clr1)
  ax2.plot(time_mnth,rmse_md, marker='o', markersize=7, color=clr1)
  ax2.plot(time_mnth,rmse_pu,'-',linewidth=1, color=clr2)
  ax2.plot(time_mnth,rmse_pl,'-',linewidth=1, color=clr2)
  ax2.set_xticks(time_mnth)
  ax2.grid('on')
  ax2.set_title(f'RMSE vs PIOMAS {varnm} {regn}, Median & IDR')
  ax2.set_xlabel('Months')

  # Time Series bias:
  ax3 = plt.axes([0.56, 0.56, 0.4, 0.38])
  ax3.plot(time_yrs,bias, '-', linewidth=2, color=clr1)
  sttl = f'NEPphys hcast vs PIOMAS Bias {varnm}\n {regn} {YRS}-{YRE}'
  ax3.set_title(sttl)

  # Monthly Bias
  ax4 = plt.axes([0.56,0.08,0.4,0.38])
  ax4.plot(time_mnth,bias_md,'-',linewidth=2, color=clr1)
  ax4.plot(time_mnth,bias_md, marker='o', markersize=7, color=clr1)
  ax4.plot(time_mnth,bias_pu,'-',linewidth=1, color=clr2)
  ax4.plot(time_mnth,bias_pl,'-',linewidth=1, color=clr2)
  ax4.set_xticks(time_mnth)
  ax4.grid('on')
  ax4.set_title(f'Bias vs PIOMAS {varnm} {regn}, Median & IDR')
  ax4.set_xlabel('Months')


  btx = 'calc_RMSE_ice_phys_hcast_PIOMAS.py'
  bottom_text(btx, pos=[0.02,0.01])

  return fig1, ax1, ax2, ax3, ax4


# -------------------
#
# Plot RMSE
#
# -------------------
clr_ber = [0,0.4,0.8]
clr2_ber = [0.7,0.8,1]
clr_arc = [0.,0.7,0.2]
clr2_arc = [0.7,1,0.9]
time_mnth = np.arange(1,13)

# Plot time series and monthly stat for RMSE and bias, Bering Sea:
fgnmb=1
time_yrs = DV[:,0]+(DV[:,1]-1)/12
fig1,ax11,ax12,ax13,ax14 = plot_rmse_bias(fgnmb,time_yrs,RMSE_Ber,R2dB_md,R2dB_pu,R2dB_pl,\
               BIAS_Ber,B2dB_md,B2dB_pu,B2dB_pl,clr_ber,clr2_ber,"Bering Sea",varnm)

# RMSE and bias for Arct.
fig3,ax31,ax32,ax33,ax34 = plot_rmse_bias(3,time_yrs,RMSE_Arc,R2dA_md,R2dA_pu,R2dA_pl,\
               BIAS_Arc,B2dA_md,B2dA_pu,B2dA_pl,clr_arc,clr2_arc,"Arctic Ocean",varnm)


# Save for plotting:
pthtmp = '/work/Dmitry.Dukhovskoy/anls_output/NEPbgc_hindcast02'
floutp = f'NEPphys_hcast_GLORYS_{varnm}_stat_{YRS}-{YRE}.npz' 
dflout = os.path.join(pthtmp,floutp)
print(f'Saving rmse bias arrays --> {dflout}')
np.savez(dflout, TM=TM, rmseB=RMSE_Ber, rmseA=RMSE_Arc, biasB=BIAS_Ber, biasA=BIAS_Arc)





