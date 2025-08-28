"""
  Calc and plot RMSE for monthly ssh
  from BGC and phys hindcasts

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
args = parser.parse_args()

YRS   = args.yrs if args.yrs else None
YRE   = args.yre if args.yre else YRS
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
  for mm in range(1,13):
    MCAL.append(mm)
    YCAL.append(yr)

len_cal = len(MCAL)  
  
icc = 0
TM  = [] 
RMSE = []
for ii in range(len(MCAL)):
  MM = MCAL[ii]
  YY = YCAL[ii] 

  # Find init date for given month, assuming hcst_time (n months) f/cast interval
  kint = np.searchsorted(hcst_interv, MM, side='right') - 1
  assert(hcst_interv[kint] <= MM < hcst_interv[kint+1]), f'Wrong time bin {kint} for {MMA}'
  MINIT = hcst_interv[kint]
  imo = MM-MINIT      # current month in the archive output

  pthhcst_bgc = (
       pthseas['MOM6_NEP']['hindcast']['pthoutp'].format(
         hindcast_name=hcast_bgc, YR=YY, MM=MINIT, DD=1
    )
  )

  pthhcst_phys = '/archive/Dmitry.Dukhovskoy/fre/NEP/2024/NEP_physics_202404_nudging-15d/' + \
                 f'gfdl.ncrc5-intel22-repro/history/{YY}-{MINIT:02d}'


  print(f'Calculating RMSE for {YY}/{MM}')
  dnmbR = mtime.datenum([YY,MM,15])

  # Read hindcast:
  docn_bgc = os.path.join(pthhcst_bgc,f'ocean_month.nc')
  print(f'Reading {docn_bgc}')
  with xarray.open_dataset(docn_bgc) as ds_bgc:
    ssh_bgc = ds_bgc['ssh'].isel(time=imo).data.squeeze()

  docn_phys = os.path.join(pthhcst_phys,f'ocean_month.nc')
  print(f'Reading {docn_phys}')
  with xarray.open_dataset(docn_phys) as ds_phys:
    ssh_phys = ds_phys['ssh'].isel(time=imo).data.squeeze()

  sqerr = (ssh_bgc-ssh_phys)**2
  nB = np.count_nonzero(~np.isnan(sqerr))
  rmseBP = np.sqrt(np.nansum(sqerr)/nB)

  RMSE.append(rmseBP)
  TM.append(dnmbR)  

  print(f"RMSE: {rmseBP:.2f}")


nyr = int(len(RMSE)/12)
assert(nyr*12==len(RMSE)), f'Check record length RMSE_Ber={len(RMSE_Ber)}'

RMSE = np.array(RMSE)
R2dB  = RMSE_Ber.reshape(nyr,12)

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
  sttl = f'NEPbgc hcast GLORYS+irlx vs PIOMAS RMSE {varnm}\n {regn} {YRS}-{YRE}'
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
  sttl = f'NEPbgc hcast vs PIOMAS Bias {varnm}\n {regn} {YRS}-{YRE}'
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


  btx = 'calc_RMSE_ice_bgc_hcast_PIOMAS.py'
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
floutp = f'NEPbgc_hcast_GLORYSirlx_{varnm}_stat_{YRS}-{YRE}.npz'
dflout = os.path.join(pthtmp,floutp)
print(f'Saving rmse bias arrays --> {dflout}')
np.savez(dflout, TM=TM, rmseB=RMSE_Ber, rmseA=RMSE_Arc, biasB=BIAS_Ber, biasA=BIAS_Arc)

