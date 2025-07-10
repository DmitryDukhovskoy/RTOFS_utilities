"""
  Plot RMSE and Bias statistics from NEPbgc and NEPphys hindcasts
  to see the impact of ice relaxation
  Statistics computed wrt PIOMAS ice fields used as target relaxation fields
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
parser.add_argument("--yrs", help="year start to calc RMSE: 1993, ..., 2020", type=int)
parser.add_argument("--yre", help="year end RMSE", type=int)
parser.add_argument("--varnm", help="field to use: ithkn or iconc", type=str)
args = parser.parse_args()

varnm = args.varnm if args.varnm else 'iconc'
YRS   = args.yrs if args.yrs else 1993
YRE   = args.yre if args.yre else 2019


# Load data:
pthtmp = '/work/Dmitry.Dukhovskoy/anls_output/NEPbgc_hindcast02'
floutp1 = f'NEPbgc_hcast_GLORYSirlx_{varnm}_stat_{YRS}-{YRE}.npz'
floutp2 = f'NEPphys_hcast_GLORYS_{varnm}_stat_{YRS}-{YRE}.npz'
dflout1 = os.path.join(pthtmp,floutp1)
dflout2 = os.path.join(pthtmp,floutp2)
print(f'Loading rmse bias arrays --> {dflout1}')
data      = np.load(dflout1)
TM1       = data['TM']
RMSE_Ber1 = data['rmseB']
RMSE_Arc1 = data['rmseA']
BIAS_Ber1 = data['biasB']
BIAS_Arc1 = data['biasA']

print(f'Loading rmse bias arrays --> {dflout2}')
data      = np.load(dflout2)
TM2       = data['TM']
RMSE_Ber2 = data['rmseB']
RMSE_Arc2 = data['rmseA']
BIAS_Ber2 = data['biasB']
BIAS_Arc2 = data['biasA']

# Check record length:
assert (
    len(TM1) == len(TM2) and TM1[0] == TM2[0]
), f'Time series mismatch: len={len(TM1)} vs {len(TM2)}, start={TM1[0]} vs {TM2[0]}'

nyr = int(len(RMSE_Ber1)/12)
assert nyr*12==len(RMSE_Ber1), f'Check record length RMSE_Ber={len(RMSE_Ber)}'

def get_median_prct(R1d,lprc):
  nyr = int(len(R1d)/12)
  assert nyr*12==len(R1d), f'Check record length {len(R1d)}'

  R2d = R1d.reshape(nyr,12)
  md = np.median(R2d, axis=0)
  pu = np.percentile(R2d, 100-lprc, axis=0)
  pl = np.percentile(R2d, lprc, axis=0)

  return md, pu, pl

def plot_rmse_bias(fgnmb, rmse1, rmse2, bias1, bias2, regn, varnm, TM):
  clr1  = [0.,0.4,0.8]
  clr11 = [0.85,0.95,1]
  clr2  = [0.8,0.4,0]
  clr21 = [1,0.95,0.85]

  DV = mtime.datevec1D(TM, fHR=False)
  DV = np.array(DV).transpose()
  time_yrs = DV[:,0]+(DV[:,1]-1)/12

  Rmd1, Ruprc1, Rlprc1 = get_median_prct(rmse1,10)
  Rmd2, Ruprc2, Rlprc2 = get_median_prct(rmse2,10)
  Bmd1, Buprc1, Blprc1 = get_median_prct(bias1,10)
  Bmd2, Buprc2, Blprc2 = get_median_prct(bias2,10)

  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  # Time series RMSE
  time_yrs = DV[:,0]+(DV[:,1]-1)/12
  ax1 = plt.axes([0.06, 0.56, 0.4, 0.38])
  ax1.plot(time_yrs,rmse1, '-', linewidth=2, color=clr1)
  ax1.plot(time_yrs,rmse2, '-', linewidth=2, color=clr2)
  ax1.grid('on')
  sttl = f'NEP hcasts vs PIOMAS RMSE {varnm}\n {regn} {YRS}-{YRE}'
  ax1.set_title(sttl)

  # Monthly RMSE
  time_mnth = np.arange(1,13)
  ax2 = plt.axes([0.06,0.08,0.4,0.38])
  ax2.plot(time_mnth,Rmd1,'-',linewidth=2, color=clr1)
  ax2.plot(time_mnth,Rmd1, marker='o', markersize=7, color=clr1)
  ax2.plot(time_mnth,Ruprc1,'-',linewidth=1, color=clr11)
  ax2.plot(time_mnth,Rlprc1,'-',linewidth=1, color=clr11)

  ax2.plot(time_mnth,Rmd2,'-',linewidth=2, color=clr2)
  ax2.plot(time_mnth,Rmd2, marker='o', markersize=7, color=clr2)
  ax2.plot(time_mnth,Ruprc2,'-',linewidth=2, color=clr21)
  ax2.plot(time_mnth,Rlprc2,'-',linewidth=2, color=clr21)

  ax2.set_xticks(time_mnth)
  ax2.grid('on')
  ax2.set_title(f'RMSE hcasts vs PIOMAS {varnm} {regn}')
  ax2.set_xlabel('Months')

  # Time Series bias:
  ax3 = plt.axes([0.56, 0.56, 0.4, 0.38])
  ax3.plot(time_yrs,bias1, '-', linewidth=2, color=clr1)
  ax3.plot(time_yrs,bias2, '-', linewidth=2, color=clr2)
  ax3.grid('on')
  sttl = f'NEPbgc hcasts vs PIOMAS Bias {varnm}\n {regn} {YRS}-{YRE}'
  ax3.set_title(sttl)

  # Monthly Bias
  ax4 = plt.axes([0.56,0.08,0.4,0.38])
  ax4.plot(time_mnth,Bmd1,'-',linewidth=2, color=clr1)
  ax4.plot(time_mnth,Bmd1, marker='o', markersize=7, color=clr1)
  ax4.plot(time_mnth,Buprc1,'-',linewidth=1, color=clr11)
  ax4.plot(time_mnth,Blprc1,'-',linewidth=1, color=clr11)

  ax4.plot(time_mnth,Bmd2,'-',linewidth=2, color=clr2)
  ax4.plot(time_mnth,Bmd2, marker='o', markersize=7, color=clr2)
  ax4.plot(time_mnth,Buprc2,'-',linewidth=2, color=clr21)
  ax4.plot(time_mnth,Blprc2,'-',linewidth=2, color=clr21)

  ax4.set_xticks(time_mnth)
  ax4.grid('on')
  ax4.set_title(f'Bias hcasts vs PIOMAS {varnm} {regn}')
  ax4.set_xlabel('Months')

  ax5 = plt.axes([0.46,0.47,0.14,0.06])
  ax5.plot([0.1,0.3],[0.2,0.2],'-',linewidth=2, color=clr1)
  ax5.text(0.35,0.18,expt1)
  ax5.plot([0.1,0.3],[0.3,0.3],'-',linewidth=2, color=clr2)
  ax5.text(0.35,0.28,expt2)
  ax5.set_xlim([0,1])
  ax5.set_ylim([0.15,0.35])
  ax5.axis('off')

  btx = 'plot_RMSE_ice_hcasts_PIOMAS.py'
  bottom_text(btx, pos=[0.02,0.01])

  return fig1, ax1, ax2, ax3, ax4


# -------------------
#
# Plot RMSE
#
# -------------------
# Plot time series and monthly stat for RMSE and bias, Bering Sea:
plt.ion()

expt1 = 'NEPbgc'
expt2 = 'NEPphys'
fgnmb=1
fig1,ax11,ax12,ax13,ax14 = plot_rmse_bias(fgnmb, RMSE_Ber1, RMSE_Ber2, \
                    BIAS_Ber1, BIAS_Ber2, "BerSea", varnm, TM1)

# RMSE and bias for Arct.
fig3,ax31,ax32,ax33,ax34 = plot_rmse_bias(2, RMSE_Arc1, RMSE_Arc2, \
                            BIAS_Arc1, BIAS_Arc2, "Arctic", varnm, TM1)



