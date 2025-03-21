"""
  Compute xcorrelation: time series of  bottom T or S at specified locaitons
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import pickle
from copy import copy
import matplotlib.colors as colors
from yaml import safe_load

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
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')
import mod_swstate as msw
import conversions as gsw

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
importlib.reload(manseas)
import pickle


# Initial date
# Look at ens run #1 - the only ens. that has 5-day av. output fields
expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
varnm    = 'salin'  # temp (potential) / salin
#dnmbS    = mtime.datenum([2015,1,1])
nensR    = 1
expt_nmb = 3   # 2 - seas f/casts with dailyOB, #3 - seas f/cast with ice relaxation
# Averaging time period:
MMS   = 4    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YRS = 1993
YRE = 1999
YAVRG = [x for x in range(YRS,YRE+1)]


if YRS == 1993 and MMS == 1:
  raise Exception("First initial month should be 4 for 1993, given MMS={MMS}")

IP = [158, 164, 148, 156, 118, 120, 128, 75,  65,  201]
JP = [698, 648, 702, 647, 714, 681, 610, 671, 744, 234]

expt_name = f'NEPphys_frcst_dailyOB-expt{expt_nmb:02d}'
run_info = f'{expt_name} init MM={MMS} e{nensR:02d}, T/S bottom: {min(YAVRG)}-{max(YAVRG)}'

print(f'Plotting {varnm} {expt_name} ')
print(f'{run_info}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

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

jdm, idm = HH.shape

pthanls = pthseas['MOM6_NEP'][expt]['pthanls'].format(expt_nmb=2)
dfnm_tts = os.path.join(pthanls,f'mnthly_Tbtm_tser_{YRS}{MMS:02d}.pkl')
dfnm_sts = os.path.join(pthanls,f'mnthly_Sbtm_tser_{YRS}{MMS:02d}.pkl')

with open(dfnm_tts, 'rb') as fid:
  TBTM2, TM = pickle.load(fid)

print(f"Loading {dfnm_sts}")
with open(dfnm_sts, 'rb') as fid:
  SBTM2, TM2 = pickle.load(fid)

pthanls = pthseas['MOM6_NEP'][expt]['pthanls'].format(expt_nmb=3)
dfnm_tts = os.path.join(pthanls,f'mnthly_Tbtm_tser_{YRS}{MMS:02d}.pkl')
dfnm_sts = os.path.join(pthanls,f'mnthly_Sbtm_tser_{YRS}{MMS:02d}.pkl')

with open(dfnm_tts, 'rb') as fid:
  TBTM3, TM = pickle.load(fid)

print(f"Loading {dfnm_sts}")
with open(dfnm_sts, 'rb') as fid:
  SBTM3, TM3 = pickle.load(fid)


Xtm = np.arange(1,len(TM)+1)
btx = 'xcorr_timesers.py'
plt.ion()
fig1 = plt.figure(1,figsize=(9,8))

ii0 = 6
S2 = SBTM2[:,ii0]
S3 = SBTM3[:,ii0+1]

S2 = S2-np.mean(S2)
S3 = S3-np.mean(S3)

ii0=9
T2 = TBTM2[:,ii0]
T3 = TBTM3[:,ii0]
#T2 = T2-np.mean(T2)
#T3 = T3-np.mean(T3)

plt.clf()

ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
ax1.plot(Xtm, T2)
ax1.plot(Xtm, T3)

mtck = [x for x in range(1,len(TM)+1,12)]
ax1.set_xticks(mtck)
ax1.grid('on')

Ib = IP[ii0]
Jb = JP[ii0]
sttl = f'Bottom T, i={Ib}, j={Jb}'
ax1.set_title(sttl)


