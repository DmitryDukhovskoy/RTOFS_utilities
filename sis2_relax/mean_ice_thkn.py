"""
  Calc monthly mean ice thkn or vol/m2 (m3/m2=m) and total vol
 
  usage: mean_ice_thkn.py --yrs=1993 --ms=4 --yre=1994 --me=3 --expt=3
  jday = day of the year to plot
  or can provide date using month/day
 
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

parser = argparse.ArgumentParser()
parser.add_argument("--yrs", help="year start: 1993, ..., 2020", type=int)
parser.add_argument("--yre", help="year end", type=int)
parser.add_argument("--ms", help="month start: 1,...,12", type=int)
parser.add_argument("--me", help="month end: 1, ..., 12", type=int)
parser.add_argument("--expt", help="experiment number: 1, ...", type=int)
args = parser.parse_args()

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
#varnm  = 'ithkn'  # iconc or ithkn
#varnm  = 'iconc'

# Start of the run - needed only for seasonal forecasts:
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1    
nens   = 1    # ens # for ensemble runs

expt_nmb = 2  

if args.yrs:
  YRS = args.yrs
if args.ms:
  MS = args.ms
if args.yre:
  YRE = args.yre
if args.me:
  ME = args.me
if args.expt:
  expt_nmb = args.expt

# Choose experiment:
expt     = 'test_ice_relax'
#expt     = "seasonal_daily"
#expt     = 'test'   # saved output during test runs
#
if expt == 'test_ice_relax':
  runname = expt
elif expt == 'test':
  runname  = 'isponge_test'
elif expt == 'seasonal_daily':
  expt_nmb = 2  # only 1 experiment 
  runname  = f"NEPphys_frcst_dailyOB-expt{expt_nmb:02d}"

#
expt_nmb0 = f"{expt_nmb:02d}"


fyaml_param='relax_expts.yaml'
with open(fyaml_param) as ff:
  param_expt = safe_load(ff)  

dt_idyn = param_expt[expt][expt_nmb0]['dt_idyn']
dt_slow = param_expt[expt][expt_nmb0]['dt_slow']
rlx_max = param_expt[expt][expt_nmb0]['rlx_max']

print(f'Expt: {expt} Run: {runname} Means for {YRS}/{MS} - {YRE}/{ME}')

if expt == 'NEP_BGCphys_GOFS':
  outfld = 'ice'
else:
  outfld = 'icem'

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
DX, DY = mmom6.dx_dy(hlon, hlat)
Acell = DX*DY

# Define mask for computing mean values:
lat_cut = 65.
LMsk = np.where(hlat>=lat_cut, 1., 0.)
area_sub = np.sum(Acell*LMsk)     # total area of the domain for mean calc. 

VolMean = []
ThknMean = []
icc = 0
for YR in range(YRS, YRE+1):
  mstart = MS
  if YR > YRS:
    mstart = 1
  mend = 12
  if YR == YRE:
    mend = ME
    
  for month in range(mstart, mend+1): 
    thkn_mnth = 0.
    vol_mnth = 0.
    imm = 0
    for mday in range(1,28,5):
      print(f"Reading {YR}/{month}/{mday}")
      dnmbR  = mtime.datenum([YR, month, mday])  
      if expt == 'seasonal_daily':
        pth1     = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=expt_nmb)
        dir_fcst = pthseas['MOM6_NEP'][expt]['dir_icefcst'].format(\
             yr_start=YRS, mo_start=MOS, ens=nens, yr_run=YR, mo_run=month)
        pthfcst = os.path.join(pth1,dir_fcst)
      elif expt == 'test_ice_relax':
        pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YY=YRS, MM=MOS, expt_nmb=expt_nmb)
      else:
        pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YY=YR, MM=month)

      # Find closest output:
      YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthfcst, dnmbR, fld=outfld)
      dv0  = mtime.datevec(dnmb0)
      YR0, MM0, DD0 = dv0[:3]


      flice_name = pthseas['MOM6_NEP'][expt]['ficename'].format(YR=YR0, jday=jday0)
      dfsis2 = os.path.join(pthfcst, flice_name)
      dset   = xarray.open_dataset(dfsis2)
      HIce = dset['sithick'].isel(time=0).data
      CIce = dset['siconc'].isel(time=0).data
      VIce = CIce*HIce  # ice vol m3/m2

      VIce = np.where(hlat < 65., 0., VIce)
      icc += 1
      imm += 1

      vol_mn  = np.nansum(VIce*Acell*LMsk)
      thkn_mn = vol_mn/area_sub
      vol_mnth  = vol_mnth + vol_mn
      thkn_mnth = thkn_mnth + thkn_mn

    vol_mnth  = vol_mnth/imm*1.e-9  # km3
    thkn_mnth = thkn_mnth/imm 
    VolMean.append(vol_mnth)
    ThknMean.append(thkn_mnth)

    print(f"YR={YR} MO={month} recs={imm} mean vol={vol_mnth:.1f}km3 thkn={thkn_mnth:.2}m")


plt.ion()

xmnths=np.arange(12)+MS

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.55, 0.8, 0.35])
ax1.plot(xmnths,ThknMean,'-o')
sttl1 = f"{expt}-{expt_nmb0} ThknMean (m), {YRS}/{MS}-{YRE}/{ME}"
ax1.set_title(sttl1)
ax1.set_xticks(xmnths)
ax1.grid('on')

ax2 = plt.axes([0.1, 0.1, 0.8, 0.35])
ax2.plot(xmnths,VolMean,'-o')
sttl1 = f"{expt}-{expt_nmb0} VolMean (km3), {YRS}/{MS}-{YRE}/{ME}"
ax2.set_title(sttl1)
ax2.set_xticks(xmnths)
ax2.grid('on')

btx = 'mean_ice_thkn.py'
bottom_text(btx, pos=[0.2, 0.01])





