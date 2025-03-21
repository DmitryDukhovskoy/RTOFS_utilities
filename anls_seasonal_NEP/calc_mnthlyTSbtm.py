"""
  Calc monthly mean and StDev  bottom T from seasonal forecasts
  Save only northern part of the domain (north of Jcut)
  Save by years

  Usage: calc_mnthlyTbtm_BerSea.py --YRS=1993 --YRE=2000 --MMI=4 --expt=3
  use --help for more information on keywargs

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
import argparse
import pickle

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

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="Start f/cast init year to save averaged data", type=int)
parser.add_argument("--YRE", help="End f/cast init year to save averaged data", type=int)
parser.add_argument("--MMI", help="Forecast initialization month ", type=int)
parser.add_argument("--expt", help="F/casts experiment number", type=int)
args = parser.parse_args()


# experiments: 2 - daily OB seasonal forecasts, 3 - same as 2 but with sea ice relaxation
# Default values: 
# Look at ens run #1 - the only ens. that has 5-day av. output fields
expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
nensR    = 1
expt_nmb = 2   # 2 - seas f/casts with dailyOB, #3 - seas f/cast with ice relaxation
# Averaging time period:
MMI   = 4    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YRS = 1993   # Start: f/cast init. year to use for monthly averaging
YRE = YRS   # End

if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI
if args.expt:
  expt_nmb = args.expt 

if YRS == 1993 and MMI == 1:
  raise Exception("First initial month should be 4 for 1993, given MMI={MMI}")

YAVRG = [x for x in range(YRS,YRE+1)]
#Jcut = 500 # chop off lower latitudes

expt_name = f'NEPphys_frcst_dailyOB-expt{expt_nmb:02d}'
run_info = f'{expt_name} init MM={MMI} e{nensR:02d}, T/S bottom: {min(YAVRG)}-{max(YAVRG)}'

print(f'Deriving T/S bottom stat for {expt_name} ')
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

#DX, DY = mmom6.dx_dy(hlon, hlat)
#Acell  = DX*DY

ocnfld = 'oceanm'
pthoutp0 = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=expt_nmb)
YR=YAVRG[0]
MM=MMI
subdir=f'oceanm_{YR}{MM:02d}'
pthfcst0 = os.path.join(pthoutp0,f'{YR}-{MM:02d}-e01','history')
list_files = manseas.list_oceanice_files(pthfcst0, prefix=ocnfld, subdir=subdir)
pthfull = os.path.join(pthfcst0,subdir)
floceanm = list_files[0]
ZM = manseas.read_oceanm3D_field(pthfull, floceanm, 'zl', notime=False)
ZM = -abs(ZM)
dP = manseas.read_oceanm3D_field(pthfull, floceanm, 'h', notime=False)
dP = np.where(dP < 1.e-3, 0., dP)
ZZ = mmom6.zm2zz(ZM)

for YR in (YAVRG):
  TM = []
  icc = 0
  for mnth in range(1,13):
    print(f"Processing {YR}/{mnth} ...")

    dnmbS    = mtime.datenum([YR,MMI,1])
    pthfcst0 = os.path.join(pthoutp0,f'{YR}-{MMI:02d}-e{nensR:02d}','history')
    Tb_mn, Tb_std, time = manseas.monthly_TSbtm_Ndaily_ocean3D(pthfcst0, YR, MMI, 'temp', \
                                     ocnfld, MAVRG=[mnth], mnth='fcst')

    Sb_mn, Sb_std, time = manseas.monthly_TSbtm_Ndaily_ocean3D(pthfcst0, YR, MMI, 'salt', \
                                     ocnfld, MAVRG=[mnth], mnth='fcst')

    time_mn = int(np.mean(np.array(time)))
    TM.append(time)

    Tb_mn  = np.expand_dims(Tb_mn, axis=0)
    Tb_std = np.expand_dims(Tb_std, axis=0)
    Sb_mn  = np.expand_dims(Sb_mn, axis=0)
    Sb_std = np.expand_dims(Sb_std, axis=0)

    if icc==0:
      TBTM_mn  = Tb_mn.copy()
      TBTM_std = Tb_std.copy()
      SBTM_mn  = Sb_mn.copy()
      SBTM_std = Sb_std.copy()
    else:
      TBTM_mn  = np.append(TBTM_mn, Tb_mn, axis=0)
      TBTM_std = np.append(TBTM_std, Tb_std, axis=0)
      SBTM_mn  = np.append(SBTM_mn, Sb_mn, axis=0)
      SBTM_std = np.append(SBTM_std, Sb_std, axis=0)

    icc += 1

  pthanls = pthseas['MOM6_NEP'][expt]['pthanls'].format(expt_nmb=expt_nmb)
  dfnm_tts = os.path.join(pthanls,f'mnthly_Tbtm_expt{expt_nmb:02d}_{YR}{MMI:02d}.pkl')
  dfnm_sts = os.path.join(pthanls,f'mnthly_Sbtm_expt{expt_nmb:02d}_{YR}{MMI:02d}.pkl')

  # Save for analysis using WOA23 or other fields:
  print(f'Dumping T bottom stat  --> {dfnm_tts}')
  with open(dfnm_tts, 'wb') as fid:
    pickle.dump([TBTM_mn, TBTM_std, TM], fid)

  print(f'Dumping S bottom stat  --> {dfnm_sts}')
  with open(dfnm_sts, 'wb') as fid:
    pickle.dump([SBTM_mn, SBTM_std, TM], fid)


