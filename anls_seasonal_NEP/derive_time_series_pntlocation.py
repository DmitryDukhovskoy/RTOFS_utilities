"""
  Derive time series of  bottom T or S at specified locaitons
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


# Initial date
# Look at ens run #1 - the only ens. that has 5-day av. output fields
expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
varnm    = 'salin'  # temp (potential) / salin
#dnmbS    = mtime.datenum([2015,1,1])
nensR    = 1
expt_nmb = 2   # 2 - seas f/casts with dailyOB, #3 - seas f/cast with ice relaxation
# Averaging time period:
MMS   = 4    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YRS = 1993
YRE = 1999
YAVRG = [x for x in range(YRS,YRE+1)]

#YAVRG = [x for x in range(2005,2006)]
#MAVRG = [1,2,3]  # months to average: Winter  JFM, Summer: JAS
#MAVRG = [4,5,6]  # months to average: Spring, AMJ
#MAVRG = [7,8,9]  # months to average: Winter  JFM, Summer: JAS
#MAVRG = [10,11,12]  # months to average: Fall

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

#DX, DY = mmom6.dx_dy(hlon, hlat)
#Acell  = DX*DY

ocnfld = 'oceanm'
pthoutp0 = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=expt_nmb)
YR=YAVRG[0]
MM=MMS
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

TM = []
icc = 0
for YR in (YAVRG):
  for mnth in range(1,13):
    dnmbS    = mtime.datenum([YR,MMS,1])
    pthfcst0 = os.path.join(pthoutp0,f'{YR}-{MMS:02d}-e{nensR:02d}','history')
    TT, time = manseas.monthly_mean_from_Ndaily_ocean3D(pthfcst0, YR, MMS, 'temp', \
                    ocnfld, MAVRG=[mnth], mnth='fcst')
    SS, time = manseas.monthly_mean_from_Ndaily_ocean3D(pthfcst0, YR, MMS, 'salt', \
                    ocnfld, MAVRG=[mnth], mnth='fcst')

    TM.append(time)

    Tb = manseas.derive_bottom_temp(TT, dP)
    Sb = manseas.derive_bottom_temp(SS, dP)

    tb1 = Tb[JP,IP]
    tb1 = np.expand_dims(tb1, axis=0)
    sb1 = Sb[JP,IP]
    sb1 = np.expand_dims(sb1, axis=0)

    if icc==0:
      TBTM = tb1.copy()
      SBTM = sb1.copy()
    else:
      TBTM = np.append(TBTM, tb1, axis=0)
      SBTM = np.append(SBTM, sb1, axis=0)

    icc += 1

btx = 'derive_time_series_pntlocation.py'


import pickle
pthanls = pthseas['MOM6_NEP'][expt]['pthanls'].format(expt_nmb=expt_nmb)
dfnm_tts = os.path.join(pthanls,f'mnthly_Tbtm_tser_{YRS}{MMS:02d}.pkl')
dfnm_sts = os.path.join(pthanls,f'mnthly_Sbtm_tser_{YRS}{MMS:02d}.pkl')

# Save for analysis using WOA23 or other fields:
print(f'Dumping Tbtm Time ser  --> {dfnm_tts}')
with open(dfnm_tts, 'wb') as fid:
  pickle.dump([TBTM,TM],fid)

print(f'Dumping Sbtm Time ser  --> {dfnm_sts}')
with open(dfnm_sts, 'wb') as fid:
  pickle.dump([SBTM,TM],fid)


