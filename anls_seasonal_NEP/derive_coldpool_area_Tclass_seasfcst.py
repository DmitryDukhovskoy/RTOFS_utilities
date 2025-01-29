"""
  Compute area cold pool on the Bering Sea shelf
  bottom water < 2C for T classes

  Plot bottom T and cold pool in the Bering Sea by seasons:
  bottom water < 2C
  e.g. On the variability of the Bering Sea Cold Pool and implications 
       for the biophysical environment
  2022
 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8979450/

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
# Averaging time period:
MMS   = 1    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YAVRG = [x for x in range(2005,2015)]

nensR    = 1
expt_nmb = 2   # 2 - seas f/casts with dailyOB

expt_name = f'NEPphys_frcst_dailyOB-expt{expt_nmb:02d}'
run_info = f'{expt_name} init MM={MMS} e{nensR:02d}, conservT bottom: {min(YAVRG)}-{max(YAVRG)}' 

print(f'Plotting {varnm} {expt_name} ')
print(f'{run_info}')

TCLASS = [-2., -1., 0., 1., 2.]
nTC    = len(TCLASS)

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

# Get indices of the polygon:
II = pthseas['ANLS_NEP']['poly_BerSea']['II']
JJ = pthseas['ANLS_NEP']['poly_BerSea']['JJ']
jdm, idm = HH.shape

DX, DY = mmom6.dx_dy(hlon, hlat)
Acell  = DX*DY
X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MS, _, _ = mmisc.inpolygon_v2(X, Y, II, JJ)  # 
JBS, IBS = np.where( (MS == 1) & (HH >= -250) & (HH < 0) ) #exclude deeep regions
MSKBS  = np.zeros((jdm,idm))
MSKBS[JBS,IBS] = 1

ocnfld = 'oceanm'
pthoutp0 = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=expt_nmb)
YR=2011
MM=4
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

# Derive pressure at ZM depths
jdm, idm = HH.shape
kdm   = len(ZM)
Z3d   = np.tile(ZM, idm*jdm).reshape((idm,jdm,kdm))
Z3d   = np.transpose(Z3d, (2, 1, 0))
PR    = np.zeros((kdm,jdm,idm))
for kk in range(kdm):
  pr_db, _ = msw.sw_press(Z3d[kk,:,:].squeeze(), hlat)
  PR[kk,:] = pr_db


CPA  = np.zeros((len(YAVRG),12,nTC-1))
Time = []
iyr  = -1
for YRS in (YAVRG):
  iyr += 1
  dnmbS    = mtime.datenum([YRS,MMS,1])
  pthfcst0 = os.path.join(pthoutp0,f'{YRS}-{MMS:02d}-e{nensR:02d}','history')
  for MM in range(1,13):
    TT, _ = manseas.monthly_mean_from_Ndaily_ocean3D(pthfcst0, YRS, MMS, 'temp', ocnfld, MAVRG=[MM])
    SS, _ = manseas.monthly_mean_from_Ndaily_ocean3D(pthfcst0, YRS, MMS, 'salt', ocnfld, MAVRG=[MM])

    Tbtm = manseas.derive_conservTbtm_from_T3d(TT, SS, PR, dP, hlon, hlat)
    # Mask outside region:
    Tbtm = np.where( (MSKBS==0) & (HH<0) , 1.e3, Tbtm)
    # Check, should be empty:
    j0,i0 = np.where( (np.isnan(Tbtm)) & (HH < -10) )
    if len(j0) > 0:
      print(f'WARNING: {len(j0)} points Bottom T is missing')

    for icl in range(1,nTC):
      t1 = TCLASS[icl-1]
      t2 = TCLASS[icl]

      JBS, IBS = np.where( (MSKBS == 1) & (Tbtm <= t2) & (Tbtm > t1))
      if len(JBS) > 0.:
        CP_area = np.sum(Acell[JBS,IBS])*1e-6 # km2
      else:
        CP_area = 0.
      print(f'TCLASS: {t1:.1f} - {t2:.1f} Area={CP_area*1e-5:.1f} x1e5 km2')

      CPA[iyr, MM-1, icl-1] = CP_area


pthanls = pthseas['MOM6_NEP'][expt]['pthanls'].format(expt_nmb=expt_nmb)
dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_fcst_{YAVRG[0]}-{YAVRG[-1]}.pkl')
print(f'Dumping monthly coldpool area --> {dflout}')
with open(dflout, 'wb') as fid:
  pickle.dump([CPA, TCLASS],fid)


