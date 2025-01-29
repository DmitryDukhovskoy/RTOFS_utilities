"""
  Calculate cold pool area region on the Bering Shelf from GOFS3.1 reanalysis
  for water T classes
  Monthly average bottom T for NEP region extracted in btmT_NEPmonthly_gofs31.py
  GOFS3.1 reanalysis 
  https://data.hycom.org/datasets/GLBv0.08/expt_53.X/data/2010/
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import time
import timeit
import pickle
from copy import copy
from yaml import safe_load

#PPTHN = '/home/Dmitry.Dukhovskoy/python'
PPTHN = []
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
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_anls_seas as manseas
import mod_plot_xsections as mxsct
import matplotlib as mtplt
import mod_regmom as mregmom
import mod_colormaps as mclrmps
import mod_gofs31 as mgofs
from mod_utils_fig import bottom_text

# Region domain
lat1 = 10.
lat2 = 81.
lon1 = 156.5
lon2 = 255.5

YAVRG = [x for x in range(2005,2015)]

pthoutp  = '/work/Dmitry.Dukhovskoy/GOFS3.1/gofs31_Tbtm_NEP/'

TCLASS = [-2., -1., 0., 1., 2.]
nTC    = len(TCLASS)

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

# Load Bering Shelf boundary, saved in plot_cold_pool_seas.py:
pthanls = pthseas['MOM6_NEP']['seasonal_daily']['pthanls'].format(expt_nmb=2)
dfnm_vrtx = os.path.join(pthanls,f'BeringShelf_region_vrtx.pkl')
if not os.path.isfile(dfnm_vrtx):
  raise Exception(f'Region boundary need to be saved in plot_cold_pool_seas.py {dfnm_vrtx}')

with open(dfnm_vrtx, 'rb') as fid:
  XYBND = pickle.load(fid)

Xbnd = XYBND[0]
Ybnd = XYBND[1]
Xbnd = np.where(Xbnd>180, Xbnd-360., Xbnd)

# Get GOFS3.1 topo, grid for NEP:
pthoutp  = '/work/Dmitry.Dukhovskoy/GOFS3.1/gofs31_Tbtm_NEP/'
pthSoutp = '/work/Dmitry.Dukhovskoy/GOFS3.1/gofs31_Sbtm_NEP/'
ftopo_regn = "gofs31_GLBv008_topo11_NEP.pkl"
dftopo_regn = os.path.join(pthoutp,ftopo_regn)
with open(dftopo_regn,'rb') as fid:
  lon1d, lat1d, HH = pickle.load(fid)
lon1d = np.where(lon1d>180, lon1d-360., lon1d)

jdm  = len(lat1d)
idm  = len(lon1d)
LONW = np.zeros((jdm,idm))
LATW = np.zeros((jdm,idm))
for ii in range(idm):
  LATW[:,ii]=lat1d
for jj in range(jdm):
  LONW[jj,:]=lon1d

# Derive Regional mask for Bering Shelf:
print('Searching indices of the Bering Shelf region ...')
II, JJ = mmisc.find_closest_indx(Xbnd, Ybnd, LONW, LATW)

DX, DY = mmom6.dx_dy(LONW, LATW)
Acell  = DX*DY
X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MS, _, _ = mmisc.inpolygon_v2(X, Y, II, JJ)  # 
JBS, IBS = np.where( (MS == 1) & (HH >= -250) & (HH < 0) ) #exclude deeep regions
MSKBS  = np.zeros((jdm,idm))
MSKBS[JBS,IBS] = 1

# use local bottom depth to compute pressure
jdm, idm = HH.shape
PR    = np.zeros((jdm,idm))
PR, _ = msw.sw_press(HH, LATW)

CPA  = np.zeros((len(YAVRG),12,nTC-1))
Time = []
iyr  = -1
for YR in (YAVRG):
  iyr += 1
  for MM in range(1,13):
    floutp  = f"gofs31_53X_Tbtm_NEP_{YR}{MM:02d}.pkl"
    flSoutp  = f"gofs31_53X_Sbtm_NEP_{YR}{MM:02d}.pkl"
    dfloutp = os.path.join(pthoutp,floutp)
    print(f'Loading {dfloutp}')
    with open(dfloutp,'rb') as fid:
      TbtmP = pickle.load(fid)
    dflSoutp = os.path.join(pthSoutp,flSoutp)
    print(f'Loading {dflSoutp}')
    with open(dflSoutp,'rb') as fid:
      Sbtm = pickle.load(fid)

    # Mask outside region:
    TbtmP = np.where( (MSKBS==0) & (HH<0) , 1.e3, TbtmP)
    # Check, should be empty:
    j0,i0 = np.where( (np.isnan(TbtmP)) & (HH < -10) )
    if len(j0) > 0:
      print(f'WARNING: {len(j0)} points Bottom T is missing')

    Sbtm = np.where( (MSKBS==0) & (HH<0) , 1.e3, Sbtm)

    SA = gsw.SA_from_SP(Sbtm, PR, LONW, LATW)

    # Compute conservative T from potential T
    print('Computing conservative T')
    Tbtm = gsw.CT_from_pt(SA, TbtmP)

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

pthanls = pthseas['MOM6_NEP']['seasonal_daily']['pthanls'].format(expt_nmb=2)
dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_gofs31_{YAVRG[0]}-{YAVRG[-1]}.pkl')
print(f'Dumping monthly coldpool area --> {dflout}')
with open(dflout, 'wb') as fid:
  pickle.dump([CPA, TCLASS],fid)

