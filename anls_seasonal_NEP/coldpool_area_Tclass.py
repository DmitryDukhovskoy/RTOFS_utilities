"""
  Compute area cold pool in the Bering Sea
  bottom water < 2C for T classes

  deep water (> 250 m) is excluded
  following:
  On the variability of the Bering Sea Cold Pool and implications 
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

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
f_insitu = True    # convert to in situ 

# Start of the run - needed only for seasonal forecasts:
YRS    = 1993 # year start of the forecast
MMS    = 4
DDS    = 1
nens   = 2    # ens # for ensemble runs

#expt    = "seasonal_fcst"
#runname = f'NEPphys_frcst_climOB_{YRS}-{MMS:02d}-e{nens:02d}'
#expt    = 'NEP_BGCphys_GOFS'
#runname = 'NEP_physics_GOFS-IC'
expt    = 'NEP_seasfcst_LZRESCALE'
runname = 'NEPphys_LZRESCALE_climOB_1993_04-e02'

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

if not expt == 'seasonal_fcst':
  YRS =  pthseas['MOM6_NEP'][expt]['year_start'] 
  MMS =  pthseas['MOM6_NEP'][expt]['month_start'] 
  DDS =  pthseas['MOM6_NEP'][expt]['day_start'] 

TCLASS = [-2., -1., 0., 1., 2.]
nTC    = len(TCLASS)

dnmbS  = mtime.datenum([YRS,MMS,DDS])  # day to plot
dnmbE  = dnmbS + 365
TPLT = mmom6.create_time_array(dnmbS, dnmbE, 5)
nT   = len(TPLT)
CPA  = np.zeros((nT,nTC-1))*np.nan
icc  = 0
for kk in range(nT):
  dnmbR = TPLT[kk][0]
  dvR = mtime.datevec(dnmbR)
  print(f'Expt: {expt} Run: {runname} Plot date: {dvR[0]}/{dvR[1]}/{dvR[2]}')

  if expt == 'seasonal_fcst':
    pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(runname=runname)
  else:
    dnmb0    = dnmbR
    dv0      = mtime.datevec(dnmb0)
    YR0, MM0, DD0 = dv0[:3]
    jday0    = int(mtime.date2jday([YR0,MM0,DD0]))
    pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YY=YR0, MM=MM0)
  
  if icc == 0:
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

    TM  = []

# Find closest output:
  if expt == 'NEP_BGCphys_GOFS':
    ocnfld = 'ocean'
  else:
    ocnfld = 'oceanm'
  if expt == 'seasonal_fcst':
    pthfcst = os.path.join(pthfcst,f'{ocnfld}_{dvR[0]}{dvR[1]:02d}')
   
  if not os.path.isdir(pthfcst):
    print(f'not exist: {pthfcst}')
    continue
 
  YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthfcst, dnmbR, fld=ocnfld)
  dv0  = mtime.datevec(dnmb0)
  YR0, MM0, DD0 = dv0[:3]

  flocn_name = pthseas['MOM6_NEP'][expt]['focname'].format(YR=YR0, jday=jday0)
  dfmom6 = os.path.join(pthfcst, flocn_name)

  if not os.path.isfile(dfmom6): 
    print(f'{dfmom6} not found ...')
    continue

# Check if date not repeated
  if (icc >0) and (dnmb0 == TM[-1]): continue
  
  TM.append(dnmb0)

  dset = xarray.open_dataset(dfmom6)
  ZM   = -dset['zl'].data
  T3d  = dset['potT'].data[0,:].squeeze()
  S3d  = dset['salt'].data[0,:].squeeze()
  dP   = dset['h'].data[0,:].squeeze()
  dP   = np.where(dP < 1.e-3, 0., dP)
  ZZ   = mmom6.zm2zz(ZM)

  # Compute absolute salinity from practical S:
  print('Computing absolute S')
  kdm   = len(ZM)
  Z3d   = np.tile(ZM, idm*jdm).reshape((idm,jdm,kdm))
  Z3d   = np.transpose(Z3d, (2, 1, 0))
  PR    = np.zeros((kdm,jdm,idm))
  for kk in range(kdm):
    pr_db, _ = msw.sw_press(Z3d[kk,:,:].squeeze(), hlon)
    PR[kk,:] = pr_db

  SA = gsw.SA_from_SP(S3d, PR, hlon, hlat)

# Compute conservative T from potential T
  print('Computing conservative T')
  CT3d = gsw.CT_from_pt(SA, T3d) 

  # Derive bottom T:
  kdm, jdm, idm = T3d.shape
  Tbtm = np.zeros((jdm,idm))*np.nan
  dpmin = 1.e-1
  for ik in range(1,kdm):
    dpup  = dP[ik-1,:].squeeze()
    dpbtm = dP[ik,:].squeeze() 
    tz    = CT3d[ik-1,:]
    if ik < kdm-1:
      Jb, Ib = np.where( (dpup > dpmin) & (dpbtm <= dpmin) )
    else:
  # Deep layers include all left:
      Jb, Ib = np.where( dpup > dpmin )
    if len(Jb) == 0: continue
    Tbtm[Jb, Ib] = tz[Jb, Ib]

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

    CPA[icc, icl-1] = CP_area

  icc += 1

TM  = np.array(TM)
# Get rid off missed data
if len(TM) < nT:
  CPA = CPA[:icc,:]


f_save = True
if f_save:
  pthdump = pthseas["MOM6_NEP"][expt]["pthdump"] 
  fldump  = f'{runname}_coldpool_areaTclass.pkl'
  drfldump = os.path.join(pthdump, fldump)
  print(f'Saving ---> {drfldump}')

  with open(drfldump, 'wb') as fid:
    pickle.dump([TM, CPA, TCLASS],fid)


# Plot time series:
plt.ion()

Tdays = TM-TM[0]
cff   = 1.e-5
CPsc  = CPA*cff

# Show Month day1:
dv1   = mtime.datevec(TM[0])
dnmb1 = mtime.datenum([dv1[0],dv1[1],1])
dv2   = mtime.datevec(TM[-1] + 31)
dnmb2 = mtime.datenum([dv2[0],dv2[1],1])
TMp   = np.arange(dnmb1,dnmb2+1)
DVp   = mtime.datevec2D(TMp)
Ip = np.where(DVp[:,2] == 1)
Tticks  = TMp[Ip] - TM[0]
DVticks = DVp[Ip,:3].squeeze()
ntcks   = len(Tticks)
tck_lbls = []
for ik in range(ntcks):
  ff = f'{DVticks[ik,1]}/{DVticks[ik,2]}'
  tck_lbls.append(ff)
  
plt.ion()

CLRS = np.array([[0., 0.2, 0.9],
                [0., 0.8, 1],
                [0.7, 0., 1],
                [0.9, 0.4, 0],
                [0.5, 0.3, 0]])


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
LNS = []
for icl in range(nTC-1):
  clr0 = CLRS[icl, :]
  CPa  = CPsc[:, icl] 
  t1 = TCLASS[icl]
  t2 = TCLASS[icl+1]
  tline = f'{t1:.1f} < t < {t2:.1f}'
  ln1, = ax1.plot(Tdays, CPa, linewidth=2, color=clr0, label=tline)
  LNS.append(ln1)
#
# Total area for all classes:
CPt = np.sum(CPsc, axis=1)
clr0 = [0,0,0]
ln1, = ax1.plot(Tdays, CPt, linewidth=2, color=clr0, label='Total')
LNS.append(ln1)

ax1.set_xticks(Tticks)
ax1.grid(True)
ax1.set_xticklabels(tck_lbls)
ax1.set_ylabel(f'km2 x {cff}')
ax1.set_xlabel(f'Months starting {DVp[0,0]}/{DVp[0,1]}/{DVp[0,2]}')
ax1.set_title(f'{runname} cold pool area km2')

ax3 = plt.axes([0.1, 0.2, 0.6, 0.22])
lgd = plt.legend(handles=LNS, loc='upper right')
ax3.axis('off')

btx = 'cold_poolarea_Tclass.py'
bottom_text(btx, pos=[0.1, 0.1])


