"""
  Calc and plot RMSE for monthly 3D fields (BGC restart fields)
  For hindcast, 3D bio fields have not been saved, thus use restart

  from BGC hindcasts 02 (corrupted ERA5) and 03 (corrected ERA5)
  BGC variables

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
parser.add_argument("--varnm", help="o2 no3 po4 sio4", type=str, required=True)
parser.add_argument("--zz", help="Depth to analyze, m, default - set of depths", type=float)
args = parser.parse_args()

varnm  = args.varnm if args.varnm else None
zz_plt = args.zz if args.zz else None
if zz_plt is not None:
  zz_plt = -abs(zz_plt)

YRS = 1993
YRE = 2000

# If zz is not specified, plot several depths:
if zz_plt is None:
  ZZP = [-10, -50, -100, -150, -500] 
else:
  ZZP = [zz_plt]

nzz = len(ZZP)

CLR = [[0., 0.4, 0.9],
       [0., 0.8, 0.5],
       [0.9, 0.6, 0],
       [0.8, 0., 0.5],
       [0.5, 0., 0.9],
       [0.5, 0.4, 0]]

hcast2_bgc = 'NEPbgc_nudged_hindcast02'
hcast3_bgc = 'NEPbgc_nudged_hindcast03'

hcst_time = 3 # f/csat time interval, months
hcst_interv = np.array([x for x in range(1,12+hcst_time,hcst_time)], dtype=int)
ilr0 = None

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


# Read vertical layers from oceanm archive file:
hcst2 = 'NEPbgc_nudged_hindcast02'
hcst3 = 'NEPbgc_nudged_hindcast03'
pthrst2 = f'/archive/Dmitry.Dukhovskoy/fre/NEP/hindcast_bgc/{hcst2}/restart'
pthrst3 = f'/archive/Dmitry.Dukhovskoy/fre/NEP/hindcast_bgc/{hcst3}/restart'

def read_restart_mom6(YRI, MMI, varnm, pthrst, iselZ=None):
  """
    3 month runs are assumed for the hindcasts
  """
  pthfull = os.path.join(pthrst,f'restdate_{YRI}{MMI:02d}01')

  match varnm:
    case 'o2' | 'po4' | 'sio4' | 'no3':
      varnc = varnm
      flnm = f'MOM_{YRI}{MMI:02d}01.res_2.nc'  # note res number may be different
    case _:
      varnc = varnm
      flnm = f'MOM_{YRI}{MMI:02d}01.res.nc'  # note res number may be different

  dflnm = os.path.join(pthfull,flnm)
  with xarray.open_dataset(dflnm, decode_times=False) as ds_mom:
    if (iselZ is not None):
      AA = ds_mom[varnc].isel(Layer=iselZ).data.squeeze()
    else:
      AA = ds_mom[varnc].data

  return AA

def find_depth_indx(zz_plt, pthrst, dnmb0):
  yr, mm = mtime.datevec(dnmb0)[:2]
  ZM = read_restart_mom6(yr,mm,'Layer',pthrst)
  ZM = -abs(ZM)
  dZ = np.abs(ZM-zz_plt)
  ilr0 = np.argmin(dZ)
  lr0  = ilr0+1
  zz0 = ZM[ilr0]  # actual depth to be plotted
 
  return ilr0, lr0, zz0

MCAL = []
YCAL = []
nmnth_file = 3   # assuming 3 months in monthly file
for yr in range(YRS,YRE+1):
  for mm in range(1,13,nmnth_file):
    if yr == 1993 and mm == 1:
      # Skip 1993/01 - no restart
      continue
    MCAL.append(mm)
    YCAL.append(yr)

len_cal = len(MCAL)  

itot = -1  # total rec counter
TM  = [] 
RMSE = np.zeros((len_cal,nzz))*np.nan
ZZ0 = []
for ii in range(len(MCAL)):
  MINIT = MCAL[ii]    # f/cast restart month
  YY = YCAL[ii] 
  dnmb0 = mtime.datenum([YY,MINIT,1])
  TM.append(dnmb0)

  print(f'Calculating RMSE from restart for {YY}/{MINIT}')
  dnmb_ref = mtime.datenum([1993,1,1])

  # Depths:
  for kk in range(nzz):
    zz_plt = ZZP[kk]
    ilr0, lr0, zz0 = find_depth_indx(zz_plt, pthrst2, dnmb0)

    # Read restart files from hindcast:
    # Convert to micro-mole/kg
    A2d_bgc2 = read_restart_mom6(YY,MINIT, varnm, pthrst2, iselZ=ilr0)*1.e6
    A2d_bgc3 = read_restart_mom6(YY,MINIT, varnm, pthrst3, iselZ=ilr0)*1.e6

    if ii == 0:
      ZZ0.append(zz0)

    sqerr = (A2d_bgc2-A2d_bgc3)**2
    nB = np.count_nonzero(~np.isnan(sqerr))
    rmseBP = np.sqrt(np.nansum(sqerr)/nB)
    RMSE[ii,kk] = rmseBP

    print(f"{zz0:.1f}m max RMSE: {rmseBP:.4f}")
    

TM = np.array(TM).flatten()
DV = mtime.datevec1D(TM, fHR=False)
DV = np.array(DV).transpose()
TMyr = TM*0.
for ik in range(len(TM)):
  dnmb = TM[ik]
  yr, jday = mtime.dnmb2jday(dnmb)
  if yr%4 == 0:
    ndays = 366
  else:
    ndays = 365
  TMyr[ik] = yr + np.floor(jday-1)/ndays


plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
# Time series RMSE
ax1 = plt.axes([0.1, 0.5, 0.85, 0.4])
hndls = []
for kk in range(nzz):
  rmse_z = RMSE[:,kk]
  clr = CLR[kk]
  zz0 = ZZ0[kk]
  ln1, = ax1.plot(TMyr,rmse_z, '-', linewidth=2, color=clr, label=f'z={zz0:.1f} m')
  hndls.append(ln1)

if nzz == 0:
  sttl = f'RMSE {varnm} (mmole/kg) z={zz0}m NEPbgc_hindcast02 vs NEPbgc_hindcast03 {YRS}-{YRE}'
else:
  sttl = f'RMSE {varnm} (mmole/kg) depths NEPbgc_hindcast02 vs NEPbgc_hindcast03 {YRS}-{YRE}'

ax1.set_title(sttl)
ax1.grid('on')

# Legend
if nzz > 0:
  ax2 = plt.axes([0.7, 0.25, 0.25, 0.18])
  ax2.legend(handles=hndls, loc='upper right')
  ax2.axis('off')

btx = 'calc_RMSE_hind2hind_3Dphys.py'
bottom_text(btx, pos=[0.05,0.35])

