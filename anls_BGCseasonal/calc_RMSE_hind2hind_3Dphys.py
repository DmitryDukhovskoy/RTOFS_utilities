"""
  Calc and plot RMSE for monthly 3D fields (ocean_month_z.nc)
  from BGC hindcasts 02 (corrupted ERA5) and 03 (corrected ERA5)
  for physical variables

  ssh, tos, sos, ssu, ssv

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
parser.add_argument("--varnm", help="field to analyze: temp, salin, uvel, vvel", type=str, required=True)
parser.add_argument("--zz", help="Depth to analyze, m, default - set of depths", type=float)
args = parser.parse_args()

varnm  = args.varnm if args.varnm else None
zz_plt = args.zz if args.zz else None
if zz_plt is not None:
  zz_plt = -abs(zz_plt)

# If zz is not specified, plot several depths:
if zz_plt is None:
  ZZP = [-10, -50, -100, -150, -500] 
else:
  ZZP = [zz_plt]

nzz = len(ZZP)

YRS = 1993
YRE = 2000

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


match varnm:
  case 'temp':
    varnc = 'thetao'
  case 'salin':
    varnc = 'so'
  case 'uvel':
    varnc = 'uo'
  case 'vvel':
    varnc = 'vo'


def find_depth_indx(ZM,zz_plt):
  ZM = ds_bgc2['z_l'].data
  ZM = -abs(ZM)
  dZ = np.abs(ZM-zz_plt)
  ilr0 = np.argmin(dZ)
  lr0  = ilr0+1
  zz0 = ZM[ilr0]

  return ilr0, lr0, zz0

MCAL = []
YCAL = []
nmnth_file = 3   # assuming 3 months in monthly file
for yr in range(YRS,YRE+1):
  for mm in range(1,13,nmnth_file):
    MCAL.append(mm)
    YCAL.append(yr)

len_cal = len(MCAL)  
len_tot = len_cal*nmnth_file  

itot = -1  # total rec counter
TM  = [] 
RMSE = np.zeros((len_tot,nzz))*np.nan
ZZ0 = []
for ii in range(len(MCAL)):
  MINIT = MCAL[ii]
  YY = YCAL[ii] 

  # Find init date for given month, assuming hcst_time (n months) f/cast interval
  #kint = np.searchsorted(hcst_interv, MM, side='right') - 1
  #assert(hcst_interv[kint] <= MM < hcst_interv[kint+1]), f'Wrong time bin {kint} for {MMA}'
  #MINIT = hcst_interv[kint]
  #imo = MM-MINIT      # current month in the archive output

  pthhcst2 = (
       pthseas['MOM6_NEP']['hindcast']['pthoutp'].format(
         hindcast_name=hcast2_bgc, YR=YY, MM=MINIT, DD=1
    )
  )

  pthhcst3 = '/archive/Dmitry.Dukhovskoy/fre/NEP/hindcast_bgc/NEPbgc_nudged_hindcast03/history/' +\
            f'{YY}{MINIT:02d}01/'


  print(f'Calculating RMSE for {YY}/{MINIT}')
  dnmb_ref = mtime.datenum([1993,1,1])

  # Depths:
  for kk in range(nzz):
    zz_plt = ZZP[kk]

    # Read hindcasts:
    docn2 = os.path.join(pthhcst2,f'ocean_month_z.nc')
    print(f'Reading {varnm} from {docn2}')
    with xarray.open_dataset(docn2, decode_times=False) as ds_bgc2:
      ZM = ds_bgc2['z_l'].data
      ilr0, lr0, zz0 = find_depth_indx(ZM, zz_plt)    
    
      A3d_bgc2 = ds_bgc2[varnc].isel(z_l=ilr0).data.squeeze()
      TIME2 = ds_bgc2['time'].data + dnmb_ref

    if ii == 0:
      ZZ0.append(zz0)

    docn3 = os.path.join(pthhcst3,f'ocean_month_z.nc')
    print(f'Reading {varnm} from {docn3}')
    with xarray.open_dataset(docn3, decode_times=False) as ds_bgc3:
      A3d_bgc3 = ds_bgc3[varnc].isel(z_l=ilr0).data.squeeze()
      TIME3 = ds_bgc3['time'].data + dnmb_ref

    # Time should be over the same time intervals:
    D = np.abs(TIME3-TIME2)
    assert np.max(D) < 1.e-19, f'ERROR: Time arrays in hindcast2 and hindcast3 are not the same' 
    TIME = TIME2.copy()

    inrec = itot
    icc = -1
    rmse_max = 0.
    for dnmb0 in TIME:
      icc += 1
      A2d_bgc2 = A3d_bgc2[icc,:,:].squeeze()
      A2d_bgc3 = A3d_bgc3[icc,:,:].squeeze()
      sqerr = (A2d_bgc2-A2d_bgc3)**2
      nB = np.count_nonzero(~np.isnan(sqerr))
      rmseBP = np.sqrt(np.nansum(sqerr)/nB)

      if kk == 0:
        TM.append(dnmb0)  
 
      inrec += 1 
      RMSE[inrec,kk] = rmseBP

      rmse_max = np.max([rmse_max,rmseBP])

    print(f"{zz0:.1f}m max RMSE: {rmse_max:.4f}")
    
  itot = inrec


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
  sttl = f'RMSE {varnm} z={zz0}m NEPbgc_hindcast02 vs NEPbgc_hindcast03 {YRS}-{YRE}'
else:
  sttl = f'RMSE {varnm} at depths  NEPbgc_hindcast02 vs NEPbgc_hindcast03 {YRS}-{YRE}'

ax1.set_title(sttl)
ax1.grid('on')

# Legend
if nzz > 0:
  ax2 = plt.axes([0.7, 0.25, 0.25, 0.18])
  ax2.legend(handles=hndls, loc='upper right')
  ax2.axis('off')

btx = 'calc_RMSE_hind2hind_3Dphys.py'
bottom_text(btx, pos=[0.05,0.35])

