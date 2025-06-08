"""
  Calc and plot time series of RMSE for monthly sea ice conc/thickness 
  from test simulations against target relaxation fields from PIOMAS

  from test experiments with relaxation and / or target fields (monthly mean)
 
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
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

YR = 1993
MMI = 4    # month of initialization of the test runs
ndav = 5
varnm = 'iconc'
EXPTS = ['irlx2','irlx4','irlx5','irlx6','irlx7','irlx8', 'irlx9']

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="year to plot: 1993, ..., 2020", default=1993, type=int)
parser.add_argument("--varnm", help="field to plot: ithkn or iconc", type=str)
args = parser.parse_args()

if args.varnm:
  varnm = args.varnm
if args.yr:
  YR = args.yr

pthrlx = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax' 


outfld  = 'icem'
prfx = ''  # 19930401 - time stamp used in SIS2 output in file names, note that find
           # closest archive file does not work for 19930401.icem*.nc file names
           # rename files using ./rename_archive_v0.sh 0 in the output dir


def read_test_run_monthly(dnmb0, outfld, pthtest, prfx, varnm):
  """
    Some test have monthly fields saved, if not - derive from N-day average
  """
  YYR, MMR = mtime.datevec(dnmb0)[:2]

  # Check if monthly file exists:
  dfmnth = os.path.join(pthtest,'ice_month.nc')
  if os.path.isfile(dfmnth):
    import pandas as pd
    # Read monthly mean from saved ice_month.nc if exists
    dset = xarray.open_dataset(dfmnth)
    tm_nep = dset['time'].data
    tmP = pd.to_datetime(tm_nep)
    nrec = len(tmP)
    TNEP = np.zeros((nrec,3), dtype=int)
    for irec in range(nrec):
      yr0 = int(tmP.year[irec])
      mo0 = int(tmP.month[irec])
      dd0 = int(tmP.day[irec])
      TNEP[irec,:] = [yr0,mo0,dd0]

    yrfcst = TNEP[:,0]
    mmfcst = TNEP[:,1]
    if YYR < np.min(yrfcst) or YYR > np.max(yrfcst):
      raise Exception(f"{YYR} is out of range for saved ice_month: {np.min(yrfcst)}/{np.max(yrfcst)}")
    ifcst  = np.where((mmfcst==MMR) & (yrfcst==YYR))[0][0]

    HIce = dset['sithick'].isel(time=ifcst).data
    CIce = dset['siconc'].isel(time=ifcst).data
    if varnm == 'iconc':
      A2d = CIce
    elif varnm == 'ithkn':
      A2d = CIce*HIce
  else:
    print(f'ice_mmonth.nc not found in {pthtest}, deriving mean ...')

    A2d = msisrlx.calc_iconc_ithkn_mnthmean(dnmb0, pthtest, varnm)

  return A2d

def read_relax_piomas(dnmb0, pthsis, varnm):
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  flthck  = f'piomas_heff{YR0}_v21.nc'
  varthck = 'heff'
  flconc  = f'piomas_area{YR0}_v21.nc'
  varconc = 'area'

  # Read saved relax. fields:
  YR1 = YR0
  YR2 = YR0+1
  flout = f'PIOMASv21_ithkn_iconc_{YR1}_{YR2}_monthly.nc'
  diclim = os.path.join(pthsis, flout)
  print(f'Reading relax fields from {diclim}')
  ds_rlx = xarray.open_dataset(diclim)
  Time = ds_rlx['time'].data
  TM = mmisc.convert_nptime_to_datenum(Time)
  dnmb0 = mtime.datenum([YR0,MM0,15,12])
  D = abs(TM-dnmb0)
  itime = np.argmin(D)
  dv0 = mtime.datevec(TM[itime])
  assert dv0[0]==YR0, f'Requested YR={YR0}, year in rlx file={dv0[0]}'
  assert dv0[1]==MM0, f'Requested month={MM0}, month in rlx file={dv0[1]}'

  match varnm:
    case('ithkn'):
      ifld = 'ithkn'
    case('iconc'):
      ifld = 'iarea'

  A2dS = ds_rlx[ifld].isel(time=itime).data

  return A2dS


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

expt       = 'seasonal_daily'
pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

YCAL, MCAL = msisrlx.cal_months_forecast(MMI, YRI=1993)

len_cal = len(MCAL)  # should be 12 for most cases
nexpt = len(EXPTS)
RMSE  = np.zeros((len_cal,nexpt))
icc   = 0
for expt in EXPTS:
  pthtest = f'/work/Dmitry.Dukhovskoy/tmp/test_{expt}'   
  print(f'Calculating RMSE for test_{expt}')

  rmse_mo = np.zeros((len_cal))
  for ii in range(len(MCAL)):
    MM = MCAL[ii]
    YY = YCAL[ii] 
    dnmbR = mtime.datenum([YY,MM,15])
    Arlx = read_relax_piomas(dnmbR, pthrlx, varnm)
    Anep = read_test_run_monthly(dnmbR, outfld, pthtest, prfx, varnm)

    sqerr = (Anep-Arlx)**2
    Jice, Iice = np.where((Anep > 1.e-10) | (Arlx > 1.e-10)) 
    Nice = len(Jice)
    assert(Nice>0), f"No ice grid found at {YY}-{MM:02d}"
    rmse_mo[ii] = np.sqrt(1./float(Nice)*np.sum(sqerr[Jice,Iice]))

  RMSE[:,icc] = rmse_mo
  icc += 1

# -------------------
#
# Plot RMSE
#
# -------------------
tck_lbls = []
for ik in range(len_cal):
  yrf = YCAL[ik]
  mof = MCAL[ik]
  ff = f'{mof:02d}\n{yrf:02d}'
  tck_lbls.append(ff)

XT = np.arange(1,len_cal+1)
plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])
lgd = []
for ii in range(nexpt): 
  err = RMSE[:,ii]
  expt = EXPTS[ii]
  ln1, = ax1.plot(XT,err, '-', linewidth=2, label=expt)

  lgd.append(ln1)

ax1.set_xticks(XT)
ax1.grid('on')
ax1.set_xlabel('Months')
ax1.set_xticklabels(tck_lbls)

sttl = f'RMSE {varnm} for test expts with ice rlx. {YR}/{MMI}'
ax1.set_title(sttl)

ax2 = plt.axes([0.2, 0.1, 0.6, 0.2])
llg = plt.legend(handles=lgd, loc='upper right')
ax2.axis('off')


btx = 'calc_RMSE_icetests.py'
bottom_text(btx, pos=[0.2, 0.1])


