"""
  Plot SST averaged over some area 
  for ensemble runs
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

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
varnm  = 'salin'  # temp (potential) / salin
sctnm  = 'xsct_EOB' 

# Start of the run 
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1    
nensR  = 1  #ens # for reference ensemble run
dnmbR  = mtime.datenum([1994,3,15])  # day/month to plot

expt    = "seasonal_fcst"
runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nensR:02d}'
dnmbS   = mtime.datenum([YRS,MOS,DDS]) 
dv_start = mtime.datevec(dnmbS)

dvR = mtime.datevec(dnmbR)
print(f'Expt: {expt} Run: {runname} Plot date: {dvR[0]}/{dvR[1]}/{dvR[2]}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

if expt == 'seasonal_fcst':
  pthfcst  = pthseas['MOM6_NEP'][expt]['pthwoutp'].format(runname=runname)
else:
  dnmb0    = dnmbR
  dv0      = mtime.datevec(dnmb0)
  YR0, MM0, DD0 = dv0[:3]
  jday0    = int(mtime.date2jday([YR0,MM0,DD0]))
  pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YY=YR0, MM=MM0)

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

# Select some region:
II = [188, 308]
JJ = [100, 320] 
lr = 5             # vertical layer to analyze

# Find closest output:
ocnfld = 'oceanm'

ENSR = [x for x in range(1,11)]
NensR = len(ENSR)
for iens in range(NensR):
  nens = ENSR[iens]
  runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
  pthfcst = pthseas['MOM6_NEP'][expt]['pthwoutp'].format(runname=runname)

  print(f"Processing ens={nens:02d} {pthfcst}")

  Fts, TM = manseas.timeser_spatavrg(pthfcst, YRS, MOS, JJ, II, lr, varnm, nens, ocnfld)

#  dF = F2d - F2dR
#  print(f"diff min/max: {np.nanmin(dF)}/{np.nanmax(dF)}")
  if iens == 0:
    dim1 = "Time"
    darr_var = xarray.DataArray(Fts, dims=(dim1), coords={dim1: TM})
    dset1D = xarray.Dataset({f"{varnm}_e{nens:02d}": darr_var}).copy()
  else:
    darr_var = xarray.DataArray(Fts, dims=(dim1), coords={dim1: TM})
    dset_var = xarray.Dataset({f"{varnm}_e{nens:02d}": darr_var})
    dset1D = xarray.merge([dset1D, dset_var])

# ===================
# Plotting
# ===================
btx = 'timeser_SST_ensmbl.py' 

plt.ion()

fgnmb = 1
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()

ax1 = plt.axes([0.1, 0.25, 0.8, 0.6])

for ie in range(NensR):
  nens = ENSR[ie]
  vards = f"{varnm}_e{nens:02d}"
  Ts = dset1D[vards].data.squeeze()
  TM = dset1D['Time'].data
  Tday = TM-TM[0]

  ax1.plot(Tday, Ts, '-')

ax1.grid('on')
ax1.set_xlabel('Forecast days')
ax1.set_ylabel('T')
if varnm == 'salin': ax1.set_ylabel('S')
dstart = f'{dv_start[0]}/{dv_start[1]}/{dv_start[2]}'
sttl = f'Seas f/cast init: {dstart}, T lr={lr} avrg j/i: {JJ[0]}:{JJ[1]}/{II[0]}:{II[1]} '
ax1.set_title(sttl)
bottom_text(btx)


