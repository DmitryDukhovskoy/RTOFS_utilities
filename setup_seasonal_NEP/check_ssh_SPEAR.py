"""
  Check daily or monthly SSH OB created from SPEAR runs
"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
#import pickle
import matplotlib.pyplot as plt
from yaml import safe_load
#import scipy.fftpack as fftpack
import scipy.fft as sfft
import datetime
from datetime import datetime

import mod_utils_ob as mutob
importlib.reload(mutob)


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
sys.path.append('./seasonal-workflow')
from boundary import Segment
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
from mod_utils_fig import bottom_text

# Climatology derived for these years, started at mstart
# Inidicate start of the SPEAR forecast:
ens_spear  = 1      # ens run used to create OB's
yr_start   = 1993
mo_start   = 4
nsgm = 3  # OB segment: 3 - South OB, 4 - West OB
dnmb_start = mtime.datenum([yr_start,mo_start,1])
dv_start   = mtime.datevec(dnmb_start)
varnm  = 'ssh'

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)

seas_yaml = 'paths_seasfcst.yaml'
with open(seas_yaml) as ff:
  fseas = safe_load(ff)


# MOM6 NEP topo/grid:
run_name   = 'seasonal_fcst_daily'
pthtopo    = gridfls['MOM6_NEP'][run_name]['pthgrid']
fgrid      = gridfls['MOM6_NEP'][run_name]['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"][run_name]["ftopo"]
outdir     = gridfls['MOM6_NEP'][run_name]['pthoutp']
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
dftopo_mom  = os.path.join(pthtopo, ftopo_mom)
LONM, LATM  = mmom6.read_mom6grid(dfgrid_mom)
HHM         = mmom6.read_mom6depth(dftopo_mom)

segments = [ Segment(1, 'north', hgrid, output_dir=outdir),
             Segment(2, 'east',  hgrid, output_dir=outdir),
             Segment(3, 'south', hgrid, output_dir=outdir),
             Segment(4, 'west',  hgrid, output_dir=outdir)]

nOB = len(segments)

# Load mapping indices exist, gmapi:
dirgmapi = config['filesystem']['spear_mom_gmapi']
flgmaph  = f'spear2mom_NEP_OB_gmapi_hpnt.nc'
flgmapu  = f'spear2mom_NEP_OB_gmapi_upnt.nc'
flgmapv  = f'spear2mom_NEP_OB_gmapi_vpnt.nc'
dflgmaph = os.path.join(dirgmapi, flgmaph)
dflgmapu = os.path.join(dirgmapi, flgmapu)
dflgmapv = os.path.join(dirgmapi, flgmapv)
# h-point indices
dsh = xarray.open_dataset(dflgmaph)

spear_dir = config['filesystem']['nep_spear_subset'].\
                   format(year=dv_start[0], ens=ens_spear)

# Ssh - 1D sections
# Load ssh daily fields for NEP subset SPEAR 
# daily ssh available in ice_daily (with "SSH" variable) and ocean_daily ("ssh") for some years
# check if varnm is ssh or SSH in the data array
flnm_spear = f'NEP_spear_{dv_start[0]}{dv_start[1]:02d}.ssh_daily.nc'
ds = xarray.open_dataset(os.path.join(spear_dir,flnm_spear))
dsvars  = list(ds.keys())
for ill in range(len(dsvars)):
  if dsvars[ill] == 'ssh':
    varnm = 'ssh'
    break
  elif dsvars[ill] == 'SSH':
    varnm = 'SSH'
    break

ds_spear = mutob.read_spear_output(spear_dir, varnm, flnm_spear, fzint=True)

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()

# Average SSH across OB section and plot time series
isgm = nsgm-1

print(f'Processing ssh OB segment={nsgm}')
INDX   = dsh[f'indx_segm{nsgm:03d}'].data
JNDX   = dsh[f'jndx_segm{nsgm:03d}'].data
dset_segm = mutob.segm_topo(nsgm, HHM, hgrid)
distOB   = dset_segm['dist_supergrid'].data
Xbtm     = dset_segm['dist_grid'].data
Hbtm     = dset_segm['topo_segm'].data
Hbtm     = np.where(Hbtm > 0, 0., Hbtm)
segm_nm  = dset_segm['segm_name'].data[0]

# Segment coordinates SPEAR
dset_sgmspear = segments[isgm]
xOB_spear = dset_sgmspear.coords.lon.data
yOB_spear = dset_sgmspear.coords.lat.data
npnts     = len(xOB_spear)
nx_spear  = dset_sgmspear.nx
ny_spear  = dset_sgmspear.ny

# Calculate distance along the OB segment:
#  distOB_spear, _ = mutob.calculate_dist_section(xOB_spear, yOB_spear)
TM = ds_spear['time'].data
ndays = len(TM)

# Average over smaller segment: close to the region where SSH was analyzed
match nsgm:
  case 3:
    regn = 'poly_south'

II = fseas['ANLS_NEP'][regn]['II']
#JJ = fseas['ANLS_NEP'][regn]['JJ']
ii1 = np.min(II)
ii2 = np.max(II)

SSH = []
for itime in range(ndays):
  ssh_spear = ds_spear[varnm].isel(time=itime).data
# 4 vertices chosen for SSH interpolation onto MOM grid:
# Pick one
  sshOB_spear = ssh_spear[JNDX[:,0],INDX[:,0]].squeeze()
#  ssh_mn = np.nanmean(sshOB_spear)
  ssh_mn = np.nanmean(sshOB_spear[ii1:ii2])
#  ssh_mn = np.nanmean(sshOB_spear[310:320])
  SSH.append(ssh_mn)
  
# COnvert into datetime object:
nrec = len(TM)
DTM = []
for ii in range(nrec):
  dd   = datetime.strptime(str(TM[ii]),'%Y-%m-%d %H:%M:%S')
  yr   = dd.year
  mo   = dd.month
  mday = dd.day 
  dnmb = mtime.datenum([yr,mo,mday])
  DTM.append(dnmb)

DTM = np.array(DTM)
DAYS = DTM-DTM[0]
dT  = DAYS[1] - DAYS[0]

SSH = np.array(SSH)
#Qfr = fftpack.fft(SSH)  # Fourier transfer
#Frw = fftpack.fftfreq(len(DAYS), dT)

# Detrend:
Pcoef = np.polyfit(DAYS,SSH,4)
Plnm  = np.poly1d(Pcoef)
Pfit = Plnm(DAYS)

SSH_dtr = SSH-Pfit

N = len(SSH_dtr)
Qfr = sfft.rfft(SSH_dtr) / N
Frw = sfft.rfftfreq(n=N, d=1./365.) # unit = 1/12 of sampling period
Qfr = np.abs(Qfr)
Qfr[0] = np.nan

# Bin averaging
# Skip the 1st bin which is nan - low-freq removed
Nav = 2
Iav = [x for x in range(1,len(Frw),Nav)]
Iav.append(N)
Frq_avg = []
Qfr_avg = []
for ii in range(len(Iav)-1):
  i1 = Iav[ii]
  i2 = Iav[ii+1]-1
  if i2 < i1: i2=i1
  Frq_avg.append(np.mean(Frw[i1:i2+1]))
  Qfr_avg.append(np.mean(Qfr[i1:i2+1]))

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))

fig1.clf()
ax1  = plt.axes([0.1, 0.54, 0.85, 0.4])
ax1.plot(SSH,'-')
ax1.plot(Pfit)
sttl = f'SSH SPEAR avg over OB={segm_nm}, {flnm_spear}'
ax1.grid('on')
ax1.set_title(sttl)
ax1.set_xlabel('Time, days')

ax2  = plt.axes([0.08, 0.08, 0.85, 0.4])
#ax2.plot(Frw, np.abs(Qfr))
ax2.plot(Frq_avg, Qfr_avg)
ax2.set_yscale('log')
#ax2.set_xscale('log')

sttl2 = 'Spectrum, m2/day2'
ticks = ax2.get_xticks()
ax2.set_xticklabels([f'{tick/N:6.2f}' if tick!=0 else '$\infty$' for tick in ticks])
ax2.set_xlabel('cyc/day')

btx = 'check_ssh_SPEAR.py'
bottom_text(btx, pos=[0.01, 0.01])



