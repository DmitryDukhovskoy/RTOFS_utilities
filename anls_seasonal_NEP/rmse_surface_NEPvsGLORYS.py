"""
  COmpute RMSE for surface 2D fields from NEP nudged simulation to compared with
  GLORYS monthly fields 

Data fields from Liz:
location for monthly GLORYS means for NEP region:
/archive/e1n/datasets/GLORYS/monthly_means/

location for padded monthly GLORYS means, concatenated by year used to generate NEP clim nudging files:
/archive/e1n/datasets/GLORYS/monthly_climatologies/

location for monthly GLORYS means, regridded to NEP for nudging as individual months:
/archive/e1n/mom6/NEP/sponge/monthly_sponge_files/

location for padded monthly GLORYS means, regridded to NEP for nudging and concatenated by year:
/archive/e1n/mom6/NEP/sponge/clims/

The last directory contains the files I used for nudging the solution to GLORYS. 

Daily GLORYS reanalysis for NEP domain prepared by Liz:
/archive/e1n/datasets/GLORYS/YYYY/nep_10

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
import datetime
from datetime import datetime


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

# Output saved at 3-month intervals
dnmb0 = mtime.datenum([2011,4,30])
expt = 'glorys_nudging'
expt_name = 'NEP_physics_202404_nudging-15d'  # GLORYS extracted for NEP domain
varnm = 'sos'  # tos = SST, sos = SSS, ssh
lr0  = 1  # ocean layers from 1, ..., 50

dnmb0 = mtime.datenum([2011,1,1])
dv0 = mtime.datevec(dnmb0)
YR0, MM0, DD0 = dv0[:3]

print(f'RMSE {expt_name} {varnm} {YR0}')

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

pthtopo    = gridfls['MOM6_NEP']['seasonal_fcst']['pthgrid']
fgrid      = gridfls['MOM6_NEP']['seasonal_fcst']['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"]["seasonal_fcst"]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

HH = dstopo_nep['depth'].data
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

match varnm:
  case "tos":
    varnm_glorys = 'thetao'
  case "sos":
    varnm_glorys = 'so'
  case "ssh":
    varnm_glorys = 'zos'

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')


ndays    = 365
TM_rmse  = np.zeros((ndays,3))
RMSE     = np.zeros((ndays))
MONTHS   = np.array([1,4,7,10])
mo_old   = 0
itime_glorys = 0
yr_old   = 0
zmin     = -500.
for iday in range(ndays):
  dnmb = dnmb0 + iday
# Find month when the simulation started:
  DV = mtime.datevec(dnmb)
  YR, MM, DD = DV[:3]
  dm = MONTHS - MM
  imm = max(np.where(dm <= 0)[0])
  mo_start = MONTHS[imm]
  yr_start = YR0
  pthdata = gridfls['MOM6_NEP'][expt]['pthoutp'].format(expt_name=expt_name, yr=yr_start, mo=mo_start)
  flnm = gridfls['MOM6_NEP'][expt]['fdata_day']
  dfl_nep = os.path.join(pthdata,flnm)
  print(f'Processing {YR}/{MM}/{DD}')

  if mo_start != mo_old: 
    print(f'Loading NEP {dfl_nep}')
    dset = xarray.open_dataset(dfl_nep)
    mo_old = mo_start

    # Create time array
    Time = dset['time'].data
    nrec = len(Time)
    TM = np.zeros((nrec))
    for irc in range(nrec):
      dmm = Time[irc].astype('datetime64[D]')
      dmm = dmm.astype(datetime)
      yr  = dmm.year
      month = dmm.month
      mday  = dmm.day
      TM[irc] = mtime.datenum([yr,month,mday])

  dtm = abs(TM-dnmb)
  
  indx0 = np.where(dtm <= 0.01)[0][0]
  A2d = dset[varnm].isel(time=indx0).data.squeeze()


  if varnm == 'zos':
    A2d = dset[varnm].isel(time=itime).data.squeeze()
    zmin = -200.
    dmm  = A2d.copy()
    dmm  = np.where(HH>zmin, np.nan, dmm)
    A2d = A2d - np.nanmean(dmm)


  # Get GLORYS monthly data for this time period
  # Note that monthly GLORYS fields are nominally referenced to mid-month
  # Select time: data padded, i.e. start from 0 and + 1 month at the end:
  #  time = -15.5, 15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 
  #    349.5, 380.5 ;

  if yr_old != YR:
    pthdata_glorys = gridfls['GLORYS_NEP']["monthly"]['pthoutp']
    flnm_glorys = gridfls['GLORYS_NEP']["monthly"]['fdata'].format(year=YR0)
    dfl_glorys = os.path.join(pthdata_glorys,flnm_glorys)
    print(f'Loading {dfl_glorys}')
    dset_glorys = xarray.open_dataset(dfl_glorys)
    yr_old = YR

  if itime_glorys != MM:
    itime_glorys = MM
    # Subtract mean:
    if varnm == 'ssh':
      G2d = dset_glorys['zos'].isel(time=itime_glorys).data.squeeze()
      zmin = -200.
      dmm  = G2d.copy()
      dmm  = np.where(HH>zmin, np.nan, dmm)
      G2d = G2d - np.nanmean(dmm)
    else:
      idepth = 0
      G2d = dset_glorys[varnm_glorys].isel(time=itime_glorys, depth=idepth).data.squeeze()

    G2d = np.where(HH >= zmin, np.nan, G2d)

  # Time stamps:
  TM_rmse[iday,0] = dnmb
  TM_rmse[iday,1] = TM[indx0]
  TM_rmse[iday,2] = itime_glorys

  # RMSE for deep regions only:
  A2d = np.where(HH >= zmin, np.nan, A2d)

  D2   = (G2d - A2d)**2
  NN   = len(np.where(~np.isnan(D2))[0])
  rmse = np.sqrt(np.nansum(D2)/NN)
  RMSE[iday] = rmse 

# rmin, rmax = mutob.minmax_clrmap(A2d, cpnt=0)
mo_old = 0
Mday = []
for iday in range(len(TM_rmse[:,0])):
  dv = mtime.datevec(TM_rmse[iday,0])
  mo = dv[1]
  if mo != mo_old:
    Mday.append(iday)
    mo_old = mo
Mday.append(Mday[-1]+31)

Tday = TM_rmse[:,0]-TM_rmse[0,0]
plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.3, 0.8, 0.6])
ax1.plot(Tday, RMSE)

sttl = f"RMSE: {expt_name} daily {varnm} vs GLORYS start:{YR0}/{MM0}/{DD0}"
ax1.set_title(sttl)
ax1.set_xticks(Mday)
ax1.grid('on')

sinfo=dfl_nep
ax3 = fig1.add_axes([0.1,0.05,0.8,0.02])
ax3.text(0,0, sinfo, fontsize=8)
ax3.axis('off')


btx = 'rmse_surface_NEPvsGLORYS.py'
bottom_text(btx, fsz=6, pos=[0.05, 0.03])










