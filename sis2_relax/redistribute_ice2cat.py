"""
  Redistribute  relax ice fields by thcikness categories
  ice categories are hrad-coded in SIS_state_initialization.F90
  these can be changed in SIS_override
  real :: hlim_dflt(8) = (/ 1.0e-10, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5 /) ! lower thickness limits 1...CatIce

  note:   nCat_dflt = 5 ; if (slab_ice) nCat_dflt = 1
  and SIS_input/ SIS_override: NCAT_ICE = 5 (if not then use default)
  so only 5 categories by default are used

"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
import pickle
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
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc
from mod_utils_fig import bottom_text
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

YR0 = 1993
MM0 = 5
ifld = 'ithkn'  # ithkn, iarea
file_type = 'monthly'  # monthly, daily, ... or clim
                       # for climatologies, do not need padded time - data will be recycled
                       # for monthly, daily, etc. need -dt and +dt at the beginn/end 

ICAT = np.array([1.0e-10, 0.1, 0.3, 0.7, 1.1])

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

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
# Hgrid lon. lat:
hlon, hlat  = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

pthsis  = gridfls['MOM6_NEP'][run_name]['pthsis']
pthdata = '/work/Dmitry.Dukhovskoy/data/PIOMAS_ice'
flthck = 'piomas20c.heff.1901.2010.v1.0.nc'
varthck = 'sit'
flconc  = 'piomas20c.area.1901.2010.v1.0.nc'
varconc = 'sic'

dflthkn = os.path.join(pthdata, flthck)
dflconc = os.path.join(pthdata, flconc)

ds_thkn = xarray.open_dataset(dflthkn)
LAT  = ds_thkn['Latitude'].data
LON  = ds_thkn['Longitude'].data

# Read saved relax. fields:
flout = f'PIOMAS_ithkn_iconc_{YR0}_{file_type}.nc'
diclim = os.path.join(pthsis, flout)
ds_rlx = xarray.open_dataset(diclim)
Time = ds_rlx['time'].data
TM = mmisc.convert_nptime_to_datenum(Time)
dnmb0 = mtime.datenum([YR0,MM0,15,12])
D = abs(TM-dnmb0)
itime = np.argmin(D)
dv0 = mtime.datevec(TM[itime])
assert dv0[0]==YR0, f'Requested YR={YR0}, year in rlx file={dv0[0]}'
assert dv0[1]==MM0, f'Requested month={MM0}, month in rlx file={dv0[1]}'

Hice = ds_rlx['ithkn'].isel(time=itime).data
Hice = np.where(HH>=0, np.nan, Hice)
Cice = ds_rlx['iarea'].isel(time=itime).data
Cice = np.where(HH>=0, np.nan, Cice)


#i0 = 261
#j0 = 728
i0 = 150
j0 = 634
hice = Hice[j0,i0]
cice = Cice[j0,i0]

# Test:
import random
ICAT0 = np.array([1.0e-10, 0.1, 0.3, 0.7, 1.1])
ncat = len(ICAT0)

nn=50
CI = np.zeros((nn))
HI = np.zeros((nn))
CC = np.zeros((nn,ncat))
HC = np.zeros((nn,ncat))
print(f"Calling redistribute_hice")
for ii in range(nn):
  cice = random.uniform(0.,1.)
  hice = random.uniform(0.,4.)
  print(f"ii={ii}, hice={hice:.4f} cice={cice:.4f}")
  hcat, ccat = msisrlx.redistribute_hice(hice, cice, ICAT=ICAT0, ck_min=1.e-2)
  CI[ii] = cice
  HI[ii] = hice
  CC[ii,:] = ccat
  HC[ii,:] = hcat 


ICATK = np.append(ICAT0,[3])

from matplotlib.patches import Polygon
plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])

ii = 10
ccat = CC[ii,:]
hcat = HC[ii,:]
hice = HI[ii]
cice = CI[ii]
#plt.bar(ICAT0, ccat, color=[0.8,0.9,1], width=0.2)
#ax1.plot(ICAT0, CC[ii,:],'-o')
dltE=0.01
clr=[0.5,0.8,1]
for kk in range(ncat):
  hmin = ICATK[kk]+dltE
  hmax = ICATK[kk+1]-dltE
  verts = [(hmin,0),(hmin,ccat[kk]),(hmax,ccat[kk]),(hmax,0)]
  poly  = Polygon(verts, facecolor=clr, edgecolor=clr, zorder=5)
  ax1.add_patch(poly)
#  ax1.plot([hmin,hmax],[ccat[kk],ccat[kk]],'-',linewidth=2, color=[0.,0.5,0.9])

ax1.set_xlim([0,1.5])

stl = f"hice={hice:.4f}, cice={cice:.4f}"
ax1.set_title(stl)
ax1.set_xticks(ICAT0)
ax1.set_xlabel('Ice Cat min Thicknesses')
ax1.set_ylabel('partial area')
ax1.grid('on')

btx = 'redistribute_ice2cat.py'
bottom_text(btx)



