"""
  Plot cold pool:
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
f_insitu = True    # convert to in situ 

# Start of the run - needed only for seasonal forecasts:
# Start of the run - needed only for seasonal forecasts:
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1
nens   = 2    # ens # for ensemble runs
dnmbR  = mtime.datenum([1993,3,1])  # day to plot


#expt    = "seasonal_fcst"
#runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
expt    = 'NEP_BGCphys_GOFS'
runname = 'NEP_physics_GOFS-IC'
#expt    = 'NEP_seasfcst_LZRESCALE'
#runname = 'NEPphys_LZRESCALE_climOB_1993_04-e02'
dnmbS   = mtime.datenum([YRS,MOS,DDS]) 

dvR = mtime.datevec(dnmbR)
print(f'Expt: {expt} Run: {runname} Plot date: {dvR[0]}/{dvR[1]}/{dvR[2]}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

if expt == 'seasonal_fcst':
  pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(runname=runname)
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

# Find closest output:
if expt == 'NEP_BGCphys_GOFS':
  ocnfld = 'ocean'
else:
  ocnfld = 'oceanm'
if expt == 'seasonal_fcst':
  pthfcst = os.path.join(pthfcst,f'{ocnfld}_{dvR[0]}{dvR[1]:02d}')
  
YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthfcst, dnmbR, fld=ocnfld)
dv0  = mtime.datevec(dnmb0)
YR0, MM0, DD0 = dv0[:3]

flocn_name = pthseas['MOM6_NEP'][expt]['focname'].format(YR=YR0, jday=jday0)
dfmom6 = os.path.join(pthfcst, flocn_name)

# Averaging period:
dnmb_av1 = dnmb0 - np.floor(ndav/2)
dnmb_av2 = dnmb_av1 + ndav-1
dv_av1 = mtime.datevec(dnmb_av1)
yrs, mms, dds = dv_av1[:3]
dv_av2 = mtime.datevec(dnmb_av2)
yre, mme, dde = dv_av2[:3]


dset = xarray.open_dataset(dfmom6)
ZM   = -dset['zl'].data
T3d  = dset['potT'].data[0,:].squeeze()
S3d  = dset['salt'].data[0,:].squeeze()
dP   = dset['h'].data[0,:].squeeze()
dP   = np.where(dP < 1.e-3, 0., dP)
ZZ   = mmom6.zm2zz(ZM)

sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')
import mod_swstate as msw
import conversions as gsw
# Compute absolute salinity from practical S:
print('Computing absolute S')
jdm, idm = HH.shape
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


CLRS = [[0.6, 0.02, 0.6],
        [0.2, 0.38, 1],
        [1, 1, 1]]

import mod_colormaps as mclrmp
clrmp = mclrmp.colormap_posneg_uneven(CLRS)
clrmp.set_bad(color=[0.6,0.6,0.6])

rmin = -2.
rmax = 2.

sttl = f"{runname} bottom Tconserv: {yrs}/{mms}/{dds}-{yre}/{mme}/{dde}"
sinfo = 'Conservative T in the near-bottom layer'

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3300*1.e3, resolution='l',\
            projection='stere', lat_ts=55, lat_0=62, lon_0=-175)

xR, yR = m(hlon, hlat)

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

img = m.pcolormesh(xR, yR, Tbtm, cmap=clrmp, vmin=rmin, vmax=rmax)
m.contour(xR, yR, HH, [-1000], colors=[(0,0,0)], linestyles='solid')
m.contour(xR, yR, Tbtm, [2.], colors=[(0,0.5,1.)], linestyles='solid')
#ax1.axis('scaled')
ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.02, 0.8, 0.05])
ax3.text(0, 0, sinfo, fontsize=10)
ax3.axis('off')

btx = 'plot_coldpool.py'
bottom_text(btx, pos=[0.2, 0.01])



