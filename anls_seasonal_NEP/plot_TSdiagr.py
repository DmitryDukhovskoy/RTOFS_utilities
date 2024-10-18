"""
  Compute and plot 
  T/S diagram - monthly means
  TO define water masses
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib as mtplt
import xarray
from copy import copy
import matplotlib.colors as colors
from yaml import safe_load
import importlib

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

from mod_utils_fig import bottom_text
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_mom6 as mmom6
import mod_utils_ob as mutob
importlib.reload(mutob)
import mod_anls_seas as manseas
import mod_solver as msolv

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# Start of the run - needed only for seasonal forecasts:
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1    
nens   = 2

# Average over time period:
dnmbS = mtime.datenum([1993,4,1])
dnmbE = mtime.datenum([1993,6,30])
dstrS = mtime.datestr(dnmbS, show_hr=False)
dstrE = mtime.datestr(dnmbE, show_hr=False)

expt    = "seasonal_fcst"
runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
#expt    = 'NEP_BGCphys_GOFS'
#runname = 'NEP_physics_GOFS-IC'
#dnmbS   = mtime.datenum([YRS,MOS,DDS]) 
#dnmbR   = dnmbS + dayrun - 1

print(f'Expt: {expt} Run: {runname} T/S for {dstrS}-{dstrE}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

if expt == 'seasonal_fcst':
  pthfcst0  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(runname=runname)
else:
  dnmb0    = dnmbR
  dv0      = mtime.datevec(dnmb0)
  YR0, MM0, DD0 = dv0[:3]
  jday0    = int(mtime.date2jday([YR0,MM0,DD0]))
  pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YY=YR0, MM=MM0)
pthtopo    = gridfls['MOM6_NEP'][expt]['pthgrid']
fgrid      = gridfls['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"][expt]["ftopo"]
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

#ndays = mtime.month_days(momean, yrmean)
dnmbP = 0
dltm  = 5
icc   = 0
dnmbS = int(dnmbS)
dnmbE = int(dnmbE)
for dnmb0 in range(dnmbS, dnmbE+1, dltm):
# Find closest output:
#  dnmbR = mtime.datenum([yrmean, momean, mday])
  dnmbR = dnmb0
  dvR = mtime.datevec(dnmbR)

  print(f'Expt: {expt} Run: {runname} Plot date: {dvR[0]}/{dvR[1]}/{dvR[2]}')

  if expt == 'NEP_BGCphys_GOFS':
    ocnfld = 'ocean'
  else:
    ocnfld = 'oceanm'

  if expt == 'seasonal_fcst':
    pthfcst = os.path.join(pthfcst0,f'{ocnfld}_{dvR[0]}{dvR[1]:02d}')

  YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthfcst, dnmbR, fld=ocnfld)
  dv0  = mtime.datevec(dnmb0)
  YR0, MM0, DD0 = dv0[:3]

  flocn_name = pthseas['MOM6_NEP'][expt]['focname'].format(YR=YR0, jday=jday0)
  dfmom6 = os.path.join(pthfcst, flocn_name)

  dset   = xarray.open_dataset(dfmom6)

  ZM  = -dset['zl'].data
  T3d = dset['potT'].data[0,:].squeeze()
  S3d = dset['salt'].data[0,:].squeeze()
  dP  = dset['h'].data[0,:].squeeze()
  dP  = np.where(dP < 1.e-3, 0., dP)
  ZZ  = mmom6.zm2zz(ZM)

  import mod_swstate as msw
  import conversions as gsw
  # Compute conservative T from potential T
  #Tcons = gsw.CT_from_pt([20, 0, 35],[20, 0, 25]) 
  # Derive potential enthalpy as Hpot = Tcons*Cp0 
  # where Cp0 is the heat capacity of seawater to be 
  # 53989.244 952 928 15 J kg-1 K-1
  # potential enthalpy [J/kg] can be used as "heat content per unit mass"
  # Compute absolute salinity from practical S:
  jdm, idm = HH.shape
  kdm   = len(ZM)
  Z3d   = np.tile(ZM, idm*jdm).reshape((idm,jdm,kdm))
  Z3d   = np.transpose(Z3d, (2, 1, 0))
  PR    = np.zeros((kdm,jdm,idm))
  for kk in range(kdm):
    pr_db, _ = msw.sw_press(Z3d[kk,:,:].squeeze(), hlon)
    PR[kk,:] = pr_db
   
  SA = gsw.SA_from_SP(S3d, PR, hlon, hlat) 
  # Compute conservative T
  #Tcons = gsw.CT_from_pt(S3d, T3d)
  # Potential enthalpy from abs. S and pot. T:
  #Hpot = gsw.pot_enthalpy_from_pt(SA, T3d)

  # Discard regions outside Bering Sea:
  # Fill bottom with nans
  for ik in range(kdm):
    a2d = T3d[ik,:]
    a2d = np.where(MSKBS < 1, np.nan, a2d)
    a2d = np.where(HH > ZZ[ik+1], np.nan, a2d)
    T3d[ik,:] = a2d

    s2d = SA[ik,:]
    s2d = np.where(MSKBS < 1, np.nan, s2d)
    s2d = np.where(HH > ZZ[ik+1], np.nan, s2d)
    SA[ik,:] = s2d

    Zpnts = np.where(np.isnan(s2d), np.nan, ZZ[ik+1])/abs(HH)  # normalize by depth
    Z3d[ik,:] = Zpnts
  if icc == 0:
    T1d = T3d[(~np.isnan(T3d))]
    S1d = SA[(~np.isnan(SA))]
    Z1d = Z3d[(~np.isnan(Z3d))]
  else:
    T1d = T1d + T3d[(~np.isnan(T3d))]
    S1d = S1d + SA[(~np.isnan(SA))]
    Z1d = Z1d + Z3d[(~np.isnan(Z3d))] # hybrid coord change in time
  
  icc += 1

T1d = T1d / float(icc)
S1d = S1d / float(icc)
Z1d = Z1d / float(icc)

# Density
SS, TT  = np.meshgrid(np.arange(15,38,0.1), np.arange(-1.8,12,0.1))
sgm0 = msw.sw_dens0(SS, TT)-1000.

#
# Colormap
import mod_colormaps as mclrmps
clrmp = mclrmps.colormap_temp(clr_ramp=[0.4,0.,0.6])
rmin  = -1.
rmax  = 0.

# Plot T/S
plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.2, 0.3, 0.6, 0.6])
img = ax1.scatter(S1d, T1d, c=Z1d,  vmin=rmin, vmax=rmax, s=8, cmap=clrmp)

ax1.axis('scaled')
ax1.set_xlim([28, 35.3])
ax1.set_ylim([-2, 11])

CS = ax1.contour(SS, TT, sgm0, [x for x in range(0,32,1)], linestyles='solid', 
            colors=[(0.6,0.6,0.6)])
ax1.clabel(CS, inline=True, fontsize=11)

#ax1.set_xticks(Tticks)
#ax1.grid(True)
#ax1.set_xticklabels(tck_lbls)
ax1.set_ylabel(r'$\theta$$^oC$')
ax1.set_xlabel('Absolute S')
ax1.set_title(f'{runname} T/S avg: {dstrS} - {dstrE}')

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


btx = 'plot_TSdiagr.py'
bottom_text(btx, pos=[0.1, 0.1])

