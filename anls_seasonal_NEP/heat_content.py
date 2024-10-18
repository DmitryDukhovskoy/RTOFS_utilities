"""
  Compute and plot 
  ocean heat conent using potential enthalpy 
  see: MCDOUGALL, T.J, "Potential Enthalpy: A Conservative Oceanic Variable 
  for Evaluating Heat Content and Heat Fluxes", JPO, 2003

  and Rainer Feistel, Eberhard Hagen "On the GIBBS thermodynamic potential 
  of seawater", Porgr. Oceanogr. 36(4), 1995

  Potential enthalpy is evaluated from absolute salinity and potential T
  conservative Temperatuer = Potential Enthalpy / Cp0, 
  Cp0 (heat capacity of sea water, = 3989.244 952 928 15 J / kg * K)

  For plotting: depth-integrate pot. enthalpy and / total depth

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
dnmbR  = mtime.datenum([1993,11,30])  # day to plot

#dnmbS   = mtime.datenum([YRS,MOS,DDS])
#dayrun  = dnmbR - dnmbS + 1 # day to plot:
dvR     = mtime.datevec(dnmbR)

expt    = "seasonal_fcst"
runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
#expt    = 'NEP_BGCphys_GOFS'
#runname = 'NEP_physics_GOFS-IC'
#dnmbS   = mtime.datenum([YRS,MOS,DDS]) 
#dnmbR   = dnmbS + dayrun - 1

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

# Find closest output:
if expt == 'seasonal_fcst':
  pth_day1=os.path.join(pthfcst,f'oceanm_{YRS}{MOS:02d}') # dir with 1st output
  YR0, MM0, DD0, jday0 = mutob.find_NEPoutput_day(pth_day1, dnmbR, \
                                                  ndav, fld='oceanm')
  dnmb0   = mtime.jday2dnmb(YR0,jday0)
  pthfcst = os.path.join(pthfcst,f'oceanm_{YR0}{MM0:02d}')

flocn_name = pthseas['MOM6_NEP'][expt]['focname'].format(YR=YR0, jday=jday0)
dfmom6 = os.path.join(pthfcst, flocn_name)

# Averaging period:
dnmb_av1 = dnmb0 - np.floor(ndav/2)
#if dnmb_av1 < dnmbS: dnmb_av1=dnmbS
dnmb_av2 = dnmb_av1 + ndav-1

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
Hpot = gsw.pot_enthalpy_from_pt(SA, T3d)

# Density
Rho = msw.sw_dens0(SA,T3d)

# Heat content = depth integral (rho(z)*Hpot(z)*dz) = J/m2
# Integrate potential enthalpie over depth and divide by the total depth
# this gives J/m2 in 1 m of water (J/m3):
coeff = 1.e-6
HpotZ = np.nansum(Rho * Hpot * dP, axis=0)/np.sum(dP, axis=0)*coeff

# Log scale:
#if np.nanmin(HpotZ) <= 0:
#  lgHpotZ = np.log10(HpotZ + abs(np.nanmin(HpotZ))+1.e-6)

#   
#  Plot
# 
import mod_colormaps as mclrmp
#rmin, rmax = mclrmp.minmax_clrmap(HpotZ)
rmin, rmax = mclrmp.minmax_clrmap(HpotZ, pmin=1, pmax=80)
rmin = -9.
rmax = 3*abs(rmin)
#rmin = 3.2
#rmax = 4.2
#cmp_cold = mclrmp.colormap_cold()
#cmp_warm = mclrmp.colormap_warm()
#clrmp = mclrmp.colormap_ssh(cpos=cmp_warm, cneg=cmp_cold, nclrs=200)
#clrmp = mtplt.colormaps.get_cmap('turbo')
# Create pos/neg colorbar with uneven pos/neg colors, 0 = white
clrmp = mclrmp.colormap_posneg_uneven([])
clrmp.set_bad(color=[0.6,0.6,0.6])

dv_av1 = mtime.datevec(dnmb_av1)
yrs, mms, dds = dv_av1[:3]
dv_av2 = mtime.datevec(dnmb_av2)
yre, mme, dde = dv_av2[:3]

sttl = f"{runname} Enth_pot*{coeff} J/m3, avrg: {yrs}/{mms}/{dds} - {yre}/{mme}/{dde}"
ss1 = 'Hpot: Potential Enthalpy from Gibbs function, integrated rho*Hpot*dp over depth and divided by total depth\n'
ss2 = 'S --> Abs. Salinity, Hpot(0) => Conserv. T = 0, thus Hpot < 0 => energy deficit\n'
sinfo = ss1 + ss2 

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

img = m.pcolormesh(xR, yR, HpotZ, cmap=clrmp, vmin=rmin, vmax=rmax)
#img = m.pcolormesh(xR, yR, lgHpotZ, cmap=clrmp, vmin=rmin, vmax=rmax)
m.contour(xR, yR, HH, [-1000], colors=[(0,0,0)], linestyles='solid')
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
 

btx = 'heat_content.py'
bottom_text(btx, pos=[0.2, 0.01])

