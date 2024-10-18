"""
  GFDL simulations (Liz's forecast)

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
MOS    = 1
DDS    = 1    
nens   = 2
yrrun  = 1994  # year to plot
morun  = 5     # month to plot
mdrun  = 15 # month day

dnmbS = mtime.datenum([YRS,MOS,DDS])
dnmbR = mtime.datenum([yrrun, morun, mdrun])
# What month number:
nmonth = (yrrun-YRS)*12 + morun  # consecutive month #

print(f'Deriving data for {yrrun}/{morun} Record # = {nmonth}')

run_nm = 'MOM6_NEP_GFDL'
expt   = 'NEP_BGC_seas'

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

pthfcst    = pthseas[run_nm][expt]['pthoutp']
pthtopo    = pthseas[run_nm][expt]['pthgrid']
fgrid      = pthseas[run_nm][expt]['fgrid']
ftopo_mom  = pthseas[run_nm][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
ndav       = pthseas[run_nm][expt]['ndav']  # # of days output averaged

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

dfmom6 = os.path.join(pthfcst,'ocean_monthly_z.199301-201912.thetao.nc')
dset   = xarray.open_dataset(dfmom6)
T3d    = dset['thetao'].data[nmonth-1,:].squeeze()
ZM     = -dset['z_l'].data
ZZ     = -dset['z_i'].data

dfmom6 = os.path.join(pthfcst,'ocean_monthly_z.199301-201912.so.nc')
dset_so= xarray.open_dataset(dfmom6)
S3d    = dset_so['so'].data[nmonth-1,:].squeeze()


import mod_swstate as msw
import conversions as gsw
# Compute conservative T from potential T
#Tcons = gsw.CT_from_pt([20, 0, 35],[20, 0, 25]) 
# Derive potential enthalpy as Hpot = Tcons*Cp0 
# where Cp0 is the heat capacity of seawater to be 
# 53989.244 952 928 15 J kg-1 K-1
# potential enthalpy [J/kg] can be used as "heat content per unit mass"
# Compute absolute salinity from practical S:
print('Computing pressure')
jdm, idm = HH.shape
kdm   = len(ZM)
ZM3d   = np.tile(ZM, idm*jdm).reshape((idm,jdm,kdm))
ZM3d   = np.transpose(ZM3d, (2, 1, 0))
PR    = np.zeros((kdm,jdm,idm))
for kk in range(kdm):
  pr_db, _ = msw.sw_press(ZM3d[kk,:,:].squeeze(), hlon)
  PR[kk,:] = pr_db
 
print('Computing abs salinity')
SA = gsw.SA_from_SP(S3d, PR, hlon, hlat) 
# Compute conservative T
#Tcons = gsw.CT_from_pt(S3d, T3d)
# Potential enthalpy from abs. S and pot. T:
print('Computing potential enthalpy')
Hpot = gsw.pot_enthalpy_from_pt(SA, T3d)

# Density
Rho = msw.sw_dens0(SA,T3d)

# layer thickness:
# derive dP3d - layer thicknesses
print('Deriving layer thickness')
dP3d = manseas.derive_dP_from_ZZtopo(ZZ, HH)


# Heat content = depth integral (rho(z)*Hpot(z)*dz) = J/m2
# Integrate potential enthalpie over depth and divide by the total depth
# this gives J/m2 in 1 m of water (J/m3):
coeff = 1.e-6
HpotZ = np.nansum(Rho * Hpot * dP3d, axis=0)/np.sum(dP3d, axis=0)*coeff


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

sttl = f"{run_nm} Enth_pot*{coeff} J/m3, from monthly mean T/S {yrrun}/{morun}"
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
 

btx = 'heat_content_gfdl.py'
bottom_text(btx, pos=[0.2, 0.01])

