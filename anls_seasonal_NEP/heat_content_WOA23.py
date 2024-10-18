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
import pdb
import netCDF4
from netCDF4 import Dataset as ncFile
import importlib
import yaml
from yaml import safe_load

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
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_anls_seas as manseas
import mod_plot_xsections as mxsct
import matplotlib as mtplt
import mod_regmom as mregmom

varnm   = 'temp'  # temp / salin
sctnm   = 'xsct_BerSea'  
YR      = 1993
MM      = 11  # for seasonal indicate month, for annual: MM = 13

grd=0.25
if grd==0.25:
  cgrd=4
woa='woa23'

seas, decade, yr1_dec, yr2_dec = manseas.season_decade_woa(YR,MM) 

woa_seas = {"13": "Jan-Mar",
            "14": "Apr-Jun",
            "15": "Jul-Spt",
            "16": "Oct-Dec",
            "0": "annual"}

urlBase = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/'
urlT    = f"{urlBase}temperature/netcdf/{decade}/0.25/"
urlS    = f"{urlBase}salinity/netcdf/{decade}/0.25/"
tfnm    = f"woa23_{decade}_t{seas:02d}_{cgrd:02d}.nc"
sfnm    = f"woa23_{decade}_s{seas:02d}_{cgrd:02d}.nc"

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)


def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  dmm = np.copy(dmm0)
  return dmm

# Get lon/lat
iz   = 0
furl = os.path.join(urlT,tfnm)
ZZ  = read_field(furl,'depth')
ZZ  = -abs(ZZ)
latW = read_field(furl,'lat')
lonW0 = read_field(furl,'lon')

# Subsample domain of interest:
lon1 = 150.
lon2 = 255.
lat1 = 10.
lat2 = 82.

# Read T/S:
furl = os.path.join(urlT,tfnm)
var_read = 't_an'
A3d = read_field(furl,var_read)
A3d = np.where(A3d > 1.e10, np.nan, A3d)
# reshaffle to have -180/180 lon inside the domain
A3d, lonW = mmisc.shuffle3D_lon180(A3d, lonW0)

ix1 = np.argmin(np.abs(lonW-lon1))
ix2 = np.argmin(np.abs(lonW-lon2))+1
jx1 = np.argmin(np.abs(latW-lat1))
jx2 = np.argmin(np.abs(latW-lat2))+1
T3d = A3d[:,jx1:jx2,ix1:ix2]

furl = os.path.join(urlS,sfnm)
var_read = 's_an'
A3d = read_field(furl,var_read)
A3d = np.where(A3d > 1.e10, np.nan, A3d)
A3d, _ = mmisc.shuffle3D_lon180(A3d, lonW0)
S3d = A3d[:,jx1:jx2,ix1:ix2]

lonW = lonW[ix1:ix2]
latW = latW[jx1:jx2]
jdm  = len(latW)
idm  = len(lonW)
LONW = np.zeros((jdm,idm))
LATW = np.zeros((jdm,idm))
for ii in range(idm):
  LATW[:,ii]=latW
for jj in range(jdm):
  LONW[jj,:]=lonW

# Convert in situ --> potential T at z=0:
print('Converting T in situ --> T potential')
kdm  = len(ZZ)
Z3d  = np.tile(ZZ, idm*jdm).reshape((idm,jdm,kdm))
Z3d  = np.transpose(Z3d, (2, 1, 0))
Tp3d = mregmom.insitu2pot_3D(T3d, S3d, Z3d, LATW)


# Compute conservative T from potential T
#Tcons = gsw.CT_from_pt([20, 0, 35],[20, 0, 25]) 
# Derive potential enthalpy as Hpot = Tcons*Cp0 
# where Cp0 is the heat capacity of seawater to be 
# 53989.244 952 928 15 J kg-1 K-1
# potential enthalpy [J/kg] can be used as "heat content per unit mass"
# Compute absolute salinity from practical S:
import mod_swstate as msw
import conversions as gsw
PR   = np.zeros((kdm,jdm,idm))
for kk in range(kdm):
  pr_db, _ = msw.sw_press(Z3d[kk,:,:].squeeze(), LONW)
  PR[kk,:] = pr_db

# Estimate layer thickness:
dP = manseas.derive_dP_WOA(ZZ, Tp3d) 

print('Computing absolute Salinity')
SA = gsw.SA_from_SP(S3d, PR, LONW, LATW)
# Compute conservative T
#Tcons = gsw.CT_from_pt(S3d, Tp3d)
# Potential enthalpy from abs. S and pot. T:
print('Computing potential enthalpy')
Hpot = gsw.pot_enthalpy_from_pt(SA, Tp3d)

# Density
Rho = msw.sw_dens0(SA,Tp3d)

# Heat content = depth integral (rho(z)*Hpot(z)*dz) = J/m2
# Integrate potential enthalpie over depth and divide by the total depth
# this gives J/m2 in 1 m of water (J/m3):
coeff = 1.e-6
HpotZ = np.nansum(Rho * Hpot * dP, axis=0)/np.sum(dP, axis=0)*coeff

#   
#  Plot
# 
import mod_colormaps as mclrmp
#rmin, rmax = mclrmp.minmax_clrmap(HpotZ)
rmin, rmax = mclrmp.minmax_clrmap(HpotZ, pmin=1, pmax=90)
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

# Indices:
#Xdist = np.arange(len(Hbtm))
# 

seas_nm = woa_seas[f"{seas}"]
sttl = f"WOA23 Pot.Enthalpy x {coeff} J/m3 ,  decade:{yr1_dec}-{yr2_dec} {seas_nm} {sctnm}"

ss1 = 'Hpot: Potential Enthalpy from Gibbs function, integrated rho*Hpot*dp over depth and divided by total depth\n'
ss2 = 'S --> Abs. Salinity, Hpot(0) => Conserv. T = 0, thus Hpot < 0 => energy deficit\n'
sinfo = ss1 + ss2

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3300*1.e3, resolution='l',\
            projection='stere', lat_ts=55, lat_0=62, lon_0=-175)
xR, yR = m(LONW, LATW)

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

img = m.pcolormesh(xR, yR, HpotZ, cmap=clrmp, vmin=rmin, vmax=rmax)
#img = m.pcolormesh(xR, yR, lgHpotZ, cmap=clrmp, vmin=rmin, vmax=rmax)
#m.contour(xR, yR, HH, [-1000], colors=[(0,0,0)], linestyles='solid')
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


btx = 'heat_content_WOA23.py'
bottom_text(btx, pos=[0.2, 0.01])



