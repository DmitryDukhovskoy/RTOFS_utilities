"""
  Plot number of observations in the grid cells for  
  Regional climatology NNEP WOA23 1/10 degree grid

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
import pickle

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
import mod_colormaps as mclrmps
import mod_swstate as msw
import conversions as gsw

YR  = 9999 # indicate any year in the decade need to plot, 9999 - all years
MM  = 11
Zplt = 999.  # depth to plot, make Zplt > 0 to show Nobs in the bottom layer
#woa='woa23'

if Zplt <= 0:
  print(f'The number of T obs in layer = {Zplt:.3f} m, MM={MM} \n')
else:
  print(f'The number of T obs in the bottom layer MM={MM}\n')

woa_seas = {"13": "Jan-Mar",
            "14": "Apr-Jun",
            "15": "Jul-Spt",
            "16": "Oct-Dec",
            "0": "annual"}

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

# Load Bering Shelf boundary, saved in plot_cold_pool_seas.py:
pthanls = pthseas['MOM6_NEP']['seasonal_daily']['pthanls'].format(expt_nmb=2)  
dfnm_rbnd = os.path.join(pthanls,f'BeringShelf_boundary_816x342.pkl')
if not os.path.isfile(dfnm_rbnd):
  raise Exception(f'Region boundary need to be saved in plot_cold_pool_seas.py {dfnm_rbnd}')

with open(dfnm_rbnd, 'rb') as fid:
  RBND = pickle.load(fid)

Xbnd = RBND[0]
Ybnd = RBND[1]

def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  dmm = np.copy(dmm0)
  return dmm

def get_lonlat_ZZ(urlT, tfnm, lon1=150., lon2=255., lat1=10., lat2=82.):
  # Get lon/lat for specified region
  furl = os.path.join(urlT,tfnm)
  ZZ  = read_field(furl,'depth')
  ZZ  = -abs(ZZ)
  latW = read_field(furl,'lat')
  lonW = read_field(furl,'lon')
#  lonW = mmisc.shuffle1D_lon180_to0360(lonW0)
  ix1 = np.argmin(np.abs(lonW-lon1))
  ix2 = np.argmin(np.abs(lonW-lon2))+1
  jx1 = np.argmin(np.abs(latW-lat1))
  jx2 = np.argmin(np.abs(latW-lat2))+1
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
  
  return [ix1,ix2,jx1,jx2], LONW, LATW, ZZ

def read_TS(urlT, tfnm, var_read):
  # Read T/S from WOA23:
  furl = os.path.join(urlT,tfnm)
#var_read = 't_an'
  A3d = read_field(furl,var_read)
  A3d = np.where(A3d > 1.e10, np.nan, A3d)
# reshaffle to have -180/180 lon inside the domain
#  lonW0 = read_field(furl,'lon')
#  A3d, _ = mmisc.shuffle3D_lon180(A3d, lonW0)

  return A3d


if YR <= 2050:
  seas, decade, yr1_dec, yr2_dec = manseas.season_decade_woa(YR, MM, month2season=False)
else:
  seas, decade, yr1_dec, yr2_dec = manseas.season_decade_woa(YR, MM, \
                                   month2season=False, decadal_clim=False)

urlBase = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa/REGCLIM/NNPv2/DATA/'
urlT    = f"{urlBase}temperature/netcdf/{decade}/0.10/"
urlS    = f"{urlBase}salinity/netcdf/{decade}/0.10/"
tfnm    = f"nnp_{decade}_t{seas:02d}_10.nc"
sfnm    = f"nnp_{decade}_s{seas:02d}_10.nc"

IJX, LONW, LATW, ZZ = get_lonlat_ZZ(urlT, tfnm)
ix1,ix2,jx1,jx2 = IJX
DX, DY = mmom6.dx_dy(LONW, LATW)
Acell  = DX*DY
jdm, idm = LONW.shape
kdm  = len(ZZ)
Z3d  = np.tile(ZZ, idm*jdm).reshape((idm,jdm,kdm))
Z3d  = np.transpose(Z3d, (2, 1, 0))

# The number of observations of sea_water_temperature in each 
# grid-square at each standard depth level.
A3d = read_TS(urlT, tfnm, 't_dd')  # temp. observations
N3d = A3d[:,jx1:jx2,ix1:ix2]
kdm, jdm, idm = N3d.shape

N3d = np.where(N3d < -100., np.nan, N3d)
import mod_swstate as msw
import conversions as gsw
PR   = np.zeros((kdm,jdm,idm))
for kk in range(kdm):
  pr_db, _ = msw.sw_press(Z3d[kk,:,:].squeeze(), LATW)
  PR[kk,:] = pr_db

# Estimate layer thickness:
dP = manseas.derive_dP_WOA(ZZ, N3d)


# Land sea mask for WOA
iz = 0
LMsk = N3d[iz,:,:].squeeze()
LMsk = np.where(np.isfinite(LMsk), 1, 0)

# Derive Regional mask for Bering Shelf:
Iregn, Jregn = mmisc.find_closest_indx(Xbnd, Ybnd, LONW, LATW)

X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MSKBS, _, _ = mmisc.inpolygon_v2(X, Y, Iregn,Jregn)  # 
Jout, Iout = np.where(MSKBS<1)

#seas, decade, yr1_dec, yr2_dec = manseas.season_decade_woa(YR, 1, month2season=False)


clrmp = mclrmps.colormap_discrete()
CLRS  = clrmp.colors
nclrs = CLRS.shape[0]
clrmp.set_bad(color=[0.3,0.3,0.3])
clrmp.set_under(color=[1,1,1])

cmp_Lmsk = mclrmps.colormap_landmask(clr0=[0.8,0.8,0.8])

rmin = 1
rmax = 11

# What depth level:
if Zplt <= 0:
  D = abs(ZZ-Zplt)
  iz0 = np.argmin(D)
  Nobs = N3d[iz0,:,:].squeeze()
  zdata = f'{Zplt:.1f}m'
else:
  Nobs = manseas.derive_bottom_temp(N3d, dP) 
  zdata = 'Bottom Lr'

Nobs = np.where(Nobs < 1, np.nan, Nobs)

plt.ion()

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3300*1.e3, resolution='l',\
            projection='stere', lat_ts=55, lat_0=62, lon_0=-175)

xR, yR = m(LONW, LATW)
xCP, yCP = m(Xbnd,Ybnd) # Cold Pool region


sttl = f"Reg. clim. 1/10 NNEPv2, nmb Tobs {zdata},  MM={MM} {yr1_dec}-{yr2_dec}"


fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

ax1.pcolormesh(xR, yR, LMsk, cmap=cmp_Lmsk, vmin=0, vmax=0.2)
#img = m.pcolormesh(xR, yR, Nobs, cmap=clrmp, vmin=rmin, vmax=rmax)
#ax1.axis('scaled')
img = ax1.scatter(xR, yR, c=Nobs, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(sttl)
ax1.plot(xCP,yCP,'-', color=[1, 0.6,0], linewidth=2)

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

sinfo = f'The number of T observations in each grid-square at each depth\n'
sinfo = sinfo + f'NOAA NCEI Regional Northern Northeast Pacific Clim. v2 1/10 deg\n'
sinfo = sinfo + f'{urlT}'
ax3 = fig1.add_axes([0.02, 0.04, 0.8, 0.05])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')

btx = 'nmbobsrv_NNEPWOA23.py'
bottom_text(btx, pos=[0.2, 0.01])

