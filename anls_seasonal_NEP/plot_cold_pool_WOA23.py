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

YR      = 2008 # indicate any year in the decade need to plot
MM      = 11  # for seasonal indicate any month in the season, for annual: MM = 13

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

# Regional high-resolution clim for North Northeast pacific:
#https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa/REGCLIM/NNPv2/DATA/temperature/netcdf/B5C2/0.10/nnp_B5C2_t13_10.nc.html
urlBase = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/'
urlT    = f"{urlBase}temperature/netcdf/{decade}/0.25/"
urlS    = f"{urlBase}salinity/netcdf/{decade}/0.25/"
tfnm    = f"woa23_{decade}_t{seas:02d}_{cgrd:02d}.nc"
sfnm    = f"woa23_{decade}_s{seas:02d}_{cgrd:02d}.nc"

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

# Land sea mask for WOA
LMsk = T3d[iz,:,:].squeeze()
LMsk = np.where(np.isfinite(LMsk), -10, 1)

# Derive Regional mask for Bering Shelf:
Iregn, Jregn = mmisc.find_closest_indx(Xbnd, Ybnd, LONW, LATW)

X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MSKBS, _, _ = mmisc.inpolygon_v2(X, Y, Iregn,Jregn)  # 
Jout, Iout = np.where(MSKBS<1)

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
  pr_db, _ = msw.sw_press(Z3d[kk,:,:].squeeze(), LATW)
  PR[kk,:] = pr_db

# Estimate layer thickness:
dP = manseas.derive_dP_WOA(ZZ, Tp3d) 

print('Computing absolute Salinity')
SA = gsw.SA_from_SP(S3d, PR, LONW, LATW)
# Compute conservative T
CT3d = gsw.CT_from_pt(S3d, Tp3d)

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

# Mask outside region:
Tbtm = np.where( (MSKBS<1.) & (LMsk<0) , 1.e3, Tbtm)
# Check, should be empty:
j0,i0 = np.where( (np.isnan(Tbtm)) & (LMsk<0) )
if len(j0) > 0:
  print(f'WARNING: {len(j0)} points Bottom T is missing')

CLRS = [[0.6, 0.02, 0.6],
        [0.2, 0.38, 1],
        [0., 0.8, 0.5],
        [0.9, 0.6, 0],
        [1, 1, 1]]

clrmp = mclrmps.colormap_posneg_uneven(CLRS)
clrmp.set_bad(color=[0.6,0.6,0.6])

rmin = -1.8
rmax = 2.
tscntrs = [x/10 for x in range(-20,80,5)]
tslabels = [x/10 for x in range(-20,80,5)]

sinfo = 'Conservative T in the near-bottom layer'

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
width  = 2200*1.e3
height = 2200*1.e3
lat0   = 62.
lon0   = -172.

m = Basemap(width=width, height=height, resolution='l',\
            projection='stere', lat_ts=55, lat_0=lat0, lon_0=lon0)

xR, yR = m(LONW,LATW)

btx = 'plot_cold_pool_WOA23.py'
seas_nm = woa_seas[f"{seas}"]
sttl = f"WOA23 Bottom conserv T,  decade:{yr1_dec}-{yr2_dec} {seas_nm}"

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()

ax1 = manseas.plot_stereogr_axis(fig1, m, xR, yR, Tbtm, clrmp, rmin, rmax, \
                       btx=btx, tscntrs=tscntrs, tslabels=tslabels, sttl=sttl)

plt.sca(ax1)
ax1.contour(xR,yR, Tbtm, [2], linestyles='solid', colors=[(1., 0.4, 0.9)])



