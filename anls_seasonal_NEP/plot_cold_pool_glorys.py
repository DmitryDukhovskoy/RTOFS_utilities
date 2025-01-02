"""
  Monthly average bottom T for NEP region 
  GLORYS reanalysis 
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import time
import timeit
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import xarray
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
from yaml import safe_load

#PPTHN = '/home/Dmitry.Dukhovskoy/python'
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

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_misc1 as mmisc
import mod_colormaps as mclrmps

YAVRG = [x for x in range(2005,2015)]
#MAVRG = [1,2,3]  # months to average: Winter  JFM, Summer: JAS
#MAVRG = [4,5,6]  # months to average: Spring, AMJ
MAVRG = [7,8,9]  # months to average: Winter  JFM, Summer: JAS
#MAVRG = [10,11,12]  # months to average: Fall


pthoutp   = '/work/Dmitry.Dukhovskoy/data/glorys_Tbtm_NEP/'
pthglorys = '/archive/e1n/datasets/GLORYS/' 
pthtopo   = '/work/Dmitry.Dukhovskoy/data/glorys_topo_NEP/'

def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  dmm = np.copy(dmm0)
  return dmm

def lookup_ncvar(nc):
  ii=0
  for var in nc.variables.values():
    ii+=1
    print('--------\n')
    print('Var # {0}'.format(ii))
    print(var)

# Get GLORYS topo:
fltopo = 'GLORYS12_topoNEP_865x1321.pkl'
dfltopo = os.path.join(pthtopo, fltopo)
with open(dfltopo, 'rb') as fid:
  HH = pickle.load(fid)

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

# Region boundaries:
# Load Bering Shelf boundary, saved in plot_cold_pool_seas.py:
pthanls = pthseas['MOM6_NEP']['seasonal_daily']['pthanls'].format(expt_nmb=2)
dfnm_vrtx = os.path.join(pthanls,f'BeringShelf_region_vrtx.pkl')
if not os.path.isfile(dfnm_vrtx):
  raise Exception(f'Region boundary need to be saved in plot_cold_pool_seas.py {dfnm_vrtx}')

with open(dfnm_vrtx, 'rb') as fid:
  XYBND = pickle.load(fid)

Xbnd = XYBND[0]
Ybnd = XYBND[1]
#Xbnd = np.where(Xbnd>180, Xbnd-360., Xbnd)

# Get lon/lat region boundaries:
pthinput = os.path.join(pthglorys,'2015','nep_10')
flnm = f'GLORYS_REANALYSIS_NEP_2015-12-01.nc'
dflnm = os.path.join(pthinput,flnm)
ds_glorys = xarray.open_dataset(dflnm)

lon = ds_glorys['longitude'].data
lat = ds_glorys['latitude'].data
ssh  = ds_glorys['zos'].data.squeeze()
Glon, Glat = np.meshgrid(lon ,lat)
jdm, idm = Glon.shape[:2]

ZM  = -ds_glorys['depth'].data
ZZ  = mmom6.zm2zz(ZM)

IIG = []
JJG = []
for ii in range(len(Xbnd)):
  x0, y0 = Xbnd[ii], Ybnd[ii]
  i0, j0 = mutil.find_indx_lonlat(x0, y0, Glon, Glat)
  IIG.append(i0)
  JJG.append(j0)

# Derive regional mask for Bering Shelf:
X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MS, _, _ = mmisc.inpolygon_v2(X, Y, IIG, JJG)  # 
JBS, IBS = np.where( (MS == 1) & (HH >= -250) & (HH < 0) ) #exclude deeep regions
MSKBS  = np.zeros((jdm,idm))
MSKBS[JBS,IBS] = 1

# Average by months & years:
icc = 0
for YR in YAVRG:
  for MM in MAVRG:
    pthglmn = '/archive/e1n/datasets/GLORYS/monthly_means/'
    flglmn = f'GLORYS_REANALYSIS_NEP_{YR}-{MM:02d}.tmp.nc' 
    dflglorys = os.path.join(pthglmn,flglmn)
    print(f'Loading {dflglorys}')
    dset = xarray.open_dataset(dflglorys)
    TT = dset['thetao'].data.squeeze()  # potential T
    SS = dset['so'].data.squeeze()
    if icc == 0:
      T3d = TT.copy()
      S3d = SS.copy()
    else:
      T3d = T3d + TT
      S3d = S3d + SS

    icc += 1

T3d = T3d/icc
S3d = S3d/icc

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
  pr_db, _ = msw.sw_press(Z3d[kk,:,:].squeeze(), Glon)
  PR[kk,:] = pr_db

# Estimate layer thickness:
dz = abs(np.diff(ZZ))
dP = np.tile(dz, idm*jdm).reshape((idm,jdm,kdm))
dP = np.transpose(dP, (2,1,0))
Jz, Iz = np.where(HH >= 0)
dP[:,Jz,Iz] = 0.

for kk in range(kdm):
  Z2d = Z3d[kk,:].squeeze()
  zz0 = ZZ[kk+1]
  dPz = dP[kk,:].squeeze()
  Jz,Iz = np.where( (HH >= zz0) & (dPz > 1.e-6) )  # hit bottom
  if len(Jz) == 0: 
    continue
# correct dp to make lower interf ZZ[+1] at the bottom
  dzBtm = abs(HH-zz0)
  dP[kk,Jz,Iz] = dP[kk,Jz,Iz] - dzBtm[Jz,Iz] 
  dP[kk+1:,Jz,Iz] = 0.


# TO speed up computation of cons. T, blank southern part of the region:
#S3d[:,:400,:] = np.nan
#T3d[:,:400,:] = np.nan

SA = gsw.SA_from_SP(S3d, PR, Glon, Glat)

# Compute conservative T from potential T
print('Computing conservative T')
CT3d = gsw.CT_from_pt(SA, T3d)

# Derive bottom T:
#kdm, jdm, idm = T3d.shape
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
Tbtm = np.where( (MSKBS==0) & (HH<0) , 1.e3, Tbtm)
# Check, should be empty:
j0,i0 = np.where( (np.isnan(Tbtm)) & (HH < -10) )
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

sinfo = 'GLORYS12v1 Conservative T in the near-bottom layer'

run_info = f'GLORYS12V1,  Conserv. bottom T, {YAVRG[0]}-{YAVRG[-1]} mo: {MAVRG[0]}-{MAVRG[-1]}'

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
width  = 2200*1.e3
height = 2200*1.e3
lat0   = 62.
lon0   = -172.

m = Basemap(width=width, height=height, resolution='l',\
            projection='stere', lat_ts=55, lat_0=lat0, lon_0=lon0)

xR, yR = m(Glon, Glat)

btx = 'plot_cold_pool_glorys.py'
sttl = run_info

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()

ax1 = manseas.plot_stereogr_axis(fig1, m, xR, yR, Tbtm, clrmp, rmin, rmax, \
                       btx=btx, tscntrs=tscntrs, tslabels=tslabels, sttl=sttl)

plt.sca(ax1)
ax1.contour(xR,yR, Tbtm, [2], linestyles='solid', colors=[(1., 0.4, 0.9)])

ax1.contour(xR,yR,HH,[x for x in range(-8000,0,500)], linestyles='solid', colors=[(0.9,0.9,0.9)], linewidths=1)
# Show region:
ax1.contour(xR,yR, MSKBS, [0.9], linestyles='solid', colors=[(0.8,0.2,0)])




 
