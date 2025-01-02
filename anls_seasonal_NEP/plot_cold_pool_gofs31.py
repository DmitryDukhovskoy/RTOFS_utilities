"""
  Plot cold pool region on the Bering Shelf from GOFS3.1 reanalysis
  Monthly average bottom T for NEP region extracted in btmT_NEPmonthly_gofs31.py
  GOFS3.1 reanalysis 
  https://data.hycom.org/datasets/GLBv0.08/expt_53.X/data/2010/
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import time
import timeit
import pickle
from copy import copy
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
import mod_misc1 as mmisc
import mod_anls_seas as manseas
import mod_plot_xsections as mxsct
import matplotlib as mtplt
import mod_regmom as mregmom
import mod_colormaps as mclrmps
import mod_gofs31 as mgofs

# Region domain
lat1 = 10.
lat2 = 81.
lon1 = 156.5
lon2 = 255.5

YAVRG = [x for x in range(2005,2015)]
#MAVRG = [1,2,3]  # months to average: Winter  JFM, Summer: JAS
#MAVRG = [4,5,6]  # months to average: Spring, AMJ
#MAVRG = [7,8,9]  # months to average: Winter  JFM, Summer: JAS
MAVRG = [10,11,12]  # months to average: Fall


pthoutp  = '/work/Dmitry.Dukhovskoy/GOFS3.1/gofs31_Tbtm_NEP/'

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

# Load Bering Shelf boundary, saved in plot_cold_pool_seas.py:
pthanls = pthseas['MOM6_NEP']['seasonal_daily']['pthanls'].format(expt_nmb=2)
dfnm_vrtx = os.path.join(pthanls,f'BeringShelf_region_vrtx.pkl')
if not os.path.isfile(dfnm_vrtx):
  raise Exception(f'Region boundary need to be saved in plot_cold_pool_seas.py {dfnm_vrtx}')

with open(dfnm_vrtx, 'rb') as fid:
  XYBND = pickle.load(fid)

Xbnd = XYBND[0]
Ybnd = XYBND[1]
Xbnd = np.where(Xbnd>180, Xbnd-360., Xbnd)

# Get GOFS3.1 topo, grid for NEP:
pthoutp  = '/work/Dmitry.Dukhovskoy/GOFS3.1/gofs31_Tbtm_NEP/'
ftopo_regn = "gofs31_GLBv008_topo11_NEP.pkl"
dftopo_regn = os.path.join(pthoutp,ftopo_regn)
with open(dftopo_regn,'rb') as fid:
  lon1d, lat1d, HH = pickle.load(fid)
lon1d = np.where(lon1d>180, lon1d-360., lon1d)

jdm  = len(lat1d)
idm  = len(lon1d)
LONW = np.zeros((jdm,idm))
LATW = np.zeros((jdm,idm))
for ii in range(idm):
  LATW[:,ii]=lat1d
for jj in range(jdm):
  LONW[jj,:]=lon1d

# Derive Regional mask for Bering Shelf:
print('Searching indices of the Bering Shelf region ...')
II, JJ = mmisc.find_closest_indx(Xbnd, Ybnd, LONW, LATW)

X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MS, _, _ = mmisc.inpolygon_v2(X, Y, II, JJ)  # 
JBS, IBS = np.where( (MS == 1) & (HH >= -250) & (HH < 0) ) #exclude deeep regions
MSKBS  = np.zeros((jdm,idm))
MSKBS[JBS,IBS] = 1

# Average monthly means:
# !!! Potential T !!! 
icc = 0
Tbtm = []
for YR in YAVRG:
  for MM in MAVRG:
    floutp  = f"gofs31_53X_Tbtm_NEP_{YR}{MM:02d}.pkl"
    dfloutp = os.path.join(pthoutp,floutp)
    print(f'Loading {dfloutp}')
    with open(dfloutp,'rb') as fid:
      T2d = pickle.load(fid)

    if icc == 0:
      Tbtm = T2d.copy()
    else:
      Tbtm = Tbtm + T2d

    icc += 1

Tbtm = Tbtm/icc

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
tscntrs = [x/10 for x in range(-20,80,10)]
tslabels = [x/10 for x in range(-20,80,10)]

run_info = f'GOFS3.1-53.X,  Pot. bottom T, {YAVRG[0]}-{YAVRG[-1]} mo: {MAVRG[0]}-{MAVRG[-1]}'

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
width  = 2200*1.e3
height = 2200*1.e3
lat0   = 62.
lon0   = -172.

m = Basemap(width=width, height=height, resolution='l',\
            projection='stere', lat_ts=55, lat_0=lat0, lon_0=lon0)

xR, yR = m(LONW, LATW)

btx = 'plot_cold_pool_gofs31.py'
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


 
