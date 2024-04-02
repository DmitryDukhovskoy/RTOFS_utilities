"""
  Average bottom velocities along the western US coast
  GOFS reanalysis
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
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap

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
import mod_mom6 as mom6util
import mod_gofs31 as mgofs

# Region domain
lat1 = 22.
lat2 = 48.
lon1 = 235.42
lon2 = 251.

YR1  = 1994
YR2  = 2015
HR   = 12
MM   = 1
MD   = 1
dday = 1

# Western coast for GLBv grid:
is1 = 650
is2 = 886
js1 = 1775
js2 = 2200


pthoutp  = '/work/Dmitry.Dukhovskoy/data/gofs31_btmu_westcoast/'

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

uvarnm = 'water_u_bottom'
vvarnm = 'water_v_bottom'

# Average by months
istart = -1e30
tic  = timeit.default_timer()
for YR in range(YR1,YR2+1):
  for MM in range(1,12+1):
    ndays = mtime.month_days(MM, YR)
    ikk = 0
    ticR = timeit.default_timer()
    for MD in range(1,ndays+1):
      dnmb = mtime.datenum([YR,MM,MD]) 

      nexpt, _  = mgofs.gofs31_expt53X(dnmb)
      expt      = int(nexpt*10)

      urlB  = 'https://tds.hycom.org/thredds/dodsC/datasets/GLBv0.08/'
      urlF  = f'{urlB}expt_53.X/data/{YR}'
      urlT  = f'{urlB}expt_53.X/topo/'
      ftopo = 'depth_GLBv0.08_11.nc'
      fnmnc = f'hycom_GLBv0.08_{expt}_{YR}{MM:02d}{MD:02d}{HR:02d}_t000.nc'
      furl  = os.path.join(urlF, fnmnc)
      fturl = os.path.join(urlT, ftopo)
#nc = ncFile(fturl)
#lookup_ncvar(nc)
      if istart < 0:
        lon1d = read_field(furl, 'lon')
        lat1d = read_field(furl, 'lat')
        LON, LAT = np.meshgrid(lon1d, lat1d)
# Bathymetry - grid starts from 0W
        HH = read_field(fturl, 'bathymetry')
        HH = np.where(HH>1.e10, np.nan, HH)
        HH = -HH
        HH = np.where(HH == np.nan, 100., HH)
        mm = HH.shape[0]
        nn = HH.shape[1]  

        lat_topo = read_field(fturl, 'Latitude')
        lon_topo = read_field(fturl, 'Longitude')
        lon_topo = np.where(lon_topo>=180., lon_topo-360., lon_topo)
        istart   = np.argwhere(lon_topo == lon1d[0])[0][0]
        lon_topo = np.concatenate((lon_topo[istart:],lon_topo[:istart]))

        HH = np.concatenate((HH[:,istart:],HH[:,:istart]), axis=1)
 
# Subsample region:
#        if lon1 > 180:
#          lon1 = lon1-360.
#        if lon2 > 180:
#          lon2 = lon2-360.
#is1, js1 = mutil.find_indx_lonlat(lon1, lat1, LON, LAT)
#is1      = is1-1
#js1      = js1-1
#is2, js2 = mutil.find_indx_lonlat(lon2, lat2, LON, LAT)
#is2      = is2-1
#js2      = js2-1
        lon_s = lon1d[is1:is2+1]
        lat_s = lat1d[js1:js2+1]
        HH_s  = HH[js1:js2+1, is1:is2+1]

        ftopo_regn = "gofs31_GLBv008_topo11_westcoast.pkl"
        dftopo_regn = os.path.join(pthoutp,ftopo_regn)
        print(f"Saving region topo, grid --> {dftopo_regn}")
        with open(dftopo_regn,'wb') as fid:
          pickle.dump([lon_s, lat_s, HH_s], fid)

# Read bottom U:
# missing_value: -30000
      try:
        ubtm  = read_field(furl, uvarnm)
        vbtm  = read_field(furl, vvarnm)
      except:
        print(f"Could not read {furl}, ikk={ikk} skipping ...")
        continue

      ubtm  = np.where(ubtm < -100., np.nan, ubtm)
      vbtm  = np.where(vbtm < -100., np.nan, vbtm)

# Subsample region:
      ubtm_s = ubtm[js1:js2+1, is1:is2+1]
      vbtm_s = vbtm[js1:js2+1, is1:is2+1]
      jdm    = ubtm_s.shape[0]
      idm    = ubtm_s.shape[1]

# Sum up:
      if ikk == 0:
        UBTM = np.zeros((jdm,idm)) 
        VBTM = np.zeros((jdm,idm)) 

      ikk  += 1
      UBTM = UBTM + ubtm_s
      VBTM = VBTM + vbtm_s      
#
# Save monthly avearge:
    UBTM = UBTM/ikk
    VBTM = VBTM/ikk
    umin = np.nanmin(UBTM)
    umax = np.nanmax(UBTM)
    vmin = np.nanmin(VBTM)
    vmax = np.nanmax(VBTM)

    floutp  = f"gofs31_53X_btmuv_westcoast_{YR}{MM:02d}.pkl"
    dfloutp = os.path.join(pthoutp,floutp)

    toc    = timeit.default_timer()
    dtic   = (toc-ticR)/60.
    tictot = (toc-tic)/60.

    print(f"Saving mean UV --> {dfloutp}")
    with open(dfloutp,'wb') as fid:
      pickle.dump([UBTM,VBTM],fid)

    print(f"# of processed files = {ikk}")
    print(f"Min/max U: {umin:6.3f}/{umax:6.3f}")
    print(f"Min/max V: {vmin:6.3f}/{vmax:6.3f}")
    print(f"Processed {YR}/{MM} dt={dtic:6.2f} min, ttot={tictot:8.1f} min")

print("ALL DONE")

f_plt = False
if f_plt:
  cmpr = mutil.colormap_ssh(nclrs=200)
  rmin = -0.2
  rmax = 0.2
  cmpr.set_bad(color=[0.1, 0.1, 0.1])

  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im1 = ax1.pcolormesh(vbtm, \
                   cmap=cmpr,\
                   vmin=rmin, \
                   vmax=rmax)

  ax1.axis('scaled')
  ax1.set_xlim([600, is2])
  ax1.set_ylim([js1, js2])


 
