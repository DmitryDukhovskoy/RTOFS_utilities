"""
  Monthly average bottom T for NEP region 
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

# Specify field: temp or salin:
varread = 'salin'

# Region domain
#lat1 = 10.
#lat2 = 81.
#lon1 = 156.5
#lon2 = 255.5

YR1  = 2004
YR2  = 2015
HR   = 12
dday = 5

match varread:
  case 'temp':
    varnm = 'water_temp_bottom'
    pthoutp  = '/work/Dmitry.Dukhovskoy/GOFS3.1/gofs31_TSbtm_NEP/'
    flprfx  = f"gofs31_53X_Tbtm_NEP"
  case 'salin':
    varnm = 'salinity_bottom'
    pthoutp  = '/work/Dmitry.Dukhovskoy/GOFS3.1/gofs31_Sbtm_NEP/'
    flprfx  = f"gofs31_53X_Sbtm_NEP"

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


# Average by months
istart = -1e30
tic  = timeit.default_timer()
for YR in range(YR1,YR2+1):
  for MM in range(1,12+1):
    ndays = mtime.month_days(MM, YR)
    ikk = 0
    ticR = timeit.default_timer()
    for MD in range(1,ndays+1, dday):
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
      if istart < 0:
        # Bathymetry - grid starts from 0W
        HH = read_field(fturl, 'bathymetry')
        HH = np.where(HH>1.e10, np.nan, HH)
        HH = -HH
        HH = np.where(np.isnan(HH), 100., HH)
        mm = HH.shape[0]
        nn = HH.shape[1]  

        lat_topo = read_field(fturl, 'Latitude')
        lon_topo = read_field(fturl, 'Longitude')
        lon_topo = np.where(lon_topo>=180., lon_topo-360., lon_topo)
# Data grid starts from -180W
# reshaffle data to start from 0
        lon1d = read_field(furl, 'lon')
        lat1d = read_field(furl, 'lat')
        istart = np.argwhere(lon1d == lon_topo[0])[0][0]
        lon1d = np.concatenate((lon1d[istart:], lon1d[:istart]))
        LON, LAT = np.meshgrid(lon1d, lat1d)
# Make grid start from -180W
#        istart   = np.argwhere(lon_topo == lon1d[0])[0][0]
#        lon_topo = np.concatenate((lon_topo[istart:],lon_topo[:istart]))
#        HH = np.concatenate((HH[:,istart:],HH[:,:istart]), axis=1)
# Subsample region:
        is1 = 1940
        is2 = 3234
        js1 = 1740
        js2 = 3100
        lon_s = lon1d[is1:is2+1]
        lat_s = lat1d[js1:js2+1]
        HH_s  = HH[js1:js2+1, is1:is2+1]

        ftopo_regn = "gofs31_GLBv008_topo11_NEP.pkl"
        dftopo_regn = os.path.join(pthoutp,ftopo_regn)
        if not os.path.isfile(dftopo_regn): 
          print(f"Saving region topo, grid --> {dftopo_regn}")
          with open(dftopo_regn,'wb') as fid:
            pickle.dump([lon_s, lat_s, HH_s], fid)

# Read bottom temp:
# missing_value: -30000
      TSbtm = []
      try:
        TSbtm  = read_field(furl, varnm)
      except:
        print(f"Could not read {furl}, ikk={ikk} skipping ...")
        continue

      TSbtm  = np.where(abs(TSbtm) > 100., np.nan, TSbtm)
      # Reshuffle TSbtm/Sbtm to match topo:
      TSbtm = np.concatenate((TSbtm[:,istart:],TSbtm[:,:istart]), axis=1)

      # Subsample region:
      TSbtm_s = TSbtm[js1:js2+1, is1:is2+1]
      jdm    = TSbtm_s.shape[0]
      idm    = TSbtm_s.shape[1]

      # Sum up:
      if ikk == 0:
        TSBTM = np.zeros((jdm,idm)) 
      ikk  += 1
      TSBTM = TSBTM + TSbtm_s
#
# Save monthly avearge:
    TSBTM = TSBTM/ikk
    tmin = np.nanmin(TSBTM)
    tmax = np.nanmax(TSBTM)

    floutp  = f"{flprfx}_{YR}{MM:02d}.pkl"
    dfloutp = os.path.join(pthoutp,floutp)

    toc    = timeit.default_timer()
    dtic   = (toc-ticR)/60.
    tictot = (toc-tic)/60.
#
#    A = STOP
#
    print(f"Saving mean {varnm} --> {dfloutp}")
    with open(dfloutp,'wb') as fid:
      pickle.dump(TSBTM,fid)

    print(f"# of processed files = {ikk}")
    print(f"Min/max {varnm}: {tmin:6.3f}/{tmax:6.3f}")
    print(f'Check {varnm}(812,420)={TSBTM[812,420]:.2f}')
    print(f"Processed {YR}/{MM} dt={dtic:6.2f} min, ttot={tictot:8.1f} min")

print("ALL DONE")

f_plt = False
if f_plt:
  cmpr = mutil.colormap_ssh(nclrs=200)
  rmin = -2
  rmax = 5
  cmpr.set_bad(color=[0.1, 0.1, 0.1])

  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im1 = ax1.pcolormesh(TSbtm_s, \
                   cmap=cmpr,\
                   vmin=rmin, \
                   vmax=rmax)

  ax1.axis('scaled')
  ax1.set_xlim([600, is2])
  ax1.set_ylim([js1, js2])


 
