# 1st baroclinic Rossby radius from WOA23 
# Plot regional Rossby / max(dx, dy) of the grid cells
# gives # of grid cells per R
#
# GOFS3.1 in the NE Pacfic region
#
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#import torch
import sys
import pdb
import netCDF4
import importlib
from netCDF4 import Dataset as ncFile
import timeit
import pickle
import yaml


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
import mod_misc1 as mmsc1

grd=0.25
if grd==0.25:
  cgrd=4
woa='woa23'
seas=15    # season: 1-12 monthly, 13-winter (Jan-Mar), 14-spring (Apr-Jun), ...

f_deriveR = False # compute R/dx or load saved
f_showregn = 'NEP'

pthout  = '/work/Dmitry.Dukhovskoy/data/Rossby_WOA/'
fout1   = pthout + f'Rrossby_num_WOA23_season{seas:02d}.pkl'
foutRdx = pthout + f'Rrossby_dx_gofs.pkl'

btx = 'rrosby_gridres_GOFS.py'

with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

nrun         = "GOFS3.1"
expt         = "93.0"
pthgrid      = dct[nrun][expt]["pthgrid"]
ftopo        = dct[nrun][expt]["ftopo"]
fgrid        = dct[nrun][expt]["fgrid"]
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
DX, DY       = mhycom.dx_dy(LON,LAT)


pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom)
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom)
LONM, LATM  = mom6util.read_mom6grid(dfgrid_mom)

# Read saved Rossby R.
with open(fout1,'rb') as fid:
  RsbNum, LONW, LATW, LMsk = pickle.load(fid)

if seas == 13:
  cseas = 'Jan-Mar'
elif seas == 14:
  cseas = 'Apr-Jun'
elif seas == 15:
  cseas = 'Jul-Sep'
elif seas == 16:
  cseas = 'Oct-Dec'
else:
  cseas = f"{seas:02d}"

ny = RsbNum.shape[0]; nx = RsbNum.shape[1]


# Calc R/dx:
# NEP region:
lon1 = 145.
lon2 = 245.
lat1 = 22.
lat2 = 80.
if np.max(LONW)<lon2:
  LONW = np.where(LONW<0, LONW+360., LONW)

Iocn = np.where( (LMsk.flatten()>0) & (LONW.flatten()>=lon1) & \
                 (LONW.flatten()<=lon2) & (LATW.flatten()>=lat1) &\
                 (LATW.flatten()<=lat2) )[0]

nocn = len(Iocn)
cc   = 0
tic  = timeit.default_timer()
ticR = timeit.default_timer()
R2dx = RsbNum.copy()*0.0-100.

if f_deriveR:
  print("Calculating R/dx ...")
  for iocn in range(nocn):
    I1 = Iocn[iocn]
    jj, ii = np.unravel_index(I1, LMsk.shape)
    cc += 1
    
    if (cc % 2000) == 0:
      toc = timeit.default_timer()
      print(' {0:5.2f}% done {1:6.2f} min tot, {2:6.2f} min, max R2dx={3:8.4f}...'.\
              format(cc/nocn*100,(toc-tic)/60,(toc-ticR)/60, np.nanmax(R2dx)))
      ticR = timeit.default_timer()

    x0 = LONW[jj,ii]
    y0 = LATW[jj,ii]
    iH, jH = mutil.find_indx_lonlat(x0, y0, LON, LAT)

    dx = DX[jH,iH]*1.e-3
    dy = DY[jH,iH]*1.e-3
    dG = max([dx,dy])
    R2dx[jj,ii] = RsbNum[jj,ii]/dG

  print(f"Saving R2dx --> {foutRdx}")
  with open(foutRdx, 'wb') as fid:
    pickle.dump(R2dx, fid)

print(f"Loading R2dx <--- {foutRdx}")
with open(foutRdx, 'rb') as fid:
  R2dx = pickle.load(fid)

# Set up orthographic projection
from mpl_toolkits.basemap import Basemap, cm

# Add extra row/col for plotting
lonw = LONW[0,:]
lonw = np.append(lonw, lonw[0]+360)
latw = LATW[:,0]
latw = np.append(latw, 89.99)

lonw, latw = np.meshgrid(lonw, latw)


lon0 = 220.
lat0 = 60.
res  = 'l'
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
xR, yR = m(lonw,latw)

PMsk = ( (xR > 1e20) | (yR > 1e20) )
data = RsbNum.copy()
data = np.insert(data, 0, data[:,-1], axis=1)
data = np.insert(data, -1, data[-1,:], axis=0)
data[PMsk] = np.nan
xR[PMsk]   = 1.e30
yR[PMsk]   = 1.e30

data = data[0:ny, 0:nx]

if f_showregn == 'NEP':
# Read MOM6 NEP domain nx=342, ny=816
  pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
  ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
  fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
  dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom)
  dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom)
  LONM, LATM  = mom6util.read_mom6grid(dfgrid_mom)
  HHM         = mom6util.read_mom6depth(dftopo_mom)
  jdm         = np.shape(HHM)[0]
  idm         = np.shape(HHM)[1]

# Find boundaries of the NEP domain:
  IBND = np.array([])
  JBND = np.array([])
  dstp = 10
# Western bndry
  IBND = LONM[0:jdm:dstp, 0]
  JBND = LATM[0:jdm:dstp, 0]
# northern bndry
  IBND = np.append(IBND, LONM[jdm-1, 0:idm:dstp])
  JBND = np.append(JBND, LATM[jdm-1, 0:idm:dstp])
# eastern bndry
  IBND = np.append(IBND, LONM[0:jdm:dstp, idm-1])
  JBND = np.append(JBND, LATM[0:jdm:dstp, idm-1])
# southern
  IBND = np.append(IBND, LONM[0, 0:idm:dstp])
  JBND = np.append(JBND, LATM[0, 0:idm:dstp])

xBND = []
yBND = []
if f_showregn == 'NEP':
  xBND, yBND = m(IBND, JBND)

ctitle = f'Rossby Radius/||dl||$_\infty$, 0.08 GOFS3.1, {cseas}'
rmin = 0.
rmax = 7.


plt.ion()
cmpS = mcmp.colormap_salin()
cmpS.set_bad(color=[0.1, 0.1, 0.1])
cmpS.set_under(color=[0.85, 0.85, 0.85])

fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
plt.clf()

ax1 = plt.axes([0.04, 0.04, 0.8, 0.8])
im1 = m.pcolormesh(xR, yR, R2dx, shading='flat', cmap=cmpS,\
                   vmin=rmin, vmax=rmax)

m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

if len(xBND) > 0:
  m.plot(xBND, yBND, 'w.')

ax1.set_title(ctitle)

ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='max')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()


ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=5)

#fig1.colorbar(im,ax=ax1,orientation='horizontal')

bottom_text(btx, pos=[0.02, 0.02])


#from mod_plot_anls import zonalavrg
#ctl2='1st Barocl Rossby R zonal avrg, {0}, {1}'.format(tfnm,sfnm)
#zonalavrg(RsbNum,ctl2,lat,LMsk,btx=btx,ifg=1)

# Test: plot eigenfunctions
f_plt=0
if f_plt>0:
  plt.ion()
  fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
  plt.clf()

  im=6

  Rrsb = RsbNum[jj,ii]
  x0 = lon[ii]
  y0 = lat[jj]
  ww = W[im]
  vvk = V[:,im]
  zzk = Z_phi[0:kbtm+1]  
  nzk = zzk.shape[0]
# 
# Add surface and bottom to eig/functions
  zzV=np.zeros(nzk+2)
  zzV[0] = 0.  
  zzV[1:nzk+1]=zzk
  zzV[nzk+1]=zbtm

# Add 0 at the ends for eig/functions:
  vvV = np.zeros(zzV.shape[0])
  vvV[1:nzk+1] = vvk

  plt.plot(vvV,zzV,'.-')
  ctl = 'Eig/vector ii={2}, jj={3}, {4:5.2f}E, {5:5.2f}N, im={0}, Rr={1:6.0f} km'.\
         format(im,Rrsb,ii,jj,x0,y0)
  plt.title(ctl)





