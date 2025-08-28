# 1st baroclinic Rossby radius from WOA23 
# Plot regional Rossby / max(dx, dy) of the grid cells
# gives # of grid cells per R
#
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#import torch
import sys
import pdb
import importlib
import timeit
import xarray
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
plot_dx2 = False  # True plot effective resolution |dx^2 + dy^2|^1/2
                  # False - plot max(dx,dy) - horiz. grid spacing

f_deriveR = False # compute R/dx or load saved


pthout  = '/work/Dmitry.Dukhovskoy/data/Rossby_WOA/'
fout1   = os.path.join(pthout, f'Rrossby_num_WOA23_season{seas:02d}.pkl')

btx = 'rrosby_gridres_ARC.py'


# Arctic domain:
pthtopo_arc = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'
dfgrid = os.path.join(pthtopo_arc,'ocean_mask_ARC12.nc')
dsgrid_arc = xarray.open_dataset(dfgrid)
LONM = dsgrid_arc['x'].data
LATM = dsgrid_arc['y'].data

dftopo = os.path.join(pthtopo_arc,'ocean_topog_ARC12.nc')
dstopo = xarray.open_dataset(dftopo)
HHM    = -(dstopo['depth'].data)
jdm    = np.shape(HHM)[0]
idm    = np.shape(HHM)[1]


# Read grid resolution:
# dX, dY are on MOM "supergrid" - half grid points
DX, DY = mhycom.dx_dy(LONM, LATM)

# Effective resolution as norm-2
#RS2 = np.sqrt(DX**2 + DY**2)*1.e-3  # km
#mm, nn = RS2.shape

# Resolution as a max of DX, DY
RS = np.maximum(DX, DY)*1e-3
#
#RS1 = np.where(HHM >= 0., np.nan, RS1)
#RS2 = np.where(HHM >= 0., np.nan, RS2)


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

ny, nx = RsbNum.shape

if plot_dx2:
  # L-2 norm = sqrt(sum(x_i^2)) = ||x||_2 
  foutRdx = pthout + f'Rrossby_diag_dist_mom6arc.pkl' 
  ctitle = fr'Rossby Radius $\|dl\|_2$, NEP MOM6, {cseas}'
else:
   # L-inf norm = max|x_i|, i.e. max of (dx,dy) in this case
  foutRdx = pthout + f'Rrossby_dx_mom6arc.pkl'
  ctitle = fr'Rossby Radius/$\|dl\|_\infty$, NEP MOM6, {cseas}'

# Calc R/dx:
Iocn = np.where(HHM.flatten() < -5.)[0]
nocn = len(Iocn)
cc   = 0
tic  = timeit.default_timer()
ticR = timeit.default_timer()
R2dx = RsbNum.copy()*0.0-100.
Indx = RsbNum.copy()*0.0
if f_deriveR:
  print("Calculating R/dx ...")
  for iocn in range(nocn):
    I1 = Iocn[iocn]
    jj, ii = np.unravel_index(I1, HHM.shape)
    cc += 1
    
    x0 = LONM[jj,ii]
    y0 = LATM[jj,ii]
    iH, jH = mutil.find_indx_lonlat(x0, y0, LONW, LATW)

    dx = DX[jj,ii]*1.e-3
    dy = DY[jj,ii]*1.e-3
    if plot_dx2:
      dG = np.sqrt(dx*dx + dy*dy)
    else: 
      dG = max([dx,dy])
    R2dx[jH,iH] = RsbNum[jH,iH]/dG
    Indx[jH,iH] = 1

    if (cc % 2000) == 0:
      toc = timeit.default_timer()
      print(' {0:5.2f}% done {1:6.2f} min tot, {2:6.2f} min, max R2dx={3:8.4f}...'.\
              format(cc/nocn*100,(toc-tic)/60,(toc-ticR)/60, np.nanmax(R2dx)))
      print(f"dG={dG:6.2f} R={RsbNum[jH,iH]:6.1f} R2dx={R2dx[jH,iH]:7.4f}")
      ticR = timeit.default_timer()
   

  # ARC12 grid is coraser then the WOA23 (in degrees), this generates some gaps
  # Fill them as average values:
  dmm = R2dx.copy()
  dmm[:600,:]=np.nan
  dltxy = 5
  Igaps = np.where(dmm.flatten()<0.)[0]
  ngaps = len(Igaps)
  print(f'Filling gaps in Rossby R. field, {ngaps} pnts ...')
  ticR = timeit.default_timer()
  cc = 0
  for igg in Igaps:
    jj, ii = np.unravel_index(igg, R2dx.shape)
    R0 = RsbNum[jj,ii]
    cc += 1
    if R0 > 0.:
      xR = LONW[jj,ii]
      yR = LATW[jj,ii]
      iM, jM = mutil.find_indx_lonlat(xR, yR, LONM, LATM,dlt_err=7000, fatal_err=False)
      if iM < 0 or jM < 0:
        # point outside ARC domain
        continue

      dx = DX[jM,iM]*1.e-3
      dy = DY[jM,iM]*1.e-3
      if plot_dx2:
        dG = np.sqrt(dx*dx + dy*dy)
      else:
        dG = max([dx,dy])

      R2dx[jj,ii] = R0/dG

      if (cc % 2000) == 0:
        toc = timeit.default_timer()
        pdone = cc/float(ngaps)*100.
        print(f' {pdone:.2f}% done {(toc-ticR)/60.:.2f} min, dG={dG:.1f} R/dx={R2dx[jj,ii]:.4f}...')
        ticR = timeit.default_timer()


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


lon0 = 180.
lat0 = 70.
res  = 'l'
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
xR, yR = m(lonw,latw)

PMsk = ( (xR > 1e20) | (yR > 1e20) )
data = R2dx.copy()
data = np.insert(data, 0, data[:,-1], axis=1)
data = np.insert(data, -1, data[-1,:], axis=0)
data[PMsk] = np.nan
xR[PMsk]   = 1.e30
yR[PMsk]   = 1.e30

data = data[0:ny, 0:nx]

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

ax1.set_title(ctitle)

ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='max')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=10)

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





