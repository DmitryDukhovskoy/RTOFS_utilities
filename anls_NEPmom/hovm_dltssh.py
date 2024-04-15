# Track anomlaies of the ssh fields 
# propagating along the US West coast
# Hovmuller diagrams
# mean ssh precomputed in mean_ssh.py
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

YR       = 1993
pfld     = 'ssh'
f_derive = False
expt     = 'NEP_BGCphys'
pthpkl   = '/work/Dmitry.Dukhovskoy/data/mom6_nep_tmp/'
fldump   = f"{pthpkl}{expt}_meanSSH_{YR}.pkl"

btx = 'hovm_dltssh.py'

with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom)
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom)
LONM, LATM  = mom6util.read_mom6grid(dfgrid_mom)
HHM         = mom6util.read_mom6depth(dftopo_mom)
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]

DX, DY = mhycom.dx_dy(LONM, LATM)

# Get coast line:
#import mod_gulfstream as mgulf
#TCNT = mgulf.derive_contour(HHM, tz0=0., xFS=329, yFS=25,\
#                            xl1=85, yl1=20, xl2=340, yl2=669)

# Derive grid indices on the shelf:
import mod_anlsnep as mnep
importlib.reload(mnep)
z0 = -400.
ShIJ, IJl, IJw = mnep.find_NWP_shelf(HHM, LONM, LATM, DX, DY, z0)
ISH  = ShIJ.I
JSH  = ShIJ.J
Ldst = ShIJ.Dst 
nij  = len(ISH)

with open(fldump, 'rb') as fid:
  Amn = pickle.load(fid)

flssh  = f"expt_dltSSHshelf_{YR}.pkl"
dflssh = os.path.join(pthpkl, flssh)

if f_derive:
  ndays = 365
  icc   = -1
  DSSH  = np.zeros((nij, ndays))*np.nan
  for jday in range(1, ndays+1):
    dnmb   = mtime.jday2dnmb(YR, jday)
    dv     = mtime.datevec(dnmb)
    ds     = mtime.datestr(dnmb)
    MM     = dv[1]
    DM     = dv[2]

    print(f'Reading {jday} {YR}/{MM}/{DM}')
    pthout = f'/work/Dmitry.Dukhovskoy/run_output/NEP_BGCphys/{YR}/{MM:02d}/'
    flmom  = f'ocean_{YR}_{jday:03d}_12.nc'
    dflmom = os.path.join(pthout,flmom)
    if not os.path.isfile(dflmom):
      print(f"File not found skipping: {dflmom}")
      continue

    nc     = ncFile(dflmom,'r')
    A2d    = nc.variables[pfld][:].data.squeeze()
    A2d    = np.where(A2d > 1.e10, np.nan, A2d)
    A2d    = A2d - Amn

    icc += 1
  # Extract ssh anomalies on the shelf:
    for ik in range(nij):
      iil  = list(ISH[ik])
      jjl  = list(JSH[ik])
      dltA = A2d[jjl,iil]
      inan = np.where(~np.isnan(dltA))[0]

      if len(inan) == 0:
        print(f"ik={ik} No shelf values - all land check indices ISH JSH")
      elif len(inan) == 1:
        mda = dltA[inan[0]]
      else:
        mda  = np.nanmedian(dltA) 
      DSSH[ik, icc] = mda
      
    if (icc%50) == 0:
      print(f"Saving shelf ssh --> {dflssh}")
      with open(dflssh, 'wb') as fid:
        pickle.dump([DSSH, Ldst], fid)
    

  print(f"Saving shelf ssh --> {dflssh}")
  with open(dflssh, 'wb') as fid:
    pickle.dump([DSSH, Ldst], fid)

else:
  print(f"Loading {dflssh}")
  with open(dflssh,'rb') as fid:
    DSSH, Ldst = pickle.load(fid)

# Orient: Time x Distance
DSSH = DSSH.T

ctitle = f'{expt} MOM6 dltSSH {YR}/{MM:02d}/{DM:02d}'
cmpr = mutil.colormap_ssh(nclrs=200)
rmin = -0.2
rmax = 0.2

#================
#    Plotting
Xdst  = np.round(np.cumsum(Ldst)*1.e-3)
Ytime = np.arange(1,367)
Xdst  = np.append(Xdst,Xdst[-1]+np.round(Ldst[-1]*1.e-3))


plt.ion()

# Function to print mouse click event coordinates
#def onclick(event):
#   print([event.xdata, event.ydata])

fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
plt.clf()

ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])
im1 = ax1.pcolormesh(Xdst, Ytime, DSSH, cmap=cmpr, \
                     vmin=rmin, vmax=rmax)
ax1.set_title(ctitle)

ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()


ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

#fig1.colorbar(im,ax=ax1,orientation='horizontal')

bottom_text(btx, pos=[0.02, 0.02])


fig2 = plt.figure(2,figsize=(8,8), constrained_layout=False)
plt.clf()
ax21 = plt.axes([0.08, 0.1, 0.8, 0.8])
ax21.contour(HHM,[0], colors=[(0,0,0)])
ax21.contour(HHM, [-400], linestyles='solid', colors=[(0,0.6,0.8)])
ax21.contour(HHM, [-8000,-7000,-6000,-5000,-4000,-3000,-2000,-1000], 
                  linestyles='solid', colors=[(0.8,0.8,0.8)])
ax21.plot(IJl[:,0],IJl[:,1],'-',color=[0.6, 0.3, 0])

dX    = 500.
xStop = 0. 
for ikk in range(nij):
  i0  = ISH[ikk][-1]
  j0  = JSH[ikk][-1]
  dst = int(Xdst[ikk])
  if dst >= xStop:
    xStop += dX
    ax21.plot(i0,j0,'.',color=[1,0,0])
    ax21.text(i0, j0, f' {dst}')  

bottom_text(btx, pos=[0.02, 0.02])


f_chck = False
if f_chck
  fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
  plt.clf()

  ax1 = plt.axes([0.05, 0.05, 0.8, 0.8])
  ax1.contour(HHM, [0], colors=[(0,0,0)])
  ax1.contour(HHM, [-400], linestyles='solid', colors=[(0,0.6,0.8)])
  ax1.plot(IJl[:,0], IJl[:,1], '-', color=[0., 0.8, 0.2])
  ax1.plot(IJw[:,0], IJw[:,1], '-', color=[0., 0.4, 1])

  nij = len(ISH)
  for ii in range(nij):
    iil = ISH[ii]
    jjl = JSH[ii]

    ax1.plot(iil,jjl,'-', color=[0.5, 0.5, 0.5])
  
  




