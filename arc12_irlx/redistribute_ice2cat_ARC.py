"""
  Redistribute  relax ice fields by thcikness categories
  ice categories are hrad-coded in SIS_state_initialization.F90
  these can be changed in SIS_override
  real :: hlim_dflt(8) = (/ 1.0e-10, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5 /) ! lower thickness limits 1...CatIce

  note:   nCat_dflt = 5 ; if (slab_ice) nCat_dflt = 1
  and SIS_input/ SIS_override: NCAT_ICE = 5 (if not then use default)
  so only 5 categories by default are used

"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
import pickle
from yaml import safe_load

PPTHN = '/home/Dmitry.Dukhovskoy/python'
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc
from mod_utils_fig import bottom_text
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

ifld = 'ithkn'  # ithkn, iarea
file_type = 'monthly'  # monthly, daily, ... or clim
                       # for climatologies, do not need padded time - data will be recycled
                       # for monthly, daily, etc. need -dt and +dt at the beginn/end 
""" 
ICAT in ARC:
 check: CatIce=10
distribute_ice2cats thkn cat 1 hLim= 0.000
distribute_ice2cats thkn cat 2 hLim= 0.100
distribute_ice2cats thkn cat 3 hLim= 0.300
distribute_ice2cats thkn cat 4 hLim= 0.700
distribute_ice2cats thkn cat 5 hLim= 1.100
distribute_ice2cats thkn cat 6 hLim= 1.500
distribute_ice2cats thkn cat 7 hLim= 2.000
distribute_ice2cats thkn cat 8 hLim= 2.500
distribute_ice2cats thkn cat 9 hLim= 3.000
distribute_ice2cats thkn cat 10 hLim= 3.500
distribute_ice2cats thkn cat 11 hLim= 4.000
"""
# ICAT in NEP:
#ICAT = np.array([1.0e-10, 0.1, 0.3, 0.7, 1.1])
ICAT = np.array([1.0e-10, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5, 3.0, 3.5])

pthdata = '/work/Dmitry.Dukhovskoy/ARC12/irlx'
pthtopo = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'
frlx = 'nudging_ice.nc'
dfrlx = os.path.join(pthdata,frlx)
ds_irlx = xarray.open_dataset(dfrlx)

MM = 3
itime = MM-1
Cice = ds_irlx['sic_rg'].isel(time=itime).data
Hice = ds_irlx['sit_rg'].isel(time=itime).data

dtopo = os.path.join(pthtopo,'ocean_topog.nc')
ds_topo = xarray.open_dataset(dtopo)
HH = -ds_topo['depth'].data
jdm, idm = HH.shape

hlon = ds_irlx['lon'].data
hlat = ds_irlx['lat'].data

Hice = np.where(HH>=0, np.nan, Hice)
Cice = np.where(HH>=0, np.nan, Cice)


#i0 = 261
#j0 = 728
i0 = 188
j0 = 0
hice = Hice[j0,i0]
cice = Cice[j0,i0]

ICAT0 = ICAT
ncat = len(ICAT0)
hice = 0.2895118
cice = 0.0415798
hcat, ccat = msisrlx.redistribute_hice(hice, cice, ICAT=ICAT0, ck_min=1.e-2)
# Check conservation:
htot = np.sum(hcat*ccat)
ctot = np.sum(ccat)

# Test:
import random

nn=500
CI = np.zeros((nn))
HI = np.zeros((nn))
CC = np.zeros((nn,ncat))
HC = np.zeros((nn,ncat))
ERRH = np.zeros((nn))
ERRC = np.zeros((nn))
print(f"Calling redistribute_hice")
for ii in range(nn):
  cice = random.uniform(0.,.5)
  hice = random.uniform(0.,.3)
  print(f"ii={ii}, hice={hice:.4f} cice={cice:.4f}")
  hcat, ccat = msisrlx.redistribute_hice(hice, cice, ICAT=ICAT0, ck_min=1.e-2)
# Check conservation:
  htot = np.sum(hcat*ccat)
  ctot = np.sum(ccat)
  eh = np.abs(hice-htot)
  ec = np.abs(cice-ctot)
  if eh>1.e-6:
    print(f'ERR: not conserved h: htot={htot} hice={hice}')
  if ec>1.e-6:
    print(f'ERR: not conserved c: ctot={ctot} cice={cice}')

  CI[ii] = cice
  HI[ii] = hice
  CC[ii,:] = ccat
  HC[ii,:] = hcat 
  ERRH[ii] = eh
  ERRC[ii] = ec


ICATK = np.append(ICAT0,[3])

from matplotlib.patches import Polygon
plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])

ii = 10
ccat = CC[ii,:]
hcat = HC[ii,:]
hice = HI[ii]
cice = CI[ii]
#plt.bar(ICAT0, ccat, color=[0.8,0.9,1], width=0.2)
#ax1.plot(ICAT0, CC[ii,:],'-o')
dltE=0.01
clr=[0.5,0.8,1]
for kk in range(ncat):
  hmin = ICATK[kk]+dltE
  hmax = ICATK[kk+1]-dltE
  verts = [(hmin,0),(hmin,ccat[kk]),(hmax,ccat[kk]),(hmax,0)]
  poly  = Polygon(verts, facecolor=clr, edgecolor=clr, zorder=5)
  ax1.add_patch(poly)
#  ax1.plot([hmin,hmax],[ccat[kk],ccat[kk]],'-',linewidth=2, color=[0.,0.5,0.9])

ax1.set_xlim([0,1.5])

stl = f"hice={hice:.4f}, cice={cice:.4f}"
ax1.set_title(stl)
ax1.set_xticks(ICAT0)
ax1.set_xlabel('Ice Cat min Thicknesses')
ax1.set_ylabel('partial area')
ax1.grid('on')

btx = 'redistribute_ice2cat.py'
bottom_text(btx)



