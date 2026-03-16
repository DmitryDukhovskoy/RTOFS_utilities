"""
  Test ice redistribution using nonlinear contrained optimization
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import matplotlib.colors as colors
import random

from mod_utils_fig import bottom_text
import mod_sis2_relax as msisrlx
import mod_cice6_utils as mc6util
importlib.reload(mc6util)


#hicat = np.array([0., 0.64, 1.39, 2.47, 4.57, 50.])
hicat = np.array([0, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5, 3.0, 3.5, 50])
ICAT = hicat[:-1]  # not last value which is not actual ice cat

ncat = len(ICAT)

puny = 1e-11
nn=50
CI = np.zeros((nn))
HI = np.zeros((nn))
CC = np.zeros((nn,ncat))
HC = np.zeros((nn,ncat))
CCn = np.zeros((nn,ncat))
HCn = np.zeros((nn,ncat))
print(f"Calling redistribute_hice")
ifx = 10
ck_min = 1.e-2     # approximate ice concentration used for filling thinner ice cats 
#ck_min = 0.8e-1     # approximate ice concentration used for filling thinner ice cats 
for ii in range(nn):
  cice = random.uniform(0.,1.)
  hice = random.uniform(0.,4.)
  # Fix values for comparison:
  #if ii == ifx:
  #  cice = 0.95
  #  hice = 2.55
  print(f"ii={ii}, hice={hice:.4f} cice={cice:.4f}")
  # Old algorithm
  hcat, ccat = msisrlx.redistribute_hice(hice, cice, ICAT=ICAT, ck_min=ck_min)
  CI[ii] = cice
  HI[ii] = hice
  CC[ii,:] = ccat
  HC[ii,:] = hcat

  # New algorithm
  # Some guess of ice conc. distr. by cats
  ain_new = np.zeros((ncat))
  #ain_new[0] = ai_new
  ain_new = ain_new + cice / ncat
  sum_ain = np.sum(ain_new)
  ain_new = ain_new / sum_ain - 1.e-12
  ain_new = np.where(ain_new < puny, 0., ain_new)

  vin_new = np.zeros((ncat))
  vtot_target = hice
  dhi_min = 0.01  # min diff between cat ice thicknesses from 2 adjacent cats
  ain_min = 1.e-8    # lower bound of ain(n) to avoid zeros
  ain_new, vin_new = mc6util.adjust_thkncats_aice(ain_new, vin_new, vtot_target, \
                         hicat, dhi_min,  bnd_min=ain_min)

  CCn[ii,:] = ain_new
  HCn[ii,:] = vin_new 


from matplotlib.patches import Polygon
plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])

ii = ifx
ccat = CC[ii,:]
hcat = HC[ii,:]
hice = HI[ii]
cice = CI[ii]


sinfo=''
for ik in range(len(ccat)):
  txt = f'cat {ik+1}: hi={hcat[ik]:.3e}, ai={ccat[ik]:.3e}\n'
  sinfo = sinfo + txt

vol_tot = np.sum(ccat*hcat)
ai_tot = np.sum(cice)
txt = f'Total: vol_ice={vol_tot:.3e} m3/m2, iconc_tot={ai_tot:.3e}'
sinfo = sinfo + txt

#plt.bar(ICAT, ccat, color=[0.8,0.9,1], width=0.2)
#ax1.plot(ICAT, CC[ii,:],'-o')
dltE=0.01
clr=[0.5,0.8,1]
for kk in range(ncat):
  hmin = hicat[kk]+dltE
  hmax = hicat[kk+1]-dltE
  if kk == ncat-1:
    hmax = hmin + 1.
  # Make visible very small ccat lines:
  #ccat_plt = np.where(ccat<0.005, 0.005, ccat)
  verts = [(hmin,0),(hmin,ccat[kk]),(hmax,ccat[kk]),(hmax,0)]
  poly  = Polygon(verts, facecolor=clr, edgecolor=clr, zorder=5)
  ax1.add_patch(poly)
#  ax1.plot([hmin,hmax],[ccat[kk],ccat[kk]],'-',linewidth=2, color=[0.,0.5,0.9])

xup = ICAT[-1]+(ICAT[-1]-ICAT[-2])
ax1.set_xlim([0,xup])

stl = f"hice={hice:.4f}, cice={cice:.4f}, ck_min={ck_min:.3e}"
ax1.set_title(stl)
ax1.set_xticks(ICAT)
ax1.set_xlabel('Ice Cat min Thicknesses')
ax1.set_ylabel('partial area')
ax1.grid('on')

ax2 = plt.axes([0.1,0.1,0.8,0.25])
ax2.text(0.1,0.1,sinfo)
ax2.axis('off')

btx = 'test_ice_redistribute.py'
bottom_text(btx)


# Optimization:
ccat = CCn[ii,:]
hcat = HCn[ii,:]

fig2 = plt.figure(2,figsize=(9,8))
plt.clf()
ax21 = plt.axes([0.1, 0.4, 0.8, 0.5])

sinfo=''
for ik in range(len(ccat)):
  txt = f'cat {ik+1}: hi={hcat[ik]:.3e}, ai={ccat[ik]:.3e}\n'
  sinfo = sinfo + txt

vol_tot = np.sum(ccat*hcat)
ai_tot = np.sum(cice)
txt = f'Total: vol_ice={vol_tot:.3e} m3/m2, iconc_tot={ai_tot:.3e}'
sinfo = sinfo + txt


dltE=0.01
clr=[1,0.4,.2]
for kk in range(ncat):
  hmin = hicat[kk]+dltE
  hmax = hicat[kk+1]-dltE
  if kk == ncat-1:
    hmax = hmin + 1.
  verts = [(hmin,0),(hmin,ccat[kk]),(hmax,ccat[kk]),(hmax,0)]
  poly  = Polygon(verts, facecolor=clr, edgecolor=clr, zorder=5)
  ax21.add_patch(poly)
#  ax21.plot([hmin,hmax],[ccat[kk],ccat[kk]],'-',linewidth=2, color=[0.,0.5,0.9])

xup = ICAT[-1]+(ICAT[-1]-ICAT[-2])
ax21.set_xlim([0,xup])

stl = f"Optmz, hice={hice:.4f}, cice={cice:.4f}, ck_min={ck_min:.3e}"
ax21.set_title(stl)
ax21.set_xticks(ICAT)
ax21.set_xlabel('Ice Cat min Thicknesses')
ax21.set_ylabel('partial area')
ax21.grid('on')

ax22 = plt.axes([0.1,0.1,0.8,0.25])
ax22.text(0.1,0.1,sinfo)
ax22.axis('off')

btx = 'test_ice_redistribute.py'
bottom_text(btx)



