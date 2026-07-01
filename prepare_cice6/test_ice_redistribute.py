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
import argparse
from matplotlib.patches import Polygon

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from MyPython.mod_utils_fig import bottom_text
#from MyPython.mod_sis2_relax import redistribute_hice
from MyPython.mod_cice6_utils import adjust_thkncats_aice
import mod_sis2_relax as msrlx
importlib.reload(msrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--bplot", help="Field to plot as bar diagram, default=aice",
                    choices=['aice','vice'], default='aice',
                    type=str)
parser.add_argument("--ncat", help="N of ice thickness categories, default=10",
                    choices=[5,10], default=10, type=int)
parser.add_argument("--ax2", help="Plot the other var as lines on top bar diagr, default=0 (no)",
                   choices=[1,0], default=0, type=int)
args = parser.parse_args()

bplot = args.bplot
Nc0   = args.ncat
plot_ax2 = args.ax2 == 1

if Nc0 == 5:
  hicat = np.array([0., 0.64, 1.39, 2.47, 4.57, 50.])
elif Nc0 == 10:
  hicat = np.array([0, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5, 3.0, 3.5, 50])
ICAT = hicat[:-1]  # not last value which is not actual ice cat

ncat = len(ICAT)


plt.ion()
def plot_bar_diagr(fig1, ccat, hcat, vcat, hice, cice, sinfo, dltE, clr, clrln, stl):
  plt.clf()
  ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])

  # Create secondary Y axis
  ax2 = ax1.twinx()

  # x locations for varadd
  xvar = np.zeros(ncat)

  for kk in range(ncat):
    hmin = hicat[kk] + dltE
    hmax = hicat[kk+1] - dltE

    if kk == ncat-1:
      hmax = hmin + 1.

    xvar[kk] = 0.5*(hmin + hmax)

    if bplot == 'aice':
      yup = 1.1 * np.nanmax(ccat)
      verts = [(hmin,0),
               (hmin,ccat[kk]),
               (hmax,ccat[kk]),
               (hmax,0)]
      ylbl = 'Partial Area'
      varadd = vcat

    elif bplot == 'vice':
      yup = 1.1 * np.nanmax(vcat)
      verts = [(hmin,0),
               (hmin,vcat[kk]),
               (hmax,vcat[kk]),
               (hmax,0)]
      ylbl = 'Ice Volume (m$^3$/m$^2$_cell)'
      varadd = ccat

    poly = Polygon(verts,
                   facecolor=clr,
                   edgecolor=clr,
                   zorder=5)
    ax1.add_patch(poly)
  xup = ICAT[-1] + 1.2*(ICAT[-1]-ICAT[-2])
  ax1.set_xlim([0, xup])
  ax1.set_ylim([0, yup])

  ax1.set_title(stl)
  ax1.set_xticks(ICAT)
  ax1.tick_params(axis='y', colors=clr, labelsize=14)
  ax1.set_xlabel('Ice Cat min Thicknesses')
  ax1.set_ylabel(ylbl, color=clr, fontsize=12)
  ax1.grid('on')

  # Plot second variable
  if plot_ax2:
    ax2.plot(xvar, varadd,
             '-o',
             linewidth=2,
             color=clrln,
             markersize=5,
             zorder=10)
    yup2 = 1.2 * np.nanmax(varadd)
    ax2.set_ylim([0, yup2])

    if bplot == 'aice':
      ax2.set_ylabel('Ice Volume (m$^3$/m$^2$_cell)', color=clrln, fontsize=12)
    else:
      ax2.set_ylabel('Partial Area', color=clrln, fontsize=12)

    ax2.tick_params(axis='y', colors=clrln, labelsize=14)



  ax3 = plt.axes([0.1,0.1,0.8,0.25])
  ax3.text(0.1,0.1,sinfo)
  ax3.axis('off')

  btx = 'test_ice_redistribute.py'
  bottom_text(btx)



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
  # Thick-to-thin algorithm
  #hcat, ccat = redistribute_hice(hice, cice, ICAT=ICAT, ck_min=ck_min)
  hcat, ccat = msrlx.redistribute_hice(hice, cice, ICAT=ICAT, ck_min=ck_min)
  
  CI[ii] = cice
  HI[ii] = hice
  CC[ii,:] = ccat
  HC[ii,:] = hcat

  # New algorithm: non-linear iterative optimization
  # Some guess of ice conc. distr. by cats
  ain_new = np.zeros((ncat))
  #ain_new[0] = ai_new
  ain_new = ain_new + cice / ncat
  sum_ain = np.sum(ain_new)
  ain_new = np.where(ain_new < puny, 0., ain_new)

  vin_new = np.zeros((ncat))
  vtot_target = hice
  dhi_min = 0.01  # min diff between cat ice thicknesses from 2 adjacent cats
  ain_min = 1.e-8    # lower bound of ain(n) to avoid zeros
  ain_new, vin_new = adjust_thkncats_aice(ain_new, vin_new, vtot_target, \
                         hicat, dhi_min,  bnd_min=ain_min)

  CCn[ii,:] = ain_new
  hin_new = np.divide(vin_new, ain_new, out=np.zeros_like(vin_new, dtype=float), where=(ain_new != 0))
  HCn[ii,:] = hin_new 


# Plot 1 example
ii = ifx

# Thick-to-thin method:
ccat = CC[ii,:]  # ice concentration by cats
hcat = HC[ii,:]  # ice thickness by cats
hice = HI[ii]
cice = CI[ii]

# Do not show empty cats:
ccat[ccat < 1e-11] = np.nan
hcat[hcat < 1e-11] = np.nan
vcat = ccat * hcat

sinfo=''
for ik in range(len(ccat)):
  txt = f'cat {ik+1}: hi={hcat[ik]:.3e}, ai={ccat[ik]:.3e}\n'
  sinfo = sinfo + txt

vol_tot = np.nansum(vcat)
ai_tot = np.nansum(ccat)
txt = f'Total: vol_ice={vol_tot:.3e} m3/m2, iconc_tot={ai_tot:.3e}'
sinfo = sinfo + txt
  
stl1 = f"hice={hice:.4f}, cice={cice:.4f}, ck_min={ck_min:.3e}"

#plt.bar(ICAT, ccat, color=[0.8,0.9,1], width=0.2)
#ax1.plot(ICAT, CC[ii,:],'-o')
dltE=0.01
clr=[0.1,0.5,1]
clrln = [0.91, 0.4, 0]


fig1 = plt.figure(1,figsize=(9,8))
plot_bar_diagr(fig1, ccat, hcat, vcat, hice, cice, sinfo, dltE, clr, clrln, stl1)


# Optimization:
ccat_opt = CCn[ii,:]
hcat_opt = HCn[ii,:]
# Do not show empty cats:
ccat_opt[ccat_opt < 1e-11] = np.nan
hcat_opt[hcat_opt < 1e-11] = np.nan

vcat_opt = ccat_opt * hcat_opt

fig2 = plt.figure(2,figsize=(9,8))
plt.clf()

sinfo=''
for ik in range(len(ccat)):
  txt = f'cat {ik+1}: hi={hcat_opt[ik]:.3e}, ai={ccat_opt[ik]:.3e}\n'
  sinfo = sinfo + txt

vol_tot_opt = np.sum(ccat_opt * hcat_opt)
ai_tot_opt = np.sum(ccat_opt)
txt = f'Total: vol_ice={vol_tot_opt:.3e} m3/m2, iconc_tot={ai_tot_opt:.3e}'
sinfo = sinfo + txt
stl2 = f"Optmz, hice={hice:.4f}, cice={cice:.4f}, ck_min={ck_min:.3e}"

dltE=0.01
#clr=[1,0.4,.2]
plot_bar_diagr(fig2, ccat_opt, hcat_opt, vcat_opt, hice, cice, sinfo, dltE, clr, clrln, stl2)




