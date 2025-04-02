"""
  Plot monthly mean characteristics of the cold pool in the Bering Sea
  area by T classes
  overall T and S
  
  stat derived in calc_coldpool_stat_seasfcst.py 
  monthly mean and StDev  bottom T/S derived in calc_mnthlyTSbtm.py

  The coldpool is defined as:
  bottom water < 2C
  deep water (> 250 m) is excluded

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import pickle
from copy import copy
import matplotlib.colors as colors
from yaml import safe_load
import argparse

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
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')
import mod_swstate as msw
import conversions as gsw

from mod_utils_fig import bottom_text
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
importlib.reload(mutob)

parser = argparse.ArgumentParser()
parser.add_argument("--expt", help="experiment number: 1, ...", type=int)
parser.add_argument("--MMI", help="init month of f/cast, 1,4,7,10", type=int)
parser.add_argument("--YRS", help="start year of f/casts to derive cold poot stat, 1993, ...", type=int)
parser.add_argument("--YRE", help="end year of f/casts to derive cold pool stat, 1993, ...", type=int)
args = parser.parse_args()


# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
f_insitu = False    # convert to in situ, for shallow Bering shelf - minor difference

# Default values: 
YRS    = 1993 # year start of the forecast
YRE    = 2002
MMI    = 1
nens   = 1    # ens # for ensemble runs - =1 for 3D ocean fields
expt_nmb1 = 3  # =2 - seas. f/casts no irelax, =3 - seas. f/casts with ice relax
expt_nmb2 = 2  # =2 - seas. f/casts no irelax, =3 - seas. f/casts with ice relax

if MMI==1:
  YRS=np.max([1994,YRS])


if args.expt:
  expt_nmb = args.expt
if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI

expt_name = "seasonal_daily"
runname1   = f"NEPphys_frcst_dailyOB-expt{expt_nmb1:02d}"
runname2   = f"NEPphys_frcst_dailyOB-expt{expt_nmb2:02d}"

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

TCLASS = [-2., -1., 0., 1., 2.]
nTC    = len(TCLASS)


pthdump = pthseas["MOM6_NEP"][expt_name]["pthdump"] 
fldump1  = f'{runname1}_coldpool_stat_MI{MMI:02d}_{YRS}_{YRE}.pkl'
drfldump1 = os.path.join(pthdump, fldump1)
print(f'Loading {drfldump1}')
with open(drfldump1, 'rb') as fid:
  CP_AREA1,CP_TSTAT1,CP_SSTAT1,CP_TCLASS1 = pickle.load(fid)


fldump2  = f'{runname2}_coldpool_stat_MI{MMI:02d}_{YRS}_{YRE}.pkl'
drfldump2 = os.path.join(pthdump, fldump2)
print(f'Loading {drfldump2}')
with open(drfldump2, 'rb') as fid:
  CP_AREA2,CP_TSTAT2,CP_SSTAT2,CP_TCLASS2 = pickle.load(fid)

cff = 1.e-5
CPAmn1  = np.mean(CP_AREA1, axis=0)*cff
CPAmin1 = np.min(CP_AREA1, axis=0)*cff
CPAmax1 = np.max(CP_AREA1, axis=0)*cff
CPAmn2  = np.mean(CP_AREA2, axis=0)*cff
CPAmin2 = np.min(CP_AREA2, axis=0)*cff
CPAmax2 = np.max(CP_AREA2, axis=0)*cff

def construct_array(MN1,llim1, ulim1, MN2, llim2, ulim2):
  AA = np.array([MN1,llim1,ulim1])
  AA = np.expand_dims(AA, axis=0)
  a = np.array([MN2,llim2,ulim2])
  a = np.expand_dims(a, axis=0)
  AA = np.append(AA, a, axis=0)

  return AA

# Pool together statistics of cold pool area from 2 experiments:
AA = construct_array(CPAmn1,CPAmin1,CPAmax1,CPAmn2,CPAmin2,CPAmax2)

#AA = np.array([CPAmn1,CPAmin1,CPAmax1])
#AA = np.expand_dims(AA, axis=0)
#a = np.array([CPAmn2,CPAmin2,CPAmax2])
#a = np.expand_dims(a, axis=0)
#AA = np.append(AA, a, axis=0)

# Pool together stat of cold pool Mean Temp from 2 experiments
Tmd1  = np.median(CP_TSTAT1, axis=0)[:,0].squeeze()
Tlprc1 = np.min(CP_TSTAT1, axis=0)[:,0].squeeze()
Tuprc1 = np.max(CP_TSTAT1, axis=0)[:,0].squeeze()
Tmd2  = np.median(CP_TSTAT2, axis=0)[:,0].squeeze()
Tlprc2 = np.min(CP_TSTAT2, axis=0)[:,0].squeeze()
Tuprc2 = np.max(CP_TSTAT2, axis=0)[:,0].squeeze()
#Tlprc1 = np.mean(CP_TSTAT1, axis=0)[:,1].squeeze()
#Tuprc1 = np.mean(CP_TSTAT1, axis=0)[:,2].squeeze()
TT = construct_array(Tmd1,Tlprc1,Tuprc1,Tmd2,Tlprc2,Tuprc2)

CLRS = [[0., 0.4, 0.8],
        [0.9, 0.6, 0.]]


# Plot stat:
plt.ion()

btx = 'calc_coldpool_stat_seasfcst.py'

from matplotlib.patches import Polygon
def plot_2bars(ax1, AA, CLRS):
  XM = np.arange(MMI, MMI+12)
  dx=0.35
  dtg=0.06
  clrln = [0,0,0]

  A1 = AA[0,0,:].squeeze()
  A1mn = AA[0,1,:].squeeze()
  A1mx = AA[0,2,:].squeeze()
  A2 = AA[1,0,:].squeeze()
  A2mn = AA[1,1,:].squeeze()
  A2mx = AA[1,2,:].squeeze()
  
  for imm in range(12):
    x0 = XM[imm]
    x1 = x0-dx
    x2 = x0+dx
    
    mn1 = A1[imm]
    mn2 = A2[imm]
    min1 = A1mn[imm]
    max1 = A1mx[imm]
    min2 = A2mn[imm]
    max2 = A2mx[imm]
    verts = [(x1, 0), (x1, mn1), (x0, mn1), (x0, 0)]
    fclr1 = CLRS[0]
    poly1 = Polygon(verts, facecolor=fclr1)
    ax1.add_patch(poly1)
    xm0 = x0-dx/2
    ax1.plot([xm0, xm0],[min1,max1],'-',linewidth=2, color=clrln)
    ax1.plot([xm0-dtg,xm0+dtg],[max1,max1],'-',linewidth=2, color=clrln)
    ax1.plot([xm0-dtg,xm0+dtg],[min1,min1],'-',linewidth=2, color=clrln)

    verts = [(x0, 0), (x0, mn2), (x2, mn2), (x2, 0)]
    fclr2 = CLRS[1]
    poly2 = Polygon(verts, facecolor=fclr2)
    ax1.add_patch(poly2)
    xm0 = x0+dx/2
    ax1.plot([xm0, xm0],[min1,max1],'-',linewidth=2, color=clrln)
    ax1.plot([xm0-dtg,xm0+dtg],[max1,max1],'-',linewidth=2, color=clrln)
    ax1.plot([xm0-dtg,xm0+dtg],[min1,min1],'-',linewidth=2, color=clrln)

  ax1.set_xticks(XM)
  ax1.set_xlim([XM[0]-1,XM[-1]+1])

  return(ax1)

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.55, 0.8, 0.35])

ax1 = plot_2bars(ax1, AA, CLRS)

ax1.set_ylim([0, 10])
ax1.grid('on')
ax1.set_xlabel('F/cast months')
ax1.set_ylabel('km2 x 1.e-5')
sttl = f'F/casts irlx(blue) & no irlx(orange), Minit={MMI:02d}, {YRS}-{YRE}, CP area'
ax1.set_title(sttl)

# Plot T stat:
ax2 = plt.axes([0.1, 0.1, 0.8, 0.35])
ax2 = plot_2bars(ax2, TT, CLRS)
ax2.set_ylim([-2, 2])
ax2.grid('on')
ax2.set_xlabel('F/cast months')
ax2.set_ylabel('T deg C')
sttl = f'F/casts irlx(blue) & no irlx(orange), Minit={MMI:02d}, {YRS}-{YRE}, coldpool T'
ax2.set_title(sttl)

bottom_text(btx, pos=[0.02,0.02])



f_showreg = False
if f_showreg:
# Stereographic Map projection:
  from mpl_toolkits.basemap import Basemap, cm
  m = Basemap(width=3300*1.e3,height=3300*1.e3, resolution='l',\
              projection='stere', lat_ts=55, lat_0=62, lon_0=-175)

  xR, yR = m(hlon, hlat)

  xcMap, ycMap = m(hlon[JJ,II], hlat[JJ,II])

  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = m.pcolormesh(xR, yR, Tbtm, cmap=clrmp, vmin=rmin, vmax=rmax)
  m.plot(xcMap, ycMap, '-', color=[0.8, 0.4, 0])
#  m.contour(xR, yR, HH, [-1000], colors=[(0,0,0)], linestyles='solid')
#  m.contour(xR, yR, Tbtm, [2.], colors=[(0,0.5,1.)], linestyles='solid')
  #ax1.axis('scaled')
  ax1.set_title(sttl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  # extend: min, max, both
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  ax3 = fig1.add_axes([0.02, 0.02, 0.8, 0.05])
  ax3.text(0, 0, sinfo, fontsize=10)
  ax3.axis('off')

  btx = 'plot_coldpool.py'
  bottom_text(btx, pos=[0.2, 0.01])



