"""
  Calc mean ice thickness in the BG
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
from copy import copy
import matplotlib.colors as colors
from yaml import safe_load
import argparse
import pandas as pd
import time 

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
import mod_interp1D as mint1d
import mod_oras as moras
importlib.reload(moras)


parser = argparse.ArgumentParser()
parser.add_argument("--yrs", help="start year", type=int, required=True)
parser.add_argument("--yre", help="end year", type=int, required=True)
args = parser.parse_args()

f_save = True  # save calculated characteristics

YRS = args.yrs if args.yrs else None
YRE = args.yre if args.yre else None

# Domain region:
xlim1 = 158
xlim2 = 670
ylim1 = 800
ylim2 = 1020

# Look for max A2d in the BG regions:
# Exclude shallow regions/shelves
#iBG1 = 441
#jBG1 = 899
#iBG2 = 521
#jBG2 = 970
IIM = [430, 430, 490, 500, 500, 490, 430]
JJM = [905, 980, 980, 950, 930, 902, 902]

# Region where to find closed Contours:
# Look for max A2d in the BG regions:
iC1 = 420
jC1 = 895
iC2 = 545
jC2 = 995
IBG = [iC1,iC1,iC2,iC2]
JBG = [jC1,jC2,jC2,jC1]

pthoras = '/work/Dmitry.Dukhovskoy/data/ORAS5/ITHKN'
pthout = '/work/Dmitry.Dukhovskoy/anls_output/oras5/BG_anls'

nyrs = YRE+1-YRS
icc = 0
Ithkn_mn = []
Ithkn_md = []
Ithkn_lp = []
Ithkn_up = []
TM = []
for YR in range(YRS,YRE+1):
  for MM in range(1,13):
    timeS = time.time()
    flnm = f'iicethic_control_monthly_highres_2D_{YR}{MM:02d}_CONS_v0.1.nc'
    if YR > 2014:
      flnm = f'iicethic_control_monthly_highres_2D_{YR}{MM:02d}_OPER_v0.1.nc'
    dfl = os.path.join(pthoras,flnm)
    print(f'Opening {dfl}')
    varnm='iicethic'
    dset = xarray.open_dataset(dfl)
    A2d = dset[varnm].data.squeeze()

    dnmb0 = mtime.datenum([YR,MM,15])
    TM.append(dnmb0)

    if icc == 0:
      jdm, idm = A2d.shape
      LMSK = np.where(np.isnan(A2d),0,1)
      X, Y     = np.meshgrid(np.arange(idm), np.arange(jdm))
      MS, _, _ = mmisc.inpolygon_v2(X, Y, IBG, JBG)  # 
      JBS, IBS = np.where( (MS == 1) & (LMSK == 1) ) #exclude deeep regions
      MSKBS  = np.zeros((jdm,idm))
      MSKBS[JBS,IBS] = 1
      # Mask for finding max_A2d in the BG:
      MSH, _, _ = mmisc.inpolygon_v2(X, Y, IIM, JJM)  # 
      #JA2d, IA2d = np.where((MSH==1) & (LMSK == 1))
      MSK = np.where((MSH==1) & (LMSK == 1), 1, 0)

      lonh = dset['nav_lon'].data
      lath = dset['nav_lat'].data
      # Mask for contouring A2d:
      # MSC, _, _, = mmisc.inpolygon_v2(X,Y,IC,JC)
      DX, DY = mmom6.dx_dy(lonh, lath)
      Acell  = DX*DY*1.e-6  # km2
      Acell_BG = Acell[jC1:jC2,iC1:iC2]
            
    # A2d gradient, BG area based on the last closed contour:
    # Use subset domain
    ithk_BG = A2d[jC1:jC2,iC1:iC2]
    MSK_BG = MSK[jC1:jC2,iC1:iC2]
    ithk_BG = np.where(MSK_BG==0, np.nan, ithk_BG)
    ithk_md = np.nanmedian(ithk_BG)
    ithk_mn = np.nanmean(ithk_BG)
    ithk_up = np.percentile(ithk_BG,90)
    ithk_lp = np.percentile(ithk_BG,10)
    print(f'BG mean thkn = {ithk_mn:.1f} m, med = {ithk_md:.1f}m')

    Ithkn_mn.append(ithk_mn)
    Ithkn_md.append(ithk_md)
    Ithkn_lp.append(ithk_lp)
    Ithkn_up.append(ithk_up)
    icc += 1
    timeE = time.time()
    print(f'Elapsed time: {(timeE-timeS)*1./60.:.3f} min ')

Ithkn_mn = np.array(Ithkn_mn)
Ithkn_md = np.array(Ithkn_md)
Ithkn_up = np.array(Ithkn_up)
Ithkn_lp = np.array(Ithkn_lp)
TM = np.array(TM)

if f_save:
  fout = f'BG_ithkn_{YRS}-{YRE}.npz'
  dfout = os.path.join(pthout,fout)
  print(f'Saving ---> {dfout}')
  np.savez(dfout, TM=TM, thkmd=Ithkn_md, thkmn=Ithkn_mn, thklp=Ithkn_lp, thkup=Ithkn_up)     

f_plt = False
if f_plt:
  clrmp = mclrmps.colormap_A2d(cpos='YlOrRd', cneg='PuBuGn_r')
  rmin = -0.5
  rmax = 0.5

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  ax1.plot(

  img = ax1.pcolormesh(A2d_dmn, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.axis('scaled')
  ax1.set_xlim([xlim1,xlim2])
  ax1.set_ylim([ylim1,ylim2])
  ax1.plot([iC1,iC2],[jC1,jC1],'-',color=[0.,0.6,0.8])
  ax1.plot([iC1,iC2],[jC2,jC2],'-',color=[0.,0.6,0.8])
  ax1.plot([iC1,iC1],[jC1,jC2],'-',color=[0.,0.6,0.8])
  ax1.plot([iC2,iC2],[jC1,jC2],'-',color=[0.,0.6,0.8])

  ax1.plot(IIM,JJM,'-')

  sttl = f'A2d demeaned, ORAS5, {YR}/{MM:02d}'
  ax1.set_title(sttl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  # extend: min, max, both
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'calc_ithknBG.py'
  bottom_text(btx) 

