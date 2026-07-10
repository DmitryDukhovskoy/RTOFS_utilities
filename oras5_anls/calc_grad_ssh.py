"""
  Calc ssh gradient that determines the BG intensity
  following Prosh & Johnson
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
parser.add_argument("--yrs", help="start year", type=int)
parser.add_argument("--yre", help="end year", type=int)
args = parser.parse_args()

f_save = False  # save calculated characteristics

if args.yrs:
  YRS = args.yrs
  YRE = YRS
if args.yre:
  YRE = args.yre

# Domain region:
xlim1 = 158
xlim2 = 670
ylim1 = 800
ylim2 = 1020

# Look for max SSH in the BG regions:
# Exclude shallow regions/shelves
#iBG1 = 441
#jBG1 = 899
#iBG2 = 521
#jBG2 = 970
IIM = [430, 430, 490, 500, 500, 490, 430]
JJM = [905, 980, 980, 950, 930, 902, 902]

# Region where to find closed Contours:
# Look for max SSH in the BG regions:
iC1 = 420
jC1 = 895
iC2 = 545
jC2 = 995
IBG = [iC1,iC1,iC2,iC2]
JBG = [jC1,jC2,jC2,jC1]

pthoras = '/work/Dmitry.Dukhovskoy/data/ORAS5/SSH'
pthout = '/work/Dmitry.Dukhovskoy/anls_output/oras5/BG_anls'

nyrs = YRE+1-YRS
icc = 0
SSHMAX = np.zeros((12*nyrs))
GHMN   = np.zeros((12*nyrs)) 
GHMD   = np.zeros((12*nyrs))
GHUPRC = np.zeros((12*nyrs))
GHLPRC = np.zeros((12*nyrs))
BGAR   = np.zeros((12*nyrs))
TMANLS = np.zeros((12*nyrs))
for YR in range(YRS,YRE+1):
  for MM in range(1,13):
    timeS = time.time()
    flnm = f'sossheig_control_monthly_highres_2D_{YR}{MM:02d}_CONS_v0.1.nc'
    if YR > 2014:
      flnm = f'sossheig_control_monthly_highres_2D_{YR}{MM:02d}_OPER_v0.1.nc'
    dfl = os.path.join(pthoras,flnm)
    print(f'Opening {dfl}')
    varnm='sossheig'
    dset = xarray.open_dataset(dfl)
    ssh = dset[varnm].data.squeeze()

    dnmb0 = mtime.datenum([YR,MM,15])
    TMANLS[icc] = dnmb0

    if icc == 0:
      jdm, idm = ssh.shape
      LMSK = np.where(np.isnan(ssh),0,1)
      X, Y     = np.meshgrid(np.arange(idm), np.arange(jdm))
      MS, _, _ = mmisc.inpolygon_v2(X, Y, IBG, JBG)  # 
      JBS, IBS = np.where( (MS == 1) & (LMSK == 1) ) #exclude deeep regions
      MSKBS  = np.zeros((jdm,idm))
      MSKBS[JBS,IBS] = 1
      # Mask for finding max_ssh in the BG:
      MSH, _, _ = mmisc.inpolygon_v2(X, Y, IIM, JJM)  # 
      #JSSH, ISSH = np.where((MSH==1) & (LMSK == 1))
      MSK_SSH = np.where((MSH==1) & (LMSK == 1), 1, 0)

      lonh = dset['nav_lon'].data
      lath = dset['nav_lat'].data
      # Mask for contouring ssh:
      # MSC, _, _, = mmisc.inpolygon_v2(X,Y,IC,JC)
      DX, DY = mmom6.dx_dy(lonh, lath)
      Acell  = DX*DY*1.e-6  # km2
      Acell_BG = Acell[jC1:jC2,iC1:iC2]
            

    ssh_dmn, amn = moras.demean_ssh(ssh)
    # 
    # Find max ssh and center of BG:
    #A = ssh_dmn[JBS,IBS]

    # SSH gradient, BG area based on the last closed contour:
    # Use subset domain
    ABG = ssh_dmn[jC1:jC2,iC1:iC2]
    MSK_BG = MSK_SSH[jC1:jC2,iC1:iC2]
    #ssh_max = np.nanmax(ssh_dmn[JSSH,ISSH])
    #ssh_min = np.nanmin(ABG)
    grdH_mn, grdH_md, grdH_lprc, grdH_uprc, area_bg, ssh_max = moras.calc_grad_sizeBG(ABG, lonh, lath, \
                                                      Acell_BG, MSK_BG)    

    cff=1.e7
    print(f'Area BG = {area_bg:.1f} km2, max ssh={ssh_max:.3f}m, grad_mean={grdH_md*cff:.4f}x1e-7') 

    iyr = YR-YRS
    imo = MM-1
    SSHMAX[icc] = ssh_max
    GHMN[icc]   = grdH_mn
    GHMD[icc]   = grdH_md
    GHUPRC[icc] = grdH_uprc
    GHLPRC[icc] = grdH_lprc
    BGAR[icc]   = area_bg    

    icc += 1
    timeE = time.time()
    print(f'Elapsed time: {(timeE-timeS)*1./60.:.3f} min ')


# Construct time array: days since the reference day:
dnmb_ref = mtime.datenum([1900,1,1])
dv_ref = mtime.datevec(dnmb_ref)
TM_ref = TMANLS - dnmb_ref
if TM_ref[0] < 0:
  TM_ref[0] = 0.

# Construct data set:
dim1 = 'time'
dim2 = 'lat'
dim3 = 'lon'

da_hmax = xarray.DataArray(SSHMAX, dims=(dim1,),
            coords={dim1: TM_ref})
da_ghmn = xarray.DataArray(GHMN, dims=(dim1,),
            coords={dim1: TM_ref})
da_ghmd = xarray.DataArray(GHMD, dims=(dim1,),
            coords={dim1: TM_ref})
da_ghuprc = xarray.DataArray(GHUPRC, dims=(dim1,),
            coords={dim1: TM_ref})
da_ghlprc = xarray.DataArray(GHLPRC, dims=(dim1,),
            coords={dim1: TM_ref})
da_bga = xarray.DataArray(BGAR, dims=(dim1,),
            coords={dim1: TM_ref})

da_hmax.attrs.update({
    'units': 'meters',
    'long_name': 'Sea Surface Height Maximum',
    'description': 'Maximum sea surface height in the Beaufort Gyre',
})
da_ghmn.attrs.update({
    'units': 'm/m',
    'long_name': 'mean SSH gradient',
    'description': 'mean SSH grad in the Beaufort Gyre from max SSH to last closed contour',
})
da_ghmd.attrs.update({
    'units': 'm/m',
    'long_name': 'median SSH gradient',
    'description': 'median SSH grad in the Beaufort Gyre',
})
da_ghuprc.attrs.update({
    'units': 'm/m',
    'long_name': '90 percentile of grad SSH',
    'description': '90 percentile of the SSH gradient',
})
da_ghlprc.attrs.update({
    'units': 'm/m',
    'long_name': '10 percentile of grad SSH',
    'description': '10 percentile of the SSH gradient',
})
da_bga.attrs.update({
    'units': 'km2',
    'long_name': 'BG area',
    'description': 'BG area within the last closed SSH contour',
})

dset_bg = xarray.Dataset({
    'ssh_max': da_hmax,
    'gradh_mean': da_ghmn,
    'gradh_med': da_ghmd,
    'gradh_uperc': da_ghuprc,
    'gradh_lperc': da_ghlprc,
    'BG_area': da_bga,
})
# Add global attributes:
dset_bg.attrs.update({
  "info": "Beaufort Gyre characteristics derived from ORAS5 ssh fields",
  "code": "calc_grad_ssh.py"
})


if f_save:
  fout = f'BG_gradH_{YRS}-{YRE}.nc'
  dfout = os.path.join(pthout,fout)
  print(f'Saving ---> {dfout}')
  dset_bg.to_netcdf(
         dfout,
         format='NETCDF4',
         engine='netcdf4',
    )
       

f_plt = False
if f_plt:
  clrmp = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBuGn_r')
  rmin = -0.5
  rmax = 0.5

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

  img = ax1.pcolormesh(ssh_dmn, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.axis('scaled')
  ax1.set_xlim([xlim1,xlim2])
  ax1.set_ylim([ylim1,ylim2])
  ax1.plot([iC1,iC2],[jC1,jC1],'-',color=[0.,0.6,0.8])
  ax1.plot([iC1,iC2],[jC2,jC2],'-',color=[0.,0.6,0.8])
  ax1.plot([iC1,iC1],[jC1,jC2],'-',color=[0.,0.6,0.8])
  ax1.plot([iC2,iC2],[jC1,jC2],'-',color=[0.,0.6,0.8])

  ax1.plot(IIM,JJM,'-')

  sttl = f'SSH demeaned, ORAS5, {YR}/{MM:02d}'
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

  btx = 'calc_grad_ssh.py'
  bottom_text(btx) 

