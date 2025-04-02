"""
  Compute monthly mean characteristics of the cold pool in the Bering Sea
  area by T classes
  overall T and S
  
  monthly mean and StDev  bottom T/S derived in calc_mnthlyTSbtm.py

  The coldpool is defined as:
  bottom water < 2C
  deep water (> 250 m) is excluded

  usage: calc_coldpool_seasfcst.py --expt=3 --MMI=4 --YRS=1993 --YRE=2008
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
YRE    = YRS
MMI    = 4
nens   = 1    # ens # for ensemble runs - =1 for 3D ocean fields
expt_nmb = 3  # =2 - seas. f/casts no irelax, =3 - seas. f/casts with ice relax

if args.expt:
  expt_nmb = args.expt
if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI

expt_name = "seasonal_daily"
runname   = f"NEPphys_frcst_dailyOB-expt{expt_nmb:02d}"
expt_nmb0 = f"{expt_nmb:02d}"

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP'][expt_name]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt_name]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt_name]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape

DX, DY = mmom6.dx_dy(hlon, hlat)
Acell  = DX*DY

# Get indices of the polygon:
II = pthseas['ANLS_NEP']['poly_BerSea']['II']
JJ = pthseas['ANLS_NEP']['poly_BerSea']['JJ']
jdm, idm = HH.shape

X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MS, _, _ = mmisc.inpolygon_v2(X, Y, II, JJ)  # 
JBS, IBS = np.where( (MS == 1) & (HH >= -250) & (HH < 0) ) #exclude deeep regions
MSKBS  = np.zeros((jdm,idm))
MSKBS[JBS,IBS] = 1

TCLASS = [-2., -1., 0., 1., 2.]
nTC    = len(TCLASS)

nyrs = YRE-YRS+1
CP_AREA   = np.zeros((nyrs,12))
CP_TCLASS = np.zeros((nyrs,12,nTC-1))
CP_TSTAT  = np.zeros((nyrs,12,3)) 
CP_SSTAT  = np.zeros((nyrs,12,3)) 
icc = 0
iyr = -1

for YRI in range(YRS,YRE+1):
  iyr += 1
  pthanls = pthseas['MOM6_NEP'][expt_name]['pthanls'].format(expt_nmb=expt_nmb)
  dflt = os.path.join(pthanls,f'mnthly_Tbtm_expt{expt_nmb:02d}_{YRI}{MMI:02d}.pkl')
  dfls = os.path.join(pthanls,f'mnthly_Sbtm_expt{expt_nmb:02d}_{YRI}{MMI:02d}.pkl')

  print(f'Loading T bottom stat  <-- {dflt}')
  with open(dflt, 'rb') as fid:
    Tbtm, Tstd, TM = pickle.load(fid)

  print(f'Loading S bottom stat  <-- {dfls}')
  with open(dfls, 'rb') as fid:
    Sbtm, Sstd, TM = pickle.load(fid)

  Tcp = 2.  # T defining cold pool
  for imo in range(12):
    Tb2D = Tbtm[imo,:,:].squeeze()
    Sb2D = Sbtm[imo,:,:].squeeze()
    JBS, IBS = np.where( (MSKBS == 1) & (Tb2D < Tcp))
    CParea = np.sum(Acell[JBS,IBS])*1e-6 # km2
    CP_AREA[iyr,imo] = CParea
    Tmd = np.nanmedian(Tb2D[JBS,IBS])
    Tlprc = np.percentile(Tb2D[JBS,IBS],10)
    Tuprc = np.percentile(Tb2D[JBS,IBS],90)
    Smd = np.nanmedian(Sb2D[JBS,IBS])
    Slprc = np.percentile(Sb2D[JBS,IBS],10)
    Suprc = np.percentile(Sb2D[JBS,IBS],90)

    CP_TSTAT[iyr,imo,:3]=[Tmd,Tlprc,Tuprc]
    CP_SSTAT[iyr,imo,:3]=[Smd,Slprc,Suprc]

    for icl in range(1,nTC):
      t1 = TCLASS[icl-1]
      t2 = TCLASS[icl]

      JBS, IBS = np.where( (MSKBS == 1) & (Tb2D <= t2) & (Tb2D > t1))
      if len(JBS) > 0.:
        CP_area = np.sum(Acell[JBS,IBS])*1e-6 # km2
      else:
        CP_area = 0.

      CP_TCLASS[iyr,imo,icl-1] = CP_area



f_save = True
if f_save:
  pthdump = pthseas["MOM6_NEP"][expt_name]["pthdump"] 
  fldump  = f'{runname}_coldpool_stat_MI{MMI:02d}_{YRS}_{YRE}.pkl'
  drfldump = os.path.join(pthdump, fldump)
  print(f'Saving ---> {drfldump}')

  with open(drfldump, 'wb') as fid:
    pickle.dump([CP_AREA,CP_TSTAT,CP_SSTAT,CP_TCLASS],fid)


# Plot time series:
plt.ion()

btx = 'calc_coldpool_stat_seasfcst.py'

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



