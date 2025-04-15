"""
  Plot monthly mean Bering Strait v.sections of T&S
  Extracted in calc_BerStrFlux.py 
  Saved indiviual years in separate files

  usage: run BerStr_UV.py --expt 3 --MMI 7 --YRS 1993 --YRE 2008 --MMS 1 --MME 3

  MMS can be > MME, e.g. for winter average: MMS=12, MME=2 to include Dec - Feb. 
  However, for MMI=1, Jan, Feb and Dec are not in cnosequtive order being 11 months apart
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import pickle
import xarray
from yaml import safe_load
import argparse
from matplotlib.patches import Polygon

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
parser.add_argument("--MMS", help="start calendar month for averaging, 1,...,12", type=int)
parser.add_argument("--MME", help="end calendar month for averaging, >=MMS 1, ..., 12", type=int)
args = parser.parse_args()

# Default values: 
YRS    = 1993 # year start of the forecast
YRE    = 2008
MMI    = 4
MMS    = 1
MME    = MMS   # one month
nens   = 1    # ens # for ensemble runs - =1 for 3D ocean fields
expt_nmb = 3  # =2 - seas. f/casts no irelax, =3 - seas. f/casts with ice relax
Tref   = 0.   # Ref t for computing heat flux
Sref   = 34.8 # Ref S for FW flux

if args.expt:
  expt_nmb = args.expt
if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI
if args.MMS:
  MMS = args.MMS
  MME = MMS
if args.MME:
  MME = args.MME

# Flags:
plot_obs = True # add obs.-based estimates from Woodgate 2018

expt_name = "seasonal_daily"
runname   = f"NEPphys_frcst_dailyOB-expt{expt_nmb:02d}"

if MMI==1 and YRS==1993:
  print(f"first start year for init month {MMI} is 1994, changing {YRS} to 1994")
  YRS = 1994

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


# Define Bering Strait:
IIb = [201, 201]
JJb = [680, 667]   # <--- this will result in positive flux to the Arctic Oc. 
#JJb = [667, 680]  # <--- this will result in positive flux to the Ber. Sea

# Define segment lengths, norms, etc:
# positive: norm vector is to the left as follow the section
# negative: to the right
import mod_mom6_valid as mom6vld
DX, DY = mmom6.dx_dy(hlon,hlat)
IJ     = np.column_stack((IIb, JJb))
SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive', check_pole=False)
II     = SGMT.I_indx
JJ     = SGMT.J_indx

# Arrange time series in calendar months from f/cast months:
Tfcst = np.arange(MMI,MMI+12)
Tcal = np.mod(Tfcst,12)
Tcal = np.where(Tcal==0,12,Tcal)

if MMS <= MME:
  mav = np.arange(MMS,MME+1)  # months for seasonal averaging
else:
  mav = np.arange(1,MME+1)
  mav = np.append(mav,np.arange(MMS,13))
  
#iorder = Tcal
#ip = np.where(Tcal==12)[0][0]+1
#if ip < 12:
#  icol = np.arange(12)
#  iorder = icol+ip
#  iorder = np.where(iorder>=12, iorder-12, iorder)

Ssum = None
Tsum = None
pthanls = pthseas['MOM6_NEP'][expt_name]['pthanls'].format(expt_nmb=expt_nmb)
nyrs = YRE-YRS+1
icc = 0
for YRI in range(YRS,YRE+1):
  dflnm = os.path.join(pthanls,f'mnthly_BerSea_Fluxes_expt{expt_nmb:02d}_{YRI}{MMI:02d}.pkl')
  # Saved fields: VFlx,FWFlx,HFlx,UV2d,TT,SS,ZM,ZZ,Hbtm,LSgm,XX,YY
  print(f'Loading Fluxes {dflnm}')
  with open(dflnm, 'rb') as fid:
    _,_,_,_,TT,SS,ZM,ZZ,Hbtm,LSgm,XX,YY = pickle.load(fid)

  #U2d = np.expand_dims(U2d, axis=0)
  # Select needed months
  for imm in range(len(mav)):
    mo = mav[imm]
    jmo = np.where(Tcal == mo)[0][0]
    if Tsum is None:
      Tsum = TT[jmo,:,:]
      Tsum = np.expand_dims(Tsum, axis=0)
      Ssum = SS[jmo,:,:]
      Ssum = np.expand_dims(Ssum, axis=0)
    else:
      dmm =  TT[jmo,:,:]
      dmm = np.expand_dims(dmm, axis=0)
      Tsum = np.append(Tsum, dmm, axis=0) 
      dmm =  SS[jmo,:,:]
      dmm = np.expand_dims(dmm, axis=0)
      Ssum = np.append(Ssum, dmm, axis=0) 

  icc += 1


Tmn = np.mean(Tsum, axis=0)
Tstd = np.std(Tsum, axis=0)
nrec = Tsum.shape[0]
Smn = np.mean(Ssum, axis=0)
Sstd = np.std(Ssum, axis=0)

Sstd[24:,:] = np.nan
Tstd[24:,:] = np.nan

# Plot vertical section of U northward:
from matplotlib.patches import Polygon
plt.ion()

btx = 'BerStr_TS.py'


CLRS = [[0.6, 0.02, 0.6],
        [0.2, 0.38, 1],
        [0., 0.8, 0.5],
        [0.2,1.,0.8],
        [0.6, 0.2, 0],
        [0.9, 0.8, 0.],
        [1, 0.98, 0.9]]

tcmp = mclrmps.colormap_posneg_uneven(CLRS)
tmin = -2.
tmax = 5.
scmp = mclrmps.colormap_haline2()
smin = 30.
smax = 33.


tcntrs = [x/10 for x in range(5,40,5)]
scntrs = [x/100 for x in range(5,70,5)]

# Patch bottom:
verts = [(np.max(XX),-6000),*zip(XX,Hbtm),(np.min(XX),-6000)]
poly = Polygon(verts, facecolor='0.3', edgecolor='0.3', zorder=5)
polyS = Polygon(verts, facecolor='0.3', edgecolor='0.3', zorder=5)

sttl = f'Bering Str T & StDev\n {runname} init M={MMI:02d} avrg: {YRS}-{YRE}, {MMS}/{MME}'
sttl2 = 'mean S & stdev'

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.55, 0.5, 0.4])
img = ax1.pcolormesh(XX,ZM, Tmn, cmap=tcmp, vmin=tmin, vmax=tmax)
ax1.set_ylim([-60,0])
ax1.set_xlim([190.2, 192.1])
CS = ax1.contour(XX,ZM,Tstd,tcntrs, linestyles='solid', colors=[(0.,0.,0)])
ax1.clabel(CS, inline=True, fontsize=10)
ax1.add_patch(poly)
ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(tmin,tmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ax3 = plt.axes([0.1, 0.08, 0.5, 0.4])
img2 = ax3.pcolormesh(XX,ZM,Smn, cmap=scmp, vmin=smin, vmax=smax)
ax3.set_ylim([-60,0])
ax3.set_xlim([190.2, 192.1])
CS2 = ax3.contour(XX,ZM,Sstd,scntrs, linestyles='solid', colors=[(0.,0.,0)])
ax3.clabel(CS2, inline=True, fontsize=10)
ax3.add_patch(polyS)
ax3.set_title(sttl2)

ax4 = fig1.add_axes([ax3.get_position().x1+0.025, ax3.get_position().y0,
                   0.02, ax3.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img2, cax=ax4, orientation='vertical', extend='both')
ax4.yaxis.set_ticks(list(np.linspace(smin,smax,11)))
ax4.set_yticklabels(ax4.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

bottom_text(btx, pos=[0.05, 0.02])






