"""
  Plot monthly mean Bering Strait heat / FW fluxes
  Extracted in calc_BerStrFlux.py 
  Saved indiviual years in separate files

  usage: run plot_BerStrFlux.py --expt=3 --MMI=7 --YRS=1993 --YRE=2008
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

# Default values: 
YRS    = 1993 # year start of the forecast
YRE    = 2008
MMI    = 4
nens   = 1    # ens # for ensemble runs - =1 for 3D ocean fields
expt_nmb = 3  # =2 - seas. f/casts no irelax, =3 - seas. f/casts with ice relax
Tref   = -1.9   # Ref t for computing heat flux
Sref   = 34.8 # Ref S for FW flux

if args.expt:
  expt_nmb = args.expt
if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI

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

DX, DY = mmom6.dx_dy(hlon, hlat)
#Acell  = DX*DY

# Define Bering Strait:
IIb = [201, 201]
JJb = [680, 667]   # <--- this will result in positive flux to the Arctic Oc. 
#JJb = [667, 680]  # <--- this will result in positive flux to the Ber. Sea

# Define segment lengths, norms, etc:
# positive: norm vector is to the left as follow the section
# negative: to the right
import mod_mom6_valid as mom6vld
IJ     = np.column_stack((IIb, JJb))
SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive', check_pole=False)
II     = SGMT.I_indx
JJ     = SGMT.J_indx
nLeg   = SGMT.Leg_number
LegNorm  = SGMT.LegNorm
#II, JJ = mmisc.xsect_indx(Is, Js)
hLsgm1 = SGMT.half_Lsgm1
hLsgm2 = SGMT.half_Lsgm2
Vnrm1  = SGMT.half_norm1
Vnrm2  = SGMT.half_norm2
XX     = hlon[JJ,II]
YY     = hlat[JJ,II]
Hbtm   = HH[JJ,II]
LSgm   = np.zeros((len(II)))  # total segment length = half1 + half2

# Derive layer thicknesses:
#import mod_mom6 as mom6util
#ocnfld = 'oceanm'
#pthoutp0 = pthseas['MOM6_NEP'][expt_name]['pthoutp'].format(expt_nmb=expt_nmb)
#subdir=f'oceanm_{YRS}{MMI:02d}'
#pthfcst0 = os.path.join(pthoutp0,f'{YRS}-{MMI:02d}-e01','history')
#list_files = manseas.list_oceanice_files(pthfcst0, prefix=ocnfld, subdir=subdir)
#pthfull = os.path.join(pthfcst0,subdir)
#floceanm = list_files[0]
#ZM = manseas.read_oceanm3D_field(pthfull, floceanm, 'zl', notime=False)
#ZM = -abs(ZM)
#nlrs = len(ZM)
#
#ZZ = mom6util.zm2zz(ZM)
#dZ = abs(np.diff(ZZ))

if plot_obs:
  fobsyaml = 'bering_fluxes.yaml'
  with open(fobsyaml) as fy:
    flxobs = safe_load(fy)

  VFobs_mn = np.array(flxobs["bering_fluxes"]["Vol"]["mean"])
  VFobs_er = np.array(flxobs["bering_fluxes"]["Vol"]["mean_err"])

  HFlobs_mn = np.array(flxobs["bering_fluxes"]["Heat_low"]["mean"])
  HFlobs_er = np.array(flxobs["bering_fluxes"]["Heat_low"]["mean_err"])

  HFhobs_mn = np.array(flxobs["bering_fluxes"]["Heat_high"]["mean"])
  HFhobs_er = np.array(flxobs["bering_fluxes"]["Heat_high"]["mean_err"])

  FWlobs_mn = np.array(flxobs["bering_fluxes"]["FW_low"]["mean"])
  FWlobs_er = np.array(flxobs["bering_fluxes"]["FW_low"]["mean_err"])

  FWhobs_mn = np.array(flxobs["bering_fluxes"]["FW_high"]["mean"])
  FWhobs_er = np.array(flxobs["bering_fluxes"]["FW_high"]["mean_err"])

Xobs = np.arange(1,13)

pthanls = pthseas['MOM6_NEP'][expt_name]['pthanls'].format(expt_nmb=expt_nmb)
nyrs = YRE-YRS+1
VFlx  = np.zeros((nyrs,12))
FWFlx = np.zeros((nyrs,12))
HFlx  = np.zeros((nyrs,12))
icc = 0
for YRI in range(YRS,YRE+1):
  dflnm = os.path.join(pthanls,f'mnthly_BerSea_Fluxes_expt{expt_nmb:02d}_{YRI}{MMI:02d}.pkl')
  # Saved fields: VFlx,FWFlx,HFlx,UV2d,TT,SS,ZM,ZZ,Hbtm,LSgm,XX,YY
  print(f'Loading Fluxes {dflnm}')
  with open(dflnm, 'rb') as fid:
    vf,fw,hf,_,_,_,_,_,_,_,_,_ = pickle.load(fid)

  VFlx[icc,:]  = vf*1.e-6    # Sv, >0 - to the AO
  FWFlx[icc,:] = fw*1.e-3    # mSv  >0 - to the AO
  HFlx[icc,:]  = -hf*1.e-12   # TJ/s or TW, -1 - make Northward (to the AO) heat flux >0  

  icc += 1

# Arrange time series in calendar months from f/cast months:
Tfcst = np.arange(MMI,MMI+12)
Tcal = np.mod(Tfcst,12)
Tcal = np.where(Tcal==0,12,Tcal)

ip = np.where(Tcal==12)[0][0]+1
if ip < 12:
  icol = np.arange(12)
  iorder = icol+ip
  iorder = np.where(iorder>=12, iorder-12, iorder)
  VFlx  = VFlx[:,iorder]
  FWFlx = FWFlx[:,iorder]
  HFlx  = HFlx[:,iorder]

def get_stat(AA):
  Amn = np.mean(AA, axis=0)
  Amin = np.min(AA, axis=0)
  Amax = np.max(AA, axis=0)

  return Amn, Amin, Amax

VFmn, VFmin, VFmax = get_stat(VFlx)
FWmn, FWmin, FWmax = get_stat(FWFlx)
HFmn, HFmin, HFmax = get_stat(HFlx)

VF_mn = np.mean(VFmn)
FW_mn = np.mean(FWmn)
HF_mn = np.mean(HFmn)

def plot_stat(Amn, Amin, Amax, ax1, clrl, Xtime, clrb=[]):
  if len(clrb) == 0:
    clrb=clrl
  ax1.plot(Xtime, Amn, '-', linewidth=2, color=clrl)
  ax1.plot(Xtime, Amn, 'o', markersize=8, markerfacecolor=clrl, markeredgecolor='none')
  dtg = 0.05
  for ik in range(len(Xtime)):
    xm0 = Xtime[ik]
    min1 = Amin[ik]
    max1 = Amax[ik]
    ax1.plot([xm0, xm0],[min1,max1],'-',linewidth=1, color=clrb)
    ax1.plot([xm0-dtg,xm0+dtg],[max1,max1],'-',linewidth=1, color=clrb)
    ax1.plot([xm0-dtg,xm0+dtg],[min1,min1],'-',linewidth=1, color=clrb)
    ax1.set_xticks(Xtime)
    ax1.set_xlim([0.5, 12.5])
    ax1.grid('on')
  return ax1

def plot_observations(ax1,ymn,yer,clr_obs,clr_err):
  """
    Add observational estimates to the existing figure
  """
  ylow = ymn-yer
  yup  = ymn+yer
  verts = [*zip(Xobs,yup),*zip(np.flip(Xobs),np.flip(ylow))]
  ax1.plot(Xobs,ymn,'-', linewidth=1, color=clr_obs)
  ax1.plot(Xobs,ymn,'o', markersize=6, markerfacecolor=clr_obs, markeredgecolor='none')
  poly = Polygon(verts, facecolor=clr_err, zorder=1)
  ax1.add_patch(poly)

  return ax1

def minmax_plot(amin, amax):
  ylim1 = min(amin)
  ylim2 = max(amax)
  dy = 0.05*abs(ylim2-ylim1)
  ylim1 = ylim1-dy
  ylim2 = ylim2+dy

  return ylim1, ylim2

# Plot stat:
from matplotlib.patches import Polygon
plt.ion()

btx = 'plot_BerStrFlux.py'

Xtime = np.arange(1,13)

sttl = f'{runname} init M={MMI:02d} {YRS}-{YRE}, Bering Str monthly Northward fluxes\n' \
       f' Vol Transport, Mean={VF_mn:.2f} Sv'
clr_obs=[0,0,0]
clr_err=[0.9,0.9,0.9]

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.7, 0.8, 0.24])
clrl = [0.8, 0., 0.4]
clrb = [0, 0, 0]
ax1 = plot_stat(VFmn, VFmin, VFmax, ax1, clrl, Xtime) 
ylim1, ylim2 = minmax_plot(VFmin, VFmax)
if plot_obs:
  ymn = VFobs_mn
  yer = VFobs_er
  ylow = ymn-yer
  yup  = ymn+yer
  ax1 = plot_observations(ax1,ymn,yer,clr_obs,clr_err)
  ylim1, ylim2 = minmax_plot([ylim1,np.min(ylow)],[ylim2,np.max(yup)])
  sttl = sttl + f'  Obs: {np.mean(ymn):.2f} Sv'
ax1.set_ylim([ylim1,ylim2])
ax1.set_title(sttl)

ax2 = plt.axes([0.1, 0.4, 0.8, 0.24])
clrfw = [0.0, 0.5, 0.8]
ax2 = plot_stat(FWmn, FWmin, FWmax, ax2, clrfw, Xtime) 
sttl2 = f'FW Flux, Sref={Sref:.2f}, Mean={FW_mn:.2f} mSv'
ylim1, ylim2 = minmax_plot(FWmin, FWmax)
if plot_obs:
  ymn = FWhobs_mn
  yer = FWhobs_er
  ylow = ymn-yer
  yup  = ymn+yer
  ax2 = plot_observations(ax2,ymn,yer,clr_obs,clr_err)
  ylim1, ylim2 = minmax_plot([ylim1,np.min(ylow)],[ylim2,np.max(yup)])
  ytot1 = np.mean(ymn)
# Low estimates:
  ymn = FWlobs_mn
  yer = FWlobs_er
  ylow = ymn-yer
  yup  = ymn+yer
  ax2 = plot_observations(ax2,ymn,yer,clr_obs,clr_err)
  ylim1, ylim2 = minmax_plot([ylim1,np.min(ylow)],[ylim2,np.max(yup)])
  ytot2 = np.min(ymn)
  sttl2 = sttl2 + f'  Obs: {ytot1:.2f}/{ytot2:.2f} mSv'
ax2.set_ylim([ylim1,ylim2])
ax2.set_title(sttl2)

ax3 = plt.axes([0.1, 0.1, 0.8, 0.24])
clrhf = [0.9, 0.5, 0.0]
ax3 = plot_stat(HFmn, HFmin, HFmax, ax3, clrhf, Xtime) 
sttl3 = f'Heat Flux, Tref={Tref:.2f}C, Mean={HF_mn:.2f} TW'
ylim1, ylim2 = minmax_plot(HFmin, HFmax)
if plot_obs:
  ymn = HFhobs_mn
  yer = HFlobs_er
  ylow = ymn-yer
  yup  = ymn+yer
  ax3 = plot_observations(ax3,ymn,yer,clr_obs,clr_err)
  ylim1, ylim2 = minmax_plot([ylim1,np.min(ylow)],[ylim2,np.max(yup)])
  ytot1 = np.mean(ymn)
# Low estimates
  ymn = HFlobs_mn
  yer = HFlobs_er
  ylow = ymn-yer
  yup  = ymn+yer
  ax3 = plot_observations(ax3,ymn,yer,clr_obs,clr_err)
  ylim1, ylim2 = minmax_plot([ylim1,np.min(ylow)],[ylim2,np.max(yup)])
  ytot2 = np.mean(ymn)
  sttl3 = sttl3 + f'  Obs: {ytot1:.2f}/{ytot2:.2f} TW'
ax3.set_ylim([ylim1,ylim2])
ax3.set_title(sttl3)

if plot_obs:
  ax4 = plt.axes([0.6, 0.02, 0.38, 0.06])
  xp=0.1
  yp=0.4
  ax4.plot([0,xp],[yp,yp],'-',linewidth=1, color=clr_obs)
  ax4.plot(0.5*xp,yp,'o',markersize=6, markerfacecolor=clr_obs, markeredgecolor='none')
  ax4.text(0.12,0.08,'Moor.Obs. 2003-2015, Woodgate 2018\n'+'FW and Heat Fluxes - low/high estimates')
  ax4.set_xlim([0,1])
  ax4.set_ylim([0,1])
  ax4.axis('off')

bottom_text(btx, pos=[0.02, 0.025])


