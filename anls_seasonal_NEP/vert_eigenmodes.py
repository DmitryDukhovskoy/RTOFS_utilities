"""
  Compute and plot vertical eigenvectors for modes = 1,...  computed from 
  Strum-Louivelle problem for the vertical modes W
  Obtained from linearized equation for long waves using
  separation of variables

  Similar approached used for deriving the 1st barocl. radius which is 
  the 1st mode 

  Based on my code for the Rossby R

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
import importlib

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
import mod_mom6 as mmom6
import mod_utils_ob as mutob
importlib.reload(mutob)
import mod_anls_seas as manseas
import mod_solver as msolv

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
sctnm  = 'xsct_25N' 
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1    
nens   = 2
mode   = 2  # 1st mode 1st barocl R or 2nd mode
dnmbR  = mtime.datenum([1993,8,15])  # day to plot


#expt    = "seasonal_fcst"
#runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
#expt    = 'NEP_BGCphys_GOFS'
#runname = 'NEP_physics_GOFS-IC'
expt    = 'NEP_seasfcst_LZRESCALE'
runname = 'NEPphys_LZRESCALE_climOB_1993_04-e02'
dnmbS   = mtime.datenum([YRS,MOS,DDS])
dayrun  = dnmbR - dnmbS + 1 # day to plot:  model forecast day run - closest output will be plotted

dvR = mtime.datevec(dnmbR)
print(f'Expt: {expt} Run: {runname} Plot date: {dvR[0]}/{dvR[1]}/{dvR[2]}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

if expt == 'seasonal_fcst':
  pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(runname=runname)
else:
  dnmb0    = dnmbR
  dv0      = mtime.datevec(dnmb0)
  YR0, MM0, DD0 = dv0[:3]
  jday0    = int(mtime.date2jday([YR0,MM0,DD0]))
  pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YY=YR0, MM=MM0)
pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
ndav       = pthseas['MOM6_NEP'][expt]['ndav']  # # of days output averaged

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

Is  = pthseas['ANLS_NEP'][sctnm]['II']
Js  = pthseas['ANLS_NEP'][sctnm]['JJ']

IJ     = np.column_stack((Is, Js))
DX, DY = mmom6.dx_dy(hlon,hlat)
SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive', check_pole=False)
II     = SGMT.I_indx
JJ     = SGMT.J_indx
nLeg   = SGMT.Leg_number
#II, JJ = mmisc.xsect_indx(Is, Js)
hLsgm1 = SGMT.half_Lsgm1
hLsgm2 = SGMT.half_Lsgm2
XX     = hlon[JJ,II]
YY     = hlat[JJ,II]
Hbtm   = HH[JJ,II]
LSgm   = np.zeros((len(II)))  # total segment length = half1 + half2
for ik in range(len(II)):
   LSgm[ik] = hLsgm1[ik] + hLsgm2[ik]

# Find closest output:
if expt == 'NEP_BGCphys_GOFS':
  ocnfld = 'ocean'
else:
  ocnfld = 'oceanm'

if expt == 'seasonal_fcst':
  pthfcst = os.path.join(pthfcst,f'{ocnfld}_{dvR[0]}{dvR[1]:02d}')

YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthfcst, dnmbR, fld=ocnfld)
dv0  = mtime.datevec(dnmb0)
YR0, MM0, DD0 = dv0[:3]

flocn_name = pthseas['MOM6_NEP'][expt]['focname'].format(YR=YR0, jday=jday0)
dfmom6 = os.path.join(pthfcst, flocn_name)

# Averaging period:
dnmb_av1 = dnmb0 - np.floor(ndav/2)
#if dnmb_av1 < dnmbS: dnmb_av1=dnmbS
dnmb_av2 = dnmb_av1 + ndav-1

dset   = xarray.open_dataset(dfmom6)

ZM  = -dset['zl'].data
T2d = dset['potT'].data[0,:,JJ,II].squeeze()
T2d = np.transpose(T2d)
S2d = dset['salt'].data[0,:,JJ,II].squeeze()
S2d = np.transpose(S2d)
ZZ  = mmom6.zm2zz(ZM)

# For plotting - fill land/bottom and smooth 2D field:
#A2di = mmom6.fill_bottom(A2d, ZZ, Hbtm)

# Convert potT --> in situ
print('Converting potential T ----> in situ')
import mod_regmom as mregmom
npnts = T2d.shape[1]
for ii in range(npnts):
  T1d = T2d[:,ii].copy()
  if np.isnan(T1d[0]): continue
  S1d = S2d[:,ii].copy()
  Ts = mregmom.pot2insitu_1D(T1d, S1d, ZM, YY[ii])
  T2d[:,ii] = Ts

# Get N2 profiles:
# 1st N2 value is below surface, 
# Use BC (phi[0] = 0) to solve e/problem
# Z_phi[0] is not 0
N2, Z_phi = manseas.calc_N2(T2d, S2d, ZM, YY)

ZZphi = np.zeros((len(Z_phi)+1))  # add surface to match Phi[0] surface
ZZphi[1:] = Z_phi
# Find layer thicknesses between phi-pnts - where density is
# either layer mid-points or layer interfaces - depending where
# T/S points are
dZlr_phi = abs(np.diff(ZM))

# Solve Sturm-Liouville:
print(f'Solving Strum-Liouville, requested mode={mode}')
npnts  = len(II)
EVct   = np.zeros((T2d.shape))*np.nan
CPhase = np.zeros((npnts))*np.nan
Rbrcl  = np.zeros((npnts))*np.nan
for ik in range(npnts):
  Hb0 = Hbtm[ik]
  if Hb0 >=0: continue
  N2z = N2[:,ik].copy()

  # Filter N-point running mean
  N2zF = msolv.runmn_dz(N2z, dZlr_phi, mnwnd=5)

  # NaNs for bottom values:
  # kbtm - near-bottom value for phi(z)
  if abs(Hb0) < abs(Z_phi[-1]): 
    D = Z_phi - Hb0
    kbtm = np.where(D <= 0.)[0][0] - 1
    N2z[kbtm+1:]  = np.nan
    N2zF[kbtm+1:] = np.nan
  else:
# Truncate profile to the deepest level, make a bottom
# to keep Phi correct size
    kbtm = len(Z_phi)-1
    N2z[-1] = np.nan
    N2zF[-1] = np.nan

# Skip shallow regions
# for small kbtm - problem creating a matrix and solving it
  if kbtm < 5:
    print(f'Shallow location: Hbtm={Hb0:.2f} m')
    continue

  Phi, Rrsb, Cphs = manseas.solve_SturmLiouville(N2zF, Hb0, Z_phi, \
                           YY[ik], mode=mode)
# For modes > 1:
# Make same +/- pattern with + in the upper ocean  
  if mode > 1:
    ineg = np.where(Phi < -1.e-23)[0][0]
    if ineg > 2: Phi = -Phi

  nphi = len(Phi)
  EVct[:nphi,ik] = Phi
  CPhase[ik] = Cphs
  Rbrcl[ik]  = Rrsb


import mod_colormaps as mclrmp

if mode == 1:
# Make sure all vectors are of the same sign:
  EVct = np.where(EVct > 0., -EVct, EVct)
  rmin, rmax = mclrmp.minmax_clrmap(EVct)

  clrmp_name = 'YlGnBu_r'
  clr_ramp   = [1, 1, 1]   # add white color at the end of the colormap
  clrmp = mclrmp.addendclr_colormap(clrmp_name, clr_ramp, nramp=0.1, ramp_start=False)
  clrmp.set_bad(color=[0.1, 0.1, 0.1])
elif mode > 1:
  rmin, rmax = mclrmp.minmax_clrmap(EVct, pmin=5, pmax=95, fsym=True)
  clrmp = mclrmp.colormap_ssh(cpos='Oranges', cneg='Blues_r', nclrs=200)
  clrmp.set_bad(color=[1,1,1])

# Indices:
#Xdist = np.arange(len(Hbtm))
# 
# Distance along the section
# normalize by the total distance
Lsection = mmisc.dist_sphcrd(YY[-1],XX[-1],YY[0],XX[0]) # total length of the section, m
Xdist = np.cumsum(LSgm)
Xdist = Xdist-Xdist[0]
Xdist = Xdist/Xdist[-1]*Lsection*1.e-3  # normalized, km

xl1  = min(Xdist)
xl2  = max(Xdist)

dv_av1 = mtime.datevec(dnmb_av1)
yrs, mms, dds = dv_av1[:3]
dv_av2 = mtime.datevec(dnmb_av2)
yre, mme, dde = dv_av2[:3]

lon1 = XX[0]
lat1 = YY[0]
lon2 = XX[-1]
lat2 = YY[-1]
stxt = f"X-axis: Distance (km) along section from {lon1:5.2f}W, {lat1:5.2f}N to {lon2:5.2f}W, {lat2:5.2f}N"

sttl = f"{runname} barocl eigenvector {mode}, avrg: {yrs}/{mms}/{dds}-{yre}/{mme}/{dde} {sctnm}"
btx = 'vert_eigenmodes.py'

if mode == 1:
  mxsct.plot_xsection(EVct, Xdist, ZZphi, Hbtm, Xdist, clrmp, rmin, rmax, \
                    xl1, xl2, sttl=sttl, stxt=stxt, \
                    plot_map=HH, I_indx=II, J_indx=JJ, btx=btx)
else:
  mxsct.plot_xsection(EVct, Xdist, ZZphi, Hbtm, Xdist, clrmp, rmin, rmax, \
                    xl1, xl2, sttl=sttl, stxt=stxt, \
                    plot_map=HH, I_indx=II, J_indx=JJ, btx=btx, contours=[0.])


#bottom_text(btx)

f_chck = False
if f_chck:
  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])

  ax1.contour(HH,[0], colors=[(0,0,0)])
  ax1.contour(HH,[-1000], linestyles='solid', colors=[(0,0.6,0.9)])
  ax1.contour(hlat,[x for x in range(10,88,5)], linestyles='solid', colors=[(0.8, 0.8, 0.8)])
  ax1.contour(hlat, [55.125], linestyles='solid', colors=[(1., 0.4, 0.)])
  ax1.axis('scaled')

# Bind the button_press_event with the onclick() method
  fig1.canvas.mpl_connect('button_press_event', onclick)


