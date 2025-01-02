"""
  Plot sections with vertical distribution of T or S
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
importlib.reload(mutob)

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
varnm  = 'salin'  # temp (potential) / salin
f_insitu = True    # convert to in situ 
#sctnm  = 'xsct_BerSea' 
sctnm = 'xsct_CalUnderCur_358N'
# Start of the run - needed only for seasonal forecasts:
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1    
nens   = 2    # ens # for ensemble runs
dnmbR  = mtime.datenum([1993,11,15])  # day to plot

if varnm == 'salin': f_insitu = False

#expt    = "seasonal_fcst"
#runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
expt    = 'NEP_BGCphys_GOFS'
runname = 'NEP_physics_GOFS-IC'
#expt    = 'NEP_seasfcst_LZRESCALE'
#runname = 'NEPphys_LZRESCALE_climOB_1993_04-e02'
dnmbS   = mtime.datenum([YRS,MOS,DDS]) 
#dnmbR   = dnmbS + dayrun - 1
dayrun  = dnmbR - dnmbS + 1 # day to plot:  model forecast day run - closest output will be plotted

dvR = mtime.datevec(dnmbR)
print(f'Expt: {expt} Run: {runname} Plot date: {dvR[0]}/{dvR[1]}/{dvR[2]}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

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

# Function to print mouse click event coordinates
def onclick(event):
   print([event.xdata, event.ydata])

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
if varnm == 'temp':
  A2d = dset['potT'].data[0,:,JJ,II].squeeze()
  if f_insitu:
    S2d = dset['salt'].data[0,:,JJ,II].squeeze()
    S2d = np.transpose(S2d)
else:
  A2d = dset['salt'].data[0,:,JJ,II].squeeze()
A2d = np.transpose(A2d)
ZZ  = mmom6.zm2zz(ZM)

# For plotting - fill land/bottom and smooth 2D field:
#A2di = mmom6.fill_bottom(A2d, ZZ, Hbtm)

# Convert to in situ to compare with obs. data
if f_insitu:
  import mod_regmom as mregmom
  Tsitu = A2d.copy()
  npnts = Tsitu.shape[1]
  for ii in range(npnts):
    T1d = A2d[:,ii].copy()
    if np.isnan(T1d[0]): continue
 
    S1d = S2d[:,ii].copy()
#    Tp = mregmom.insitu2pot_1D(T1d, S1d, ZM, YY[ii], printT=True)
    Ts = mregmom.pot2insitu_1D(T1d, S1d, ZM, YY[ii])
    Tsitu[:,ii] = Ts

  A2d = Tsitu.copy()
  

if varnm == 'salin':
#  clrmp = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  clrmp = mclrmps.colormap_salin2()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = pthseas['ANLS_NEP'][sctnm]['smin']
  rmax = pthseas['ANLS_NEP'][sctnm]['smax']
elif varnm == 'temp':
  clrmp = mclrmps.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[1,1,1])
  rmin = pthseas['ANLS_NEP'][sctnm]['tmin']
  rmax = pthseas['ANLS_NEP'][sctnm]['tmax']
  if sctnm == 'xsct_BerSea' and (MM0>6 and MM0<=12):
    rmax = 10.

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

if f_insitu and varnm=='temp':
  sttl = f"{runname} Tinsitu avrg: {yrs}/{mms}/{dds}-{yre}/{mme}/{dde} {sctnm}"
else:
  sttl = f"{runname} {varnm} avrg: {yrs}/{mms}/{dds}-{yre}/{mme}/{dde} {sctnm}"

btx = 'plot_xsectTS.py'
mxsct.plot_xsection(A2d, Xdist, ZM, Hbtm, Xdist, clrmp, rmin, rmax, \
                    xl1, xl2, sttl=sttl, stxt=stxt, fgnmb=1, \
                    plot_map=HH, I_indx=II, J_indx=JJ, btx=btx)

#bottom_text(btx)

f_chck = False
if f_chck:
  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])

  ax1.contour(HH,[0], colors=[(0,0,0)])
  ax1.contour(HH,[-1000], linestyles='solid', colors=[(0,0.6,0.9)])
  ax1.contour(HH,[-500], linestyles='solid', colors=[(1,0.6,0.)])
  ax1.contour(hlat,[x for x in range(10,88,5)], linestyles='solid', colors=[(0.8, 0.8, 0.8)])
  ax1.contour(hlat, [55.125], linestyles='solid', colors=[(1., 0.4, 0.)])
  ax1.axis('scaled')

# Bind the button_press_event with the onclick() method
  fig1.canvas.mpl_connect('button_press_event', onclick)


