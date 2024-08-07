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
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_mom6 as mmom6

# experiment: year start, month start, ...
varnm = 'temp'
YRS   = 1993 
MOS   = 4
nens  = 2
expt  = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
sctnm = 'xsct_NOB' 
jday  = 118
dnmb0 = mtime.jday2dnmb(1993,jday)
ndav  = 5  # # of days averaged
dnmbS = dnmb0-ndav+1


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

pthfcst    = pthseas['MOM6_NEP']['seasonal_fcst']['pthoutp'].format(expt=expt)
pthtopo    = gridfls['MOM6_NEP']['seasonal_fcst']['pthgrid']
fgrid      = gridfls['MOM6_NEP']['seasonal_fcst']['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"]["seasonal_fcst"]["ftopo"]
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

Is  = pthseas['ANLS_NEP'][sctnm]['II']
Js  = pthseas['ANLS_NEP'][sctnm]['JJ']

nlegs   = len(Is) - 1
IJ      = np.zeros((nlegs+1,2))
IJ[:,0] = Is
IJ[:,1] = Js


JJ, II = mmisc.xsect_indx(Is, Js)

DX, DY = mmom6.dx_dy(hlon,hlat)
SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive', check_pole=False)
II     = SGMT.I_indx
JJ     = SGMT.J_indx
nLeg   = SGMT.Leg_number
XX     = hlon[JJ,II]
YY     = hlat[JJ,II]
Hbtm   = HH[JJ,II]
LSgm   = np.zeros((len(II)))  # total segment length = half1 + half2
for ik in range(len(II)):
   LSgm[ik] = hLsgm1[ik] + hLsgm2[ik]

DV     = mtime.datevec(dnmb0)
YR     = DV[0]
MM     = DV[1]
DD     = DV[2]
dfmom6 = os.path.join(pthfcst,f'oceanm_{YRS}{MM:02d}',f'oceanm_{YRS}_{jday:03d}.nc')
dset   = xarray.open_dataset(dfmom6)

XX  = hlon[JJ,II]
YY  = hlat[JJ,II]
ZM  = -dset['zl'].data
if varnm == 'temp':
  A2d = dset['potT'].data[0,:,JJ,II].squeeze()
A2d = np.transpose(A2d)
ZZ  = mmom6.zm2zz(ZM)

# For plotting - fill land/bottom and smooth 2D field:
A2di = mmom6.fill_bottom(A2d, ZZ, Hbtm)

if varnm == 'salin':
  clrmp = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = pthseas['ANLS_NEP'][sctnm]['smin']
  rmax = pthseas['ANLS_NEP'][sctnm]['smax']
elif varnm == 'temp':
  clrmp = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[1,1,1])
  rmin = pthseas['ANLS_NEP'][sctnm]['tmin']
  rmax = pthseas['ANLS_NEP'][sctnm]['tmax']

Xdist = np.arange(len(Hbtm))

xl1  = min(Xdist)
xl2  = max(Xdist)

dvS = mtime.datevec(dnmbS)
yrs, mms, dds = dvS[:3]

stxt = 'grid index'
sttl = f"{expt} {varnm} avrg: {yrs}/{mms}/{dds} - {YR}/{MM}/{DD} {sctnm}"
mutob.plot_xsection(A2d, Xdist, ZM, Hbtm, Xdist, clrmp, rmin, rmax, xl1, xl2, sttl=sttl, stxt=stxt, fgnmb=1)

btx = 'plot_xsectTS.py'
bottom_text(btx)

f_chck = False
if f_chck:
  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])

  ax1.contour(HH,[0], colors=[(0,0,0)])
  ax1.axis('scaled')

