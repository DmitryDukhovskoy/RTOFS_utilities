"""
  Plot sections with vertical distribution of U normal to section
  Quick plots - just plot U or V, works fine for sections along I /J axis
  for slanted sections - need to project onto overall norm vector

  Also note that symmetric grid in MOM6 has extra J point for V and I point for U 
  components, so that V[0] is on the S boundary, and V[dim2+1] on the N. boundary
  similar for U
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
import mod_colormaps as mclrmp
#import mod_valid_utils as mvutil
import mod_mom6 as mmom6
import mod_utils_ob as mutob
importlib.reload(mutob)

# experiment: year start, month start, ...
varnm  = 'Unorm'
sctnm  = 'xsct_WOB' 
YRS    = 1993 # year start
MOS    = 4
DDS    = 1    
nens   = 2
dayrun = 21 # day to plot:  model forecast day run - closest output will be plotted
ndav   = 5  # # of days output averaged

expt  = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
dnmbS = mtime.datenum([YRS,MOS,DDS]) 
dnmbR = dnmbS + dayrun - 1

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

# Find closest output:
pth_day1=os.path.join(pthfcst,f'oceanm_{YRS}{MOS:02d}') # dir with 1st output
yr0, mo0, md0, jday0 = mutob.find_NEPoutput_day(pth_day1, dnmbR, ndav, fld='oceanm')
dnmb0 = mtime.jday2dnmb(yr0,jday0)
dnmb_av1 = dnmb0 - ndav + 1
if dnmb_av1 < dnmbS: dnmb_av1=dnmbS

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

DV     = mtime.datevec(dnmb0)
YR     = DV[0]
MM     = DV[1]
DD     = DV[2]
dfmom6 = os.path.join(pthfcst,f'oceanm_{YRS}{MM:02d}',f'oceanm_{YRS}_{jday0:03d}.nc')
dset   = xarray.open_dataset(dfmom6)

ZM  = -dset['zl'].data
U2d  = dset['u'].data[0,:,JJ,II].squeeze()
U2d  = np.transpose(U2d)
V2d  = dset['v'].data[0,:,JJ,II].squeeze()
V2d  = np.transpose(V2d)
#DP2d = dset['h'].data[0,:,JJ,II].squeeze() # layer thickness
#DP2d = np.transpose(DP2d)
ZZ  = mmom6.zm2zz(ZM)

# Construct Unorm based on local norm vector:
Unrm = U2d.copy()
Unrm[:,0] = abs(II[1]-II[0])*V2d[:,0] + abs(JJ[1]-JJ[0])*U2d[:,0]
for ii in range(1,len(XX)):
  Unrm[:,ii] = abs(II[ii]-II[ii-1])*V2d[:,ii] + abs(JJ[ii]-JJ[ii-1])*U2d[:,ii] 

# For plotting - fill land/bottom and smooth 2D field:
A2di = mmom6.fill_bottom(Unrm, ZZ, Hbtm)

clrmp = mclrmp.colormap_ssh(cpos='Reds')
clrmp.set_bad(color=[1,1,1])
rmax = 0.1
rmin = -rmax

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

lon1 = XX[0]
lat1 = YY[0]
lon2 = XX[-1]
lat2 = YY[-1]
stxt = f"X-axis: Distance (km) along section from {lon1:5.2f}W, {lat1:5.2f}N to {lon2:5.2f}W, {lat2:5.2f}N"

sttl = f"{expt} {varnm} avrg: {yrs}/{mms}/{dds} - {YR}/{MM}/{DD} {sctnm}"
btx = 'plot_xsectUV.py'
mxsct.plot_xsection(A2di, Xdist, ZM, Hbtm, Xdist, clrmp, rmin, rmax, \
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
  ax1.axis('scaled')

