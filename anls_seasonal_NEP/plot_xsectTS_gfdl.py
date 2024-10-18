"""
  Plot sections with vertical distribution of T or S
  GFDL simulations (Liz's)
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
import mod_mom6 as mmom6
import mod_utils_ob as mutob
importlib.reload(mutob)

# experiment: year start, month start, ...
# monthly mean fields saved
# change month to plot desired date output - # of days since start date
varnm  = 'temp'  # temp/ salin
sctnm  = 'xsct_BerSea' 
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1    
nens   = 2
yrrun  = 1993  # year to plot
morun  = 11     # month to plot
mdrun  = 15 # month day

#expt  = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
dnmbS = mtime.datenum([YRS,MOS,DDS]) 
dnmbR = mtime.datenum([yrrun, morun, mdrun])
# What month number:
nmonth = (yrrun-YRS)*12 + morun  # consecutive month #

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

run_nm = 'MOM6_NEP_GFDL'
expt   = 'NEP_BGC_seas'
pthfcst    = pthseas[run_nm][expt]['pthoutp']
pthtopo    = pthseas[run_nm][expt]['pthgrid']
fgrid      = pthseas[run_nm][expt]['fgrid']
ftopo_mom  = pthseas[run_nm][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

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


if varnm == 'temp':
  dfmom6 = os.path.join(pthfcst,'ocean_monthly_z.199301-201912.thetao.nc')
  dset   = xarray.open_dataset(dfmom6)
  A2d = dset['thetao'].data[nmonth-1,:,JJ,II].squeeze()
else:
  A2d = dset['so'].data[nmonth-1,:,JJ,II].squeeze()
A2d = np.transpose(A2d)

ZM  = -dset['z_l'].data
ZZ  = mmom6.zm2zz(ZM)

day_time = dset['time'].data[nmonth-1]
day_str  = np.datetime_as_string(day_time)
YR0    = int(day_str.split("-")[0])
MM0    = int(day_str.split("-")[1])
DD0    = int(day_str.split("-")[2].split("T")[0])

print(f"Plotting monthly mean {varnm} for {YR0}/{MM0}/{DD0}")

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

lon1 = XX[0]
lat1 = YY[0]
lon2 = XX[-1]
lat2 = YY[-1]
stxt = f"X-axis: Distance (km) along section from {lon1:5.2f}W, {lat1:5.2f}N to {lon2:5.2f}W, {lat2:5.2f}N"

sttl = f"{expt} {varnm} monthly avrg {YR0}/{MM0} {sctnm}"
btx = 'plot_xsectTS_gfdl.py'
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
  ax1.contour(hlat,[x for x in range(10,88,5)], linestyles='solid', colors=[(0.8, 0.8, 0.8)])
  ax1.contour(hlat, [55.125], linestyles='solid', colors=[(1., 0.4, 0.)])
  ax1.axis('scaled')

# Bind the button_press_event with the onclick() method
  fig1.canvas.mpl_connect('button_press_event', onclick)


