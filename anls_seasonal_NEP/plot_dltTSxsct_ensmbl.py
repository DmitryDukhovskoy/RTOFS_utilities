"""
  Plot sections with vertical distribution of T or S
  differences wrt to a reference ensmble run
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
varnm  = 'temp'  # temp (potential) / salin
sctnm  = 'xsct_EOB' 

# Start of the run 
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1    
nensR  = 1  #ens # for reference ensemble run
month_av = True  # monthly average data
dnmbR  = mtime.datenum([1994,3,15])  # day/month to plot

expt    = "seasonal_fcst"
runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nensR:02d}'
dnmbS   = mtime.datenum([YRS,MOS,DDS]) 
dv_start = mtime.datevec(dnmbS)

dvR = mtime.datevec(dnmbR)
print(f'Expt: {expt} Run: {runname} Plot date: {dvR[0]}/{dvR[1]}/{dvR[2]}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

if expt == 'seasonal_fcst':
  pthfcst  = pthseas['MOM6_NEP'][expt]['pthwoutp'].format(runname=runname)
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
ocnfld = 'oceanm'
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


dsetR   = xarray.open_dataset(dfmom6)

ZM  = -dsetR['zl'].data
ZZ  = mmom6.zm2zz(ZM)

# Monthly mean flds:
yrR, moR = dvR[:2]
if month_av:
  F2dR = manseas.monthly_avrg_vertxsect(pthfcst, yrR, moR, JJ, II, varnm)
else:
  if varnm == 'temp' or varm == 'potT':
    F2dR = dsetR['potT'].data[0,:,JJ,II].squeeze()
  elif varnm == 'salin' or varnm == 'salt':
    F2dR = dsetR['salt'].data[0,:,JJ,II].squeeze()
  F2dR = np.transpose(F2dR)

# Distance along the section
# normalize by the total distance
Lsection = mmisc.dist_sphcrd(YY[-1],XX[-1],YY[0],XX[0]) # total length of the section, m
Xdist = np.cumsum(LSgm)
Xdist = Xdist-Xdist[0]
Xdist = Xdist/Xdist[-1]*Lsection*1.e-3  # normalized, km

xl1  = min(Xdist)
xl2  = max(Xdist)


# For plotting - fill land/bottom and smooth 2D field:
#A2di = mmom6.fill_bottom(A2d, ZZ, Hbtm)
ENSR = [2,3,4,5,6,7]
NensR = len(ENSR)
for iens in range(NensR):
  nens = ENSR[iens]
  runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
  pthfcst = pthseas['MOM6_NEP'][expt]['pthwoutp'].format(runname=runname)
  pthfcst = os.path.join(pthfcst,f'{ocnfld}_{dvR[0]}{dvR[1]:02d}')

  print(f"Processing ens={nens:02d} {pthfcst}")
  if month_av:
    F2d = manseas.monthly_avrg_vertxsect(pthfcst, yrR, moR, JJ, II, varnm)
  else:
    YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthfcst, dnmbR, fld=ocnfld)
    flocn_name = pthseas['MOM6_NEP'][expt]['focname'].format(YR=YR0, jday=jday0)
    dfmom6 = os.path.join(pthfcst, flocn_name)
    dset   = xarray.open_dataset(dfmom6)
    if varnm == 'temp' or varm == 'potT':
      F2d = dset['potT'].data[0,:,JJ,II].squeeze()
    elif varnm == 'salin' or varnm == 'salt':
      F2d = dset['salt'].data[0,:,JJ,II].squeeze()
    F2d = np.transpose(F2d)

#  dF = F2d - F2dR
#  print(f"diff min/max: {np.nanmin(dF)}/{np.nanmax(dF)}")
  if iens == 0:
    dim1 = "depth"
    dim2 = "distance"
    darr_var = xarray.DataArray(F2d, dims=(dim1,dim2), 
                coords={dim1: ZM, dim2: Xdist})
    dset2D = xarray.Dataset({f"{varnm}_e{nens:02d}": darr_var}).copy()
  else:
    darr_var = xarray.DataArray(F2d, dims=(dim1,dim2), 
                coords={dim1: ZM, dim2: Xdist})
    dset_var = xarray.Dataset({f"{varnm}_e{nens:02d}": darr_var})
    dset2D = xarray.merge([dset2D, dset_var])

# ===================
# Plotting
# ===================

dltx = 0.06
dlty = 0.07
dx = 0.27
dy = 0.27
xl = 0.06
yb = 0.15
FPOS = [[xl, yb+dy+dlty, dx, dy],
        [xl+dx+dltx, yb+dy+dlty, dx, dy],
        [xl+2*(dx+dltx), yb+dy+dlty, dx, dy],
        [xl, yb, dx, dy],
        [xl+dx+dltx, yb, dx, dy],
        [xl+2*(dx+dltx), yb, dx, dy]]

clrmp  = mutil.colormap_ssh(nclrs=200)
rmin = -1.
rmax = 1.

Xsgm = dset2D['distance'].data
Hbtm = np.where(Hbtm > 0, 0., Hbtm)
segm_nm = sctnm

match sctnm:
  case 'xsct_SOB':
    xl1 = 0
  case 'xsct_NOB':
    xl1 = 200.
  case 'xsct_EOB':
    xl1 = 0.
  case 'xsct_WOB':
    xl1 = 0.
  case _:
    raise Exception(f"{sctnm} not recognized")
xl2  = max(Xsgm)

dnmb0 = dnmbR
dv0   = mtime.datevec(dnmb0)

btx = 'plot_dltTSxsct_ensmbl.py' 


from matplotlib.patches import Polygon

plt.ion()

fgnmb = 1
fig1 = plt.figure(fgnmb,figsize=(12,9))
plt.clf()

for ie in range(NensR):
  nens = ENSR[ie]
  ax1 = plt.axes(FPOS[ie])
  vards = f"{varnm}_e{nens:02d}"
  F2d = dset2D[vards].data.squeeze()
  dF  = F2d - F2dR


  if month_av: 
    dstr = f'Month avrg {dvR[0]}/{dvR[1]}'
  else:
    dstr = f'5day avrg {dvR[0]}/{dvR[1]}/{dvR[2]}'
    
  dstart = f'{dv_start[0]}/{dv_start[1]}/{dv_start[2]}'
  sttl = f'ens={nens:02d}'

  if ie == 0:
    btx0     = btx
    stxt     = f'NEP seas fcst climOB init: {dstart} diff {varnm} {sctnm} wrt ens={nensR}, {dstr}'
    clrb_pos = [0.1, 0.07, 0.8, 0.04]
  else:
    btx0     = ''
    stxt     = ''
    clrb_pos = []

  mutob.plot_Nxsections(dF, Xsgm, ZM, Hbtm, Xsgm, clrmp, \
            rmin, rmax, xl1, xl2, fig1, ax1, sttl=sttl, stxt=stxt,\
            btx=btx0, clrb_ornt='horiz', \
            clrb_pos=[0.1, 0.07, 0.8, 0.02], txt_pos=[0.1, 0.9, 0.8, 0.09])


