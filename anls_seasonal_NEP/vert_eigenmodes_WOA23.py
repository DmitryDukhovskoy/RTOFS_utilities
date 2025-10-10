"""
  Compute and plot vertical eigenmodes computed from 
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
import netCDF4
from netCDF4 import Dataset as ncFile
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
YR      = 1993
MM      = 8  # for seasonal indicate month, for annual: MM = 13
mode    = 2 # baroclinic mode to plot

grd=0.25
if grd==0.25:
  cgrd=4
woa='woa23'

seas, decade, yr1_dec, yr2_dec = manseas.season_decade_woa(YR,MM)

woa_seas = {"13": "Jan-Mar",
            "14": "Apr-Jun",
            "15": "Jul-Spt",
            "16": "Oct-Dec",
            "0": "annual"}

urlBase = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/'
urlT    = f"{urlBase}temperature/netcdf/{decade}/0.25/"
urlS    = f"{urlBase}salinity/netcdf/{decade}/0.25/"
tfnm    = f"woa23_{decade}_t{seas:02d}_{cgrd:02d}.nc"
sfnm    = f"woa23_{decade}_s{seas:02d}_{cgrd:02d}.nc"
# All period average:
#tfnm=f'{woa}_decav_t{seas:02d}_{cgrd:02d}.nc'
#sfnm=f'{woa}_decav_s{seas:02d}_{cgrd:02d}.nc'

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  dmm = np.copy(dmm0)
  return dmm


# Get MOM6 topo for filling bottom in WOA data:
#fyaml = 'pypaths_gfdlpub.yaml'
#with open(fyaml) as ff:
#  gridfls = safe_load(ff)
#pthtopo    = gridfls['MOM6_NEP']['seasonal_fcst']['pthgrid']
#fgrid      = gridfls['MOM6_NEP']['seasonal_fcst']['fgrid']
#ftopo_mom  = gridfls["MOM6_NEP"]["seasonal_fcst"]["ftopo"]
#hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
#dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
## Hgrid lon. lat:
#hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
#HH = dstopo_nep['depth'].data
#HH = np.where(HH < 1.e-20, np.nan, HH)
#HH = -HH
#HH = np.where(np.isnan(HH), 1., HH)

# Get lon/lat
iz   = 0
furl = os.path.join(urlT,tfnm)
ZZ  = read_field(furl,'depth')
ZZ  = -abs(ZZ)
latW = read_field(furl,'lat')
lonW0 = read_field(furl,'lon')
jdm  = len(latW)
idm  = len(lonW0)


# Read T/S:
furl = os.path.join(urlT,tfnm)
var_read = 't_an'
A3d = read_field(furl,var_read)
A3d = np.where(A3d > 1.e10, np.nan, A3d)
# reshaffle to have -180/180 lon inside the domain
A3d, lonW = mmisc.shuffle3D_lon180(A3d, lonW0)

# Find WOA indices to match NEP section
# Get bottom profile for plotting
IsWOA, JsWOA, dsetBtm = manseas.xsct_segments_woa(sctnm, lonW, latW)

LONW = np.zeros((jdm,idm))
LATW = np.zeros((jdm,idm))
for ii in range(idm):
  LATW[:,ii]=latW
for jj in range(jdm):
  LONW[jj,:]=lonW

IJ      = np.column_stack((IsWOA, JsWOA))
DX, DY  = mmom6.dx_dy(LONW,LATW)
SGMT    = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive', check_pole=False)
II      = SGMT.I_indx
JJ      = SGMT.J_indx
nLeg    = SGMT.Leg_number
#II, JJ = mmisc.xsect_indx(Is, Js)
hLsgm1  = SGMT.half_Lsgm1
hLsgm2  = SGMT.half_Lsgm2
XX      = lonW[II]
YY      = latW[JJ]
LSgm    = np.zeros((len(II)))  # total segment length = half1 + half2
for ik in range(len(II)):
   LSgm[ik] = hLsgm1[ik] + hLsgm2[ik]

# Get section:
T2d = A3d[:,JJ,II].squeeze()

furl = os.path.join(urlS,sfnm)
var_read = 's_an'
A3d = read_field(furl,var_read)
A3d = np.where(A3d > 1.e10, np.nan, A3d)
A3d, _ = mmisc.shuffle3D_lon180(A3d, lonW0)
S2d = A3d[:,JJ,II].squeeze()

# Land-sea mask
LMsk = A3d[0,:,:].squeeze()
LMsk = np.where(np.isnan(LMsk), 1, -10)

#ZM  = mmom6.zz2zm(ZZ)

# Get N2 profiles:
# N2 profiles are at new depth levels: half-way between T/S points
# 1st N2 value is below surface, 
# Use BC (phi[0] = 0) to solve e/problem
# Z_phi[0] is not 0
ZZ[0] = 0.
N2, Z_phi = manseas.calc_N2(T2d, S2d, ZZ, YY)

ZZphi = np.zeros((len(Z_phi)+1))  # add surface to match Phi[0] surface
ZZphi[1:] = Z_phi
# Find layer thicknesses between phi-pnts - where density is
# either layer mid-points or layer interfaces - depending where
# T/S points are
dZlr_phi = abs(np.diff(ZZ))

# Solve Sturm-Liouville:
print(f'Solving Strum-Liouville, requested mode={mode}')
npnts  = len(II)
EVct   = np.zeros((T2d.shape))*np.nan
CPhase = np.zeros((npnts))*np.nan
Rbrcl  = np.zeros((npnts))*np.nan
for ik in range(npnts):
  N2z = N2[:,ik].copy()
  if np.isnan(N2z[0]): continue

# Find bottom:
# Escape coastal regions < 20 m
  k1 = np.where(np.isnan(N2z))[0]
  if k1.size:
    kbtm = k1.min()-1
# Fill below bottom:
    N2z[kbtm+1:]=N2z[kbtm]
    zbtm = ZZ[kbtm+1]  # bottom depth = bottom interface of the near-bottom gr.cell
  else:
# Case when Bottom is deeper than WOA last depth level
# extend last grid cell to some depth > WOA deepest level
# Or use MOM6 topo - see example for R rossby from WOA data
# Make sure that zbtm does not = Z_phi[kbtm]
# to avoid singularities in matrix AA
#    print('Deep site: ii={0}, jj={1}, Hbtm={2:6.1f}'.\
#          format(ii,jj,HH[jj,ii]))
    kbtm = N2z.shape[0]-1
    zbtm = ZZ[-1]-100.

    if abs(Z_phi[kbtm] - zbtm) < 1.e-1:
      zbtm = Z_phi[kbtm] - 1.e-1

  # Filter N-point running mean
  N2zF = msolv.runmn_dz(N2z, dZlr_phi, mnwnd=5)

  # add NaNs for bottom values:
  # kbtm - near-bottom value for phi(z)
  Hb0 = zbtm
  if abs(Hb0) < abs(Z_phi[-1]): 
    N2zF[kbtm+1:] = np.nan
  else:
# Truncate profile to the deepest level, make a bottom
# to keep Phi correct size
    N2zF[kbtm] = np.nan

# Skip shallow regions
# for small kbtm - problem creating a matrix and solving it
# Skip shallow regions
# for small kbtm - problem creating a matrix and solving it
  if kbtm < 5:
    print(f'Shallow location: ik={ik} kbtm={kbtm} Hbtm={Hb0:.2f} m')
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
  clrmp.set_bad(color=[0.1, 0.1, 0.1])

# Indices:
#Xdist = np.arange(len(Hbtm))
# 
# Distance along the section
# normalize by the total distance
Lsection = mmisc.dist_sphcrd(YY[-1],XX[-1],YY[0],XX[0]) # total length of the section, m
Xdist = np.cumsum(LSgm)*1.e-3
Hbtm  = dsetBtm['Hbtm_section'].data
Xdist_mom = dsetBtm['Dist_section'].data
Xdist = Xdist/Xdist[-1]*Xdist_mom[-1]    # normalize distance to match MOM6

xl1  = min(Xdist)
xl2  = max(Xdist)

lon1 = XX[0]
lat1 = YY[0]
lon2 = XX[-1]
lat2 = YY[-1]
stxt = f"X-axis: Distance (km) along section from {lon1:5.2f}W, {lat1:5.2f}N to {lon2:5.2f}W, {lat2:5.2f}N"

seas_nm = woa_seas[f"{seas}"]
sttl = f"WOA23 eigenvector {mode} decade:{yr1_dec}-{yr2_dec} {seas_nm} {sctnm}"

btx = 'vert_eigenmodes_WOA23.py'
if mode == 1:
  mxsct.plot_xsection(EVct, Xdist, ZZphi, Hbtm, Xdist, clrmp, rmin, rmax, \
                    xl1, xl2, sttl=sttl, stxt=stxt,  \
                    plot_map=LMsk, I_indx=II, J_indx=JJ, btx=btx, patch_btm=False)
else:
  mxsct.plot_xsection(EVct, Xdist, ZZphi, Hbtm, Xdist, clrmp, rmin, rmax, \
                    xl1, xl2, sttl=sttl, stxt=stxt, \
                    plot_map=LMsk, I_indx=II, J_indx=JJ, btx=btx, \
                    patch_btm=False, contours=[0.])

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


