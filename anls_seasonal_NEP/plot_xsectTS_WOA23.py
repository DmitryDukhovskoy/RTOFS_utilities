"""
  Plot T/S vertical sections from WOA23
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
from netCDF4 import Dataset as ncFile
import importlib
import yaml
from yaml import safe_load

PPTHN = []
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
import mod_read_hycom as mhycom
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_anls_seas as manseas
import mod_plot_xsections as mxsct

varnm   = 'salin'  # temp / salin
sctnm   = 'xsct_BerSea'  
YR      = 1993
MM      = 11  # for seasonal indicate month, for annual: MM = 13

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

def lookup_ncvar(nc):
  ii=0
  for var in nc.variables.values():
    ii+=1
    print('--------\n')
    print('Var # {0}'.format(ii))
    print(var)
 
iz   = 0
if varnm == 'temp':
  furl = os.path.join(urlT,tfnm)
  var_read = 't_an' 
else:
  furl = os.path.join(urlS,sfnm)
  var_read = 's_an'

A3d = read_field(furl,var_read)
A3d = np.where(A3d > 1.e10, np.nan, A3d) 
kdm, jdm, idm = A3d.shape

ZZ  = read_field(furl,'depth')
ZZ  = -abs(ZZ)
latW = read_field(furl,'lat')
lonW = read_field(furl,'lon')

# reshaffle to have -180/180 lon inside the domain
A3d, lonW = mmisc.shuffle3D_lon180(A3d, lonW) 

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


# Land-sea mask
LMsk = A3d[iz,:,:].squeeze()
LMsk = np.where(np.isfinite(LMsk), -10, 1)

A2d = A3d[:,JJ,II].squeeze()
if varnm == 'salin':
#  clrmp = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  clrmp = mclrmps.colormap_salin2()
  rmin = pthseas['ANLS_NEP'][sctnm]['smin']
  rmax = pthseas['ANLS_NEP'][sctnm]['smax']
elif varnm == 'temp':
  clrmp = mclrmps.colormap_temp(clr_ramp=[0.9,0.8,1])
  rmin = pthseas['ANLS_NEP'][sctnm]['tmin']
  rmax = pthseas['ANLS_NEP'][sctnm]['tmax']
  if sctnm == 'xsct_BerSea' and (MM>6 and MM<=12):
    rmax = 10.
  
clrmp.set_bad(color=[0.5, 0.5, 0.5])
XX    = lonW[II]
YY    = latW[JJ]
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
sttl = f"WOA23 {varnm} decade:{yr1_dec}-{yr2_dec} {seas_nm} {sctnm}"
btx = 'plot_xsectTS_WOA23.py'
mxsct.plot_xsection(A2d, Xdist, ZZ, Hbtm, Xdist_mom, clrmp, rmin, rmax, \
                    xl1, xl2, sttl=sttl, stxt=stxt, \
                    plot_map=LMsk, I_indx=II, J_indx=JJ, btx=btx, patch_btm=False)





