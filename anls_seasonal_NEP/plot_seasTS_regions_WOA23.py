"""
  Plot T/S fields for North East Pac. region
  for different seasons
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
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_anls_seas as manseas
import mod_plot_xsections as mxsct
import matplotlib as mtplt
import mod_regmom as mregmom
import mod_colormaps as mclrmps

varnm     = 'salin'  # temp / salin  Note: WOA temp is in situ !!!
#regn_name = 'CalCur'  
regn_name = 'BeringChuk'
YR        = 2011
MM        = 8  # to find season, indicate month, for annual: MM = 13

lr0  = 11  # ocean layers from 1, ..., 102
          # lr 11 =-50, lr 21 =-102 m, lr 25 = -200 m

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

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

# Find lon/lat of the domain:
II = pthseas['ANLS_NEP'][regn_name]['II']
JJ = pthseas['ANLS_NEP'][regn_name]['JJ']
xlim1 = min(II)
xlim2 = max(II)
ylim1 = min(JJ)
ylim2 = max(JJ)

pthtopo    = pthseas['MOM6_NEP']['seasonal_daily']['pthgrid']
fgrid      = pthseas['MOM6_NEP']['seasonal_daily']['fgrid']
dfgrid_mom = os.path.join(pthtopo, fgrid)
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

lon1 = hlon[ylim2, xlim1]
lon2 = hlon[ylim1, xlim2]
lat1 = hlat[ylim1, xlim1]
lat2 = hlat[ylim2, xlim2]


def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  dmm = np.copy(dmm0)
  return dmm

def read_2Dfield_subsample(furl, var_read, lonW0, latW0, lon1, lon2, lat1, lat2, iz0):
  A3d = read_field(furl,var_read)
  A3d = np.where(A3d > 1.e10, np.nan, A3d)
  # reshaffle to have -180/180 lon inside the domain
  A3d, lonWr = mmisc.shuffle3D_lon180(A3d, lonW0)

  # Subsample domain of interest:
  ix1 = np.argmin(np.abs(lonWr-lon1))
  ix2 = np.argmin(np.abs(lonWr-lon2))+1
  jx1 = np.argmin(np.abs(latW0-lat1))
  jx2 = np.argmin(np.abs(latW0-lat2))+1
  A2d = A3d[iz0,jx1:jx2,ix1:ix2].squeeze()

  lonW = lonWr[ix1:ix2]
  latW = latW0[jx1:jx2]

  return A2d, lonW, latW, ix1, ix2, jx1, jx2

# Get lon/lat
iz0 = lr0-1
furl = os.path.join(urlT,tfnm)
ZM  = read_field(furl,'depth')
ZM  = -abs(ZM)
zz0 = ZM[iz0]
latW0 = read_field(furl,'lat')
lonW0 = read_field(furl,'lon')

# Read T/S:
# For in situ T, need S to convert it to potential T
furl = os.path.join(urlT,tfnm)
if varnm == 'temp':
  var_read = 't_an'
  furl = os.path.join(urlT,tfnm)
  furlS = os.path.join(urlS,sfnm)
elif varnm == 'salin':
  var_read = 's_an'
  furl = os.path.join(urlS,sfnm)

S2d = []
A2d, lonW, latW, ix1, ix2, jx1, jx2 = read_2Dfield_subsample(furl, var_read, \
                                      lonW0, latW0, lon1, lon2, lat1, lat2, iz0)
if varnm == 'temp' and abs(zz0) > 10.:
  S2d, _, _, _, _, _, _ = read_2Dfield_subsample(furlS, 's_an', \
                                      lonW0, latW0, lon1, lon2, lat1, lat2, iz0)

jdm  = len(latW)
idm  = len(lonW)
LONW = np.zeros((jdm,idm))
LATW = np.zeros((jdm,idm))
for ii in range(idm):
  LATW[:,ii]=latW
for jj in range(jdm):
  LONW[jj,:]=lonW

# Convert in situ --> potential T at z=0:
if varnm == 'temp' and abs(zz0) > 10.:
  import mod_regmom as mregmom
  print('Converting T in situ --> T potential')
  Tpot = mregmom.insitu2pot_2D(A2d, S2d, zz0, LATW)
  A2d = Tpot.copy()

# Define NEP domain, find in WOA:
lon_mom_s = hlon[ylim1, xlim1:xlim2+1]
lat_mom_s = hlat[ylim1, xlim1:xlim2+1]
idx_woa_s, jdx_woa_s = mmisc.find_closest_indx(lon_mom_s, lat_mom_s, lonW, latW)

lon_mom_n = hlon[ylim2, xlim1:xlim2+1]
lat_mom_n = hlat[ylim2, xlim1:xlim2+1]
idx_woa_n, jdx_woa_n = mmisc.find_closest_indx(lon_mom_n, lat_mom_n, lonW, latW)

lon_mom_w = hlon[ylim1:ylim2+1, xlim1]
lat_mom_w = hlat[ylim1:ylim2+1, xlim1]
idx_woa_w, jdx_woa_w = mmisc.find_closest_indx(lon_mom_w, lat_mom_w, lonW, latW)

lon_mom_e = hlon[ylim1:ylim2+1, xlim2]
lat_mom_e = hlat[ylim1:ylim2+1, xlim2]
idx_woa_e, jdx_woa_e = mmisc.find_closest_indx(lon_mom_e, lat_mom_e, lonW, latW)

IIG, JJG = mmisc.connect_segments([idx_woa_w, idx_woa_n, idx_woa_e, idx_woa_s], \
                                  [jdx_woa_w, jdx_woa_n, jdx_woa_e, jdx_woa_s])

#X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
#MS, _, _ = mmisc.inpolygon_v2(X, Y, IIG, JJG)  # 
#Jout, Iout = np.where(MS<1)

xNEP = LONW[JJG,IIG]
yNEP = LATW[JJG,IIG]

rmin, rmax, tscntrs, tslabels = manseas.colormap_params(regn_name, varnm, zz0=zz0)

seas_nm = woa_seas[f"{seas}"]
sttl = f"WOA23 {varnm} decade:{yr1_dec}-{yr2_dec} {seas_nm} z={zz0:8.1f} m"

if varnm == 'salin' or varnm == 'salt':
  clrmp = mclrmps.colormap_haline2()
  clrmp.set_bad(color=[0., 0., 0.])
#  clrmp.set_under(color=[0.6, 0.6, 0.6])
elif varnm == 'temp' or varnm == 'potT':
  clrmp = mclrmps.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[0.,0.,0.])
  sttl = f"WOA23 Tpotential decade:{yr1_dec}-{yr2_dec} {seas_nm} z={zz0:8.1f} m"
elif varnm == 'ssh':
  clrmp = mclrmps.colormap_ssh(nclrs=200)
  rmin = -0.5
  rmax = 0.5

# Stereographic projection:
from mpl_toolkits.basemap import Basemap, cm
match regn_name:
  case 'CalCur':
    width  = 4000*1.e3
    height = 4000*1.e3
    lat0   = 33.5
    lon0   = -128.
  case 'BeringChuk':
    width  = 3300*1.e3
    height = 3700*1.e3
    lat0   = 65.
    lon0   = -175.
 
m = Basemap(width=width, height=height, resolution='l',\
            projection='stere', lat_ts=55, lat_0=lat0, lon_0=lon0)

xR, yR = m(LONW, LATW)

btx = 'plot_seasTS_regions_WOA23.py'

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
 
ax1 = manseas.plot_stereogr_axis(fig1, m, xR, yR, A2d, clrmp, rmin, rmax, \
                       btx=btx, tscntrs=tscntrs, tslabels=tslabels, sttl=sttl)

# Plot NEP domain:
plt.sca(ax1)
xdom, ydom = m(xNEP, yNEP)
m.plot(xdom, ydom, 'w-')

#match regn_name:
#jdim, idim = A2d.shape
#ilim1 = 0
#ilim2 = idim
#jlim1 = 0
#jlim2 = jdim
#
#  case 'CalCur':
#    manseas.plot2D_CalCur(A2d, clrmp, rmin, rmax, ilim1, ilim2, jlim1, jlim2, \
#                  fgnmb=1, btx=btx, tscntrs=tscntrs, tslabels=tslabels, sttl=sttl, \
#                  hlon=LONW, hlat=LATW, HH=[])





