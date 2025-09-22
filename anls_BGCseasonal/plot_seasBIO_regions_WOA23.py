"""
  Plot BGC fields for North East Pac. regions
  for different seasons
  from WOA23

  BGC fields are available on 1 and 5-dgr grids
  by seasons / months
  for 1965-2022 
  and 1971-2000 

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
import argparse

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
importlib.reload(manseas)

parser = argparse.ArgumentParser()
parser.add_argument("--varnm", help="o2, po4, no3, sio4", type=str, required=True)
parser.add_argument("--regnm", help="CalCur BeringChuk GulfAlaska", type=str, required=True)
parser.add_argument("--mseas", help="Month/season WOA: 1,..,12, 13-Winter, ..., 0-annual", \
                    type=int, required=True)
parser.add_argument("--zz", help="Aprx depth to plot, m ", type=float, required=True)
args = parser.parse_args()

varnm = args.varnm if args.varnm else None
regn_name = args.regnm if args.regnm else None
mseas = args.mseas if args.mseas else None
zz_plt = args.zz if args.zz is not None else None # ocean depth to plot

zz_plt = -abs(zz_plt)

cgrd = 1
woa = 'woa23'

woa_seas = {"13": "Jan-Mar",
            "14": "Apr-Jun",
            "15": "Jul-Spt",
            "16": "Oct-Dec",
            "00": "annual"}

urlBase = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/'
#'oxygen/netcdf/all/1.00/woa23_all_o00_01.nc'
url_o2    = os.path.join(f"{urlBase}","oxygen/netcdf/all/1.00/")
flnm_o2   = f"woa23_all_o{mseas:02d}_01.nc"
#silicate/netcdf/all/1.00/woa23_all_i01_01.nc
url_sio4  = os.path.join(f"{urlBase}","silicate/netcdf/all/1.00/")
flnm_sio4 = f"woa23_all_i{mseas:02d}_01.nc"
# nitrate/netcdf/all/1.00/woa23_all_n06_01.nc
url_no3   = os.path.join(f"{urlBase}","nitrate/netcdf/all/1.00/")
flnm_no3  = f"woa23_all_n{mseas:02d}_01.nc"
# phosphate/netcdf/all/1.00/woa23_all_p02_01.nc
url_po4   = os.path.join(f"{urlBase}","phosphate/netcdf/all/1.00/")
flnm_po4  = f"woa23_all_p{mseas:02d}_01.nc"


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
furl = os.path.join(url_o2, flnm_o2)
ZM  = read_field(furl,'depth')
ZM  = -abs(ZM)
ZM = -abs(ZM)
dZ = np.abs(ZM-zz_plt)
iz0 = np.argmin(dZ)
lr0  = iz0+1
zz0 = ZM[iz0]  # actual depth to be plotted
latW0 = read_field(furl,'lat')
lonW0 = read_field(furl,'lon')

# Read fields from WOA23
# Read objectively analyzed means
# also available - mean of unflagged fields within the grid cell
match varnm:
  case 'o2':
    furl = os.path.join(url_o2, flnm_o2)
    varnc = 'o_an'
  case 'po4':
    furl = os.path.join(url_po4, flnm_po4)
    varnc = 'p_an'
  case 'no3':
    furl = os.path.join(url_no3, flnm_no3)
    varnc = 'n_an'
  case 'sio4':
    furl = os.path.join(url_sio4, flnm_sio4)
    varnc = 'i_an'

A2d, lonW, latW, ix1, ix2, jx1, jx2 = read_2Dfield_subsample(furl, varnc, \
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

if 0 < mseas < 13:
  seas_nm = f"{mseas:02d}"
else:
  seas_nm = woa_seas[f"{seas}"]


rmin, rmax, tscntrs, tslabels = manseas.colormap_params(regn_name, varnm, zz0=zz0)

Ncmp = 200
log_scale = False
log_str = ''
match varnm:
  case 'o2':
    cff = 1.
    unts = 'mcromol/kg'
    clrmp = mclrmps.colormap_haline2(end_clr=[0.7,0.1,0], start_clr=[0.8,0.8,1])

  case 'po4':
    cff = 1.
    unts = 'mcromol/kg'
    clrmp = mclrmps.colormap_conc()

  case 'sio4':
    cff = 1.
    unts = 'mcromol/kg'
    clrmp = mclrmps.colormap_ice_thkn()

  case 'no3':
    cff = 1.
    log_scale = True
    log_str = 'log'
    unts = 'mcromol/kg'
    CLRS = [[1, 1, 1],
        [0.6, 0.02, 0.6],
        [0.2, 0.38, 1],
        [0., 0.8, 0.8],
        [0.4, 0.8, 0],
        [1, 1, 0.5],
        [1, 0.8, 0.6],
        [1, 0.6, 0],
        [0.7, 0.1, 0.1]]

    #clrmp = mclrmps.colormap_temp(clr_ramp=[1,1,1])
    clrmp = mclrmps.colormap_posneg_uneven(CLRS)


clrmp.set_bad(color=[0., 0., 0.])

if log_scale:
  JJ0,II0 = np.where(A2d <= 1.e-32)
  if len(JJ0) > 0:
    A2d[JJ0,II0] = np.nan
  lA2d = np.log(A2d)
  if len(JJ0) > 0:
    lA2d[JJ0,II0] = 0.

  A2d = lA2d.copy()



sttl = f"WOA23 {log_str} {varnm} obj.mean 1965-2022 {seas_nm} z={zz0:8.1f} m"


# Stereographic projection:
from mpl_toolkits.basemap import Basemap, cm

lon0, lat0, height, width = manseas.stereogr_params_regions(regn_name)

m = Basemap(width=width, height=height, resolution='l',\
            projection='stere', lat_ts=55, lat_0=lat0, lon_0=lon0)

xR, yR = m(LONW, LATW)

btx = 'plot_seasBIO_regions_WOA23.py'

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





