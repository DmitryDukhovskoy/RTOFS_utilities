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
import argparse
import xarray

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

parser = argparse.ArgumentParser()
parser.add_argument("--varnm", help="salin or temp", type=str, required=True)
parser.add_argument("--yr", help="Year within a decade to plot", type=int, required=True)
parser.add_argument("--mm", help="Month within WOA season, 1,..., 12, 13-annual", type=int, required=True)
#parser.add_argument("--lr", help="WOA vert layer 1,...,102", type=int, required=True)
parser.add_argument("--zz", help="Depth to plot, m >0", type=float, required=True)
args = parser.parse_args()

# Note WOA tempis in situ !!!
varnm = args.varnm if args.varnm else None
YR = args.yr if args.yr else None
MM = args.mm if args.mm else None
zz_plt = args.zz if args.zz else None
if zz_plt is not None:
  zz_plt = -abs(zz_plt)

#lr0 = args.lr if args.lr else None  # ocean layers from 1, ..., 102
                                    # lr 11 =-50, lr 21 =-102 m, lr 25 = -200 m

regn_name = 'ArcticOcean'

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

# Lon/lat of the domain:
lon1 = -180.
lon2 = 180.
lat1 = 55.
lat2 = 90.

# Arctic domain:
pthtopo_arc = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'
dfgrid = os.path.join(pthtopo_arc,'ocean_mask_ARC12.nc')
dsgrid_arc = xarray.open_dataset(dfgrid)
LONM = dsgrid_arc['x'].data
LATM = dsgrid_arc['y'].data


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

  # Subsample domain of interest:
  ix1 = np.argmin(np.abs(lonW0-lon1))
  ix2 = np.argmin(np.abs(lonW0-lon2))+1
  jx1 = np.argmin(np.abs(latW0-lat1))
  jx2 = np.argmin(np.abs(latW0-lat2))+1
  A2d = A3d[iz0,jx1:jx2,ix1:ix2].squeeze()

  lonW = lonW0[ix1:ix2]
  latW = latW0[jx1:jx2]

  return A2d, lonW, latW, ix1, ix2, jx1, jx2

# Get lon/lat, depths:
#iz0 = lr0-1
furl = os.path.join(urlT,tfnm)
ZM  = read_field(furl,'depth')
ZM  = -abs(ZM)
dZ = abs(ZM-zz_plt)
iz0 = np.argmin(dZ)
lr0  = iz0+1
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


rmin, rmax, tscntrs, tslabels = manseas.colormap_params(regn_name, varnm, zz0=zz0)

seas_nm = woa_seas[f"{seas}"]
sttl = f"WOA23 {varnm} decade:{yr1_dec}-{yr2_dec} {seas_nm} z={zz0:8.1f} m"

lnd_clr = [0.4,0.4,0.4]

if varnm == 'salin' or varnm == 'salt':
  #clrmp = mclrmps.colormap_haline2()
  clrmp = mclrmps.colormap_haline()
  clrmp.set_bad(color=lnd_clr)
  rmin = 24.
  rmax = 36.
#  clrmp.set_under(color=[0.6, 0.6, 0.6])
elif varnm == 'temp' or varnm == 'potT':
  clrmp = mclrmps.colormap_temp(clr_ramp=[0.6,0.4,.9])
  #clrmp = mclrmps.colormap_temp2()
  clrmp.set_bad(color=lnd_clr)
  sttl = f"WOA23 Tpotential decade:{yr1_dec}-{yr2_dec} {seas_nm} z={zz0:8.1f} m"
  rmin = -2.
  rmax = 18.
elif varnm == 'ssh':
  clrmp = mclrmps.colormap_ssh(nclrs=200)
  rmin = -0.5
  rmax = 0.5

# Stereographic projection:
from mpl_toolkits.basemap import Basemap, cm
lon0 = 0.
lat0 = 55.
res  = 'l'
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
xR, yR = m(LONW, LATW)

btx = 'plot_seasTS_Arctic_WOA23.py'

print('Plotting ...')
plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()

ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

m.fillcontinents(color=lnd_clr, lake_color='white')
#m.drawcoastlines(color='w')
m.drawparallels(np.arange(-40,88.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))
cntr_clr = [0.5, 0.5, 0.5]

img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)

#if len(tscntrs) > 0:
#  CS = m.contour(xR, yR, A2d, tscntrs, colors=[cntr_clr], linestyles='solid', linewidths=1)
#  if len(tslabels) > 0:
#    ax1.clabel(CS, tslabels, inline=1, fontsize=10)

ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


xl1 = 2.7e6
xl2 = 10.2e6
yl1 = 3.6e6
yl2 = 11.e6
ax1.axis('scaled')
#ax1.set_xlim([xl1,xl2])
#ax1.set_ylim([yl1,yl2])

ax3 = plt.axes([0.2,0.1,0.1,0.05])
ax3.text(0.5, 0.5, 'WOA23',
         fontsize=12,
         ha='center', va='center')
ax3.set_xticks([])
ax3.set_yticks([])

bottom_text(btx, fsz=8, pos=[0.05, 0.03])


