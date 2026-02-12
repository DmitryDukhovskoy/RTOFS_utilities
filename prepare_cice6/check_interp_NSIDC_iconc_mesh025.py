"""
  Check Interpolated NSIDC sea ice concentration 
  to MOM6/CICE6 mesh025 grid

  gmapi indices: get_gmapi_NSIDC_to_mesh025.py

  NSIDC fields from 
  https://noaadata.apps.nsidc.org/NOAA/G02202_V6/north/daily/2025/

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
from mpl_toolkits.basemap import Basemap, cm
import argparse
                   
# Append custom module paths
PPTHN = None
if 'PPTHN' not in locals() or PPTHN is None:
  cwd = os.getcwd()    
  parts = cwd.split(os.sep)
  if 'python' in parts:
    idx = parts.index('python')
    PPTHN = os.sep + os.path.join(*parts[:idx + 1])
  else:
    raise RuntimeError("Directory 'python' not found in current working directory path.")

sys.path.extend([
    os.path.join(PPTHN, 'MyPython', 'hycom_utils'),
    os.path.join(PPTHN, 'MyPython', 'draw_map'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])


from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

YR = 2025
MM = 1 
DD = 15

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", type=str, required=True)
parser.add_argument("--yr", help=f"year of NSIDC data, default={YR}", type=int)
parser.add_argument("--mm", help=f"month of NSIDCS data to plot, default={MM}", type=int)
parser.add_argument("--dd", help="month day to plot, default={DD}", type=int)
args = parser.parse_args()
  
regn = args.regn if args.regn else None
YR   = args.yr if args.yr else YR
MM   = args.mm if args.mm else MM
DD   = args.dd if args.dd else DD
  
syst_info = os.uname() 
machine = syst_info.nodename
  
if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
else:
  print("Unknown machine:", machine)
    
fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)
    
# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
    
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape

def read_NSIDC(YR,MM,DD,regn,pthnsidc,varnm):
  if regn == 'south':
    flnsidc = f"sic_pss25_{YR}{MM:02d}{DD:02d}_am2_v06r00.nc"  
  else:
    flnsidc = f"sic_psn25_{YR}{MM:02d}{DD:02d}_am2_v06r00.nc"  

  with xarray.open_dataset(os.path.join(pthnsidc,flnsidc)) as ds_nsidc:
    A = ds_nsidc[varnm].data.squeeze()

  return A

# Get lon/lat for NSIDC data
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")
Xnsidc = read_NSIDC(YR,MM,1,regn,pthnsidc,'x')
Ynsidc = read_NSIDC(YR,MM,1,regn,pthnsidc,'y')

XX, YY = np.meshgrid(Xnsidc, Ynsidc, indexing='xy')
# Determine ellipsoid parameters from NSIDC information
# Note that Radius of ellipsoid WGS84 is typically referred to major semi-axis (equatorial radius)
if regn == 'south':
  slat0 = -70.  # latitude of 0 distortion, standard lat. 
  lon0  = -90.   # orientation of the 0 longitude wth to X axis on polar grid, not NSIDC is fliped upside down
  flat_inv = 298.279411123064
  ax_maj = 6378273.
  flat = 1./flat_inv  # flattening
  eccentr = np.sqrt(2*flat - flat**2)
  R_polar = ax_maj*(1.-flat)   # polar radius or semi-minor axis


  LON, LAT = mmisc.convert_polarXY_lonlat(XX,YY, North=False, E=eccentr, RE=ax_maj, SLAT=slat0, LON0_dir=lon0)
  assert np.max(LAT) < 0., f"For southern hemisphere latitudes should be < 0"
  LON = -LON   # nor sure why but this makes sign of the longitudes right
else:
  # Get gmapi 4 NSIDC grid points for interpolation
  pthdump  = os.path.join(pthdata,"gmapi_NSIDC")
  fgmapi  = f'NSIDC_NRTice_MOM6_gmapi_{idm}x{jdm}_{regn}.nc'
  dfgmapi = os.path.join(pthdump, fgmapi)
  print(f'Loading gmapi --> {dfgmapi}')
  with xarray.open_dataset(dfgmapi) as dgmapi:
    LON = dgmapi['longit'].data
    LAT = dgmapi['latit'].data
 

print(f"Processing {YR}/{MM}/{DD} {regn} ...")
AA = read_NSIDC(YR, MM, DD, regn, pthnsidc, 'cdr_seaice_conc')

# Interpolated fields:
fliceout = f'NSIDC_iconc_interp_mesh025_{jdm}x{idm}_{YR}{MM:02d}_{regn}.nc'
dfliceout = os.path.join(pthnsidc,fliceout)
print(f'Loading interpolated ice conc {dfliceout}')
with xarray.open_dataset(dfliceout) as dsint:
  AI = dsint['ice_conc'].isel(time=DD-1).data.squeeze()


clrmp = mclrmps.colormap_conc()
rmin = 0.
rmax = 1.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

LON1 = LON.copy()
LON2 = LON.copy()
LON1 = np.where(LON1 < -175, np.nan, LON1)
LON2 = np.where(LON2 > 172, np.nan, LON2)
LON3 = np.where(LON < 0, LON+360., LON)
LON3 = np.where(LON3 > 350., np.nan, LON3)
lon_cntr1 = [x for x in range(-180,0,45)]  # grey -180:0
lon_cntr2 = [x for x in range(45,178,45)]  # blue: 0 to 180 E
lat_cntr = [x for x in range(-80,-20,10)]

print("Plotting ...")

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.05, 0.3, 0.4, 0.4])

ax1.pcolormesh(AA, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.contour(LON1, lon_cntr1, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
ax1.contour(LON2, lon_cntr2, linestyles='solid', colors=[(0.,0.5,0.9)], linewidths=1)
ax1.contour(LON1,[0], linestyles='solid', colors=[(1,0.,0.)], linewidths=1)
ax1.contour(LON3,[180], linestyles='solid', colors=[(0.6,0.,0.9)], linewidths=1)
ax1.contour(LAT, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
ax1.invert_yaxis()
ax1.axis('scaled')
ax1.set_title(f'NSIDC iconc {YR}/{MM:02d}/{DD:02d}')


# Interpolated iconc
if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere', boundinglat=50, lon_0=-45, resolution='l')
  parallels = np.arange(40, 89, 10.)
  meridians = np.arange(-360, 359., 45.)

#lons, lats = m.makegrid(idim, jdim) # get lat/lons of ny by nx evenly spaced grid.
xh, yh = m(hlon,hlat) # GFS coords

ax2 = plt.axes([0.55, 0.3, 0.4, 0.4])
# draw parallels.
if regn == 'south':
  parallels = np.arange(-80,-10,10.)
else:
  parallels = np.arange(50,89,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-360,359.,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

img = ax2.pcolormesh(xh, yh, AI, cmap=clrmp, vmin=rmin, vmax=rmax)
ax2.set_title('NSIDC iconc interp to mesh025')

ax3 = fig1.add_axes([0.2, 0.2, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'check_interp_NSIDC_mesh025.py'
bottom_text(btx)


f_ij = False
if f_ij:
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  img1 = ax1.pcolormesh(AI, cmap=clrmp, vmin=rmin, vmax=rmax)







