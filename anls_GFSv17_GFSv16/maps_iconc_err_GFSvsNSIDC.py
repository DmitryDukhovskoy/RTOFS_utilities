"""
  Error and dfference maps of ice conc
  derived from datmUFS and NSIDC NRT 

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

init_date = 20250604
init_hr = 6
regn = 'south'

runnm = 'retrov17_01' 
strnm = '4'
hrS = 0      # 1st forecast, hr, =0 - initial state
hrE = 384    # last forecast, hr
dltHR = 6    # output time freq, hrs

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", type=str, required=True)
parser.add_argument("--runnm", help=f"Name of the run, default={runnm}", type=str)
parser.add_argument("--strnm", help=f"Stream, default={strnm}", choices=['1a','4'], type=str)
parser.add_argument("--init", help=f"init date, default {init_date}", type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fhr", help=f"f/cast hr to plot: {hrS}:{hrE} default={hrE}", type=int)
args = parser.parse_args()
  
regn      = args.regn if args.regn else None
run_name  = args.runnm if args.runnm else runnm
str_name  = args.strnm if args.strnm else strnm
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fcst_hr   = args.fhr if args.fhr is not None else hrE

assert hrS <= fcst_hr <= hrE, f"Requested f/cast hour {fcst_hr} is outside time range: {hrS}/{hrE}"
 
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

run_full = f"{run_name}_stream{str_name}"
if init_date is not None and init_hr is not None:
  RUNS = [f"{init_date}{init_hr:02d}"]
else:
  RUNS = mgfscice.gfs_retro_runs(run_full, node_nm, model='ice')

fyaml = 'gfs17_paths.yaml'
with open(fyaml) as ff:
  pths_gfs = safe_load(ff)

# Get MOM6 grid
pthgrid = pths_gfs[node_nm]["MOM6"]["pthgrid"]
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
  elif regn == 'north':
    flnsidc = f"sic_psn25_{YR}{MM:02d}{DD:02d}_am2_v06r00.nc"  

  with xarray.open_dataset(os.path.join(pthnsidc,flnsidc)) as ds_nsidc:
    A = ds_nsidc[varnm].data.squeeze()

  return A

# Get date:
plot_init = fcst_hr == 0  # initial conditions

# Output to plot:
# Assumed all runs have same forecast duration
hrs_outp = np.array([x for x in range(hrS, hrE+1, dltHR)])
dhr = np.abs(hrs_outp - fcst_hr)
iplot = np.argmin(dhr)
fcst_hr = hrs_outp[iplot]  # correct f/cast hour if needed

dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb

# Date to plot:
dnmb0 = dnmbI + fcst_hr/24
yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
YR,MM,DD = mtime.datevec(dnmb0)[:3] 
nsec0 = hr0*3600


# Get lon/lat for NSIDC data
pthdata = pths_gfs[node_nm]["CICE6"]["pthdata"]
pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")

RMsk = np.where(HH>=0, 0, 1)
if regn == 'south':
  RMsk = np.where(hlat > -50., 0, RMsk)
elif regn == 'north':
  RMsk = np.where(hlat < 40, 0., RMsk)


comroot = pths_gfs[node_nm]["CICE6"]["comroot"].format(
  run_name=run_name,
  stream=str_name,
  init=init_date,
  ihr=init_hr
  )

pthout_cice = os.path.join(pths_gfs[node_nm]["CICE6"]["pthcice"].format(
  comroot=comroot
))
print(f"cice dir: {pthout_cice}")


if plot_init:
  flinp = f"gfs.t{init_hr:02d}z.ic.nc"
else:
  flinp = f"gfs.t{init_hr:02d}z.{dltHR}hr_avg.f{fcst_hr:03d}.nc"

dflice = os.path.join(pthout_cice,flinp)

print(f"Processing {YR}/{MM}/{DD}:{fcst_hr:02d} run {init_date}:{init_hr:02d}...")
print(f"Processing {dflice}")
with xarray.open_dataset(dflice) as dcice:
  AA = dcice['aice_h'].data.squeeze()

# Interpolated fields:
fliceout = f'NSIDC_iconc_interp_mesh025_{jdm}x{idm}_{YR}{MM:02d}_{regn}.nc'
dfliceout = os.path.join(pthnsidc,fliceout)
print(f'Loading interpolated ice conc {dfliceout}')
with xarray.open_dataset(dfliceout) as dsint:
  AI = dsint['ice_conc'].isel(time=DD-1).squeeze()

AA = np.where(RMsk == 0, np.nan, AA)
AI = np.where(RMsk == 0, np.nan, AI)
sqerr = (AA-AI)**2
abserr = np.abs(AA-AI)
diff  = (AA-AI)

plt.ion()


clrmp = mclrmps.colormap_conc()
rmin = 0.
rmax = 1.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

clrmp_dlt = mclrmps.colormap_uv()
dmin = -1
dmax = 1
clrmp_dlt.set_bad(color=[0.2, 0.2, 0.2])

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
  xl1 = -8.5e6
  xl2 = -1.e6
  yl1 = xl1
  yl2 = xl2 
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=50,lon_0=0,resolution='l')
  parallels = np.arange(40,89,10.)
  meridians = np.arange(-360,359.,45.)
  xl1 = None
  xl2 = None
  yl1 = xl1
  yl2 = xl2 

runname_date =f"{run_name}_stream{str_name} {init_date}{init_hr:02d}"

xh, yh = m(hlon,hlat) # GFS coords
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.05, 0.55, 0.4, 0.4])
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
img1 = ax1.pcolormesh(xh,yh,AA, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(f'{runname_date}  iconc {YR}/{MM:02d}/{DD:02d}')
if regn == 'south':
  ax1.set_xlim([xl1, xl2]) 
  ax1.set_ylim([yl1, yl2]) 
  ax1.invert_yaxis()
  ax1.invert_xaxis()

# Interpolated iconc
ax2 = plt.axes([0.55, 0.55, 0.4, 0.4])
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
img2 = ax2.pcolormesh(xh, yh, AI, cmap=clrmp, vmin=rmin, vmax=rmax)
ax2.set_title('NSIDC iconc interp to mesh025')
if regn == 'south':
  ax2.set_xlim([xl1, xl2])      
  ax2.set_ylim([yl1, yl2])      
  ax2.invert_yaxis()
  ax2.invert_xaxis()

ax21 = plt.axes([0.05, 0.1, 0.4, 0.4])
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ax21.pcolormesh(xh,yh,np.abs(AA-AI), cmap=clrmp, vmin=rmin, vmax=rmax)
ax21.set_title(f'|err| iconc GFS vs  NSIDC {YR}/{MM:02d}/{DD:02d}')
if regn == 'south':
  ax21.set_xlim([xl1, xl2])      
  ax21.set_ylim([yl1, yl2])      
  ax21.invert_yaxis()
  ax21.invert_xaxis()

ax22 = plt.axes([0.55, 0.1, 0.4, 0.4])
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
img2 = ax22.pcolormesh(xh,yh,(AA-AI), cmap=clrmp_dlt, vmin=dmin, vmax=dmax)
ax22.set_title(f'diff iconc GFS vs NSIDC {YR}/{MM:02d}/{DD:02d}')
if regn == 'south':
  ax22.set_xlim([xl1, xl2])      
  ax22.set_ylim([yl1, yl2])      
  ax22.invert_yaxis()
  ax22.invert_xaxis()


# Colorbars
ax3 = fig1.add_axes([0.05, 0.05, 0.4, 0.02])
clb = plt.colorbar(img1, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ax4 = fig1.add_axes([0.55, 0.05, 0.4, 0.02])
clb = plt.colorbar(img2, cax=ax4, orientation='horizontal', extend='both')
ax4.xaxis.set_ticks(list(np.linspace(dmin,dmax,11)))
ax4.set_xticklabels(ax4.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'maps_iconc_err_GFSvsNSIDC.py'
bottom_text(btx, pos=[0.1,0.01])



