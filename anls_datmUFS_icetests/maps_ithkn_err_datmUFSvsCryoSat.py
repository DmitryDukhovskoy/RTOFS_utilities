"""
  Error and dfference maps of ice thickness
  derived from datmUFS and monthly clim CryoSat fields

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

expt = 'ufs_datm_mx025_v02'
#init_date = 20250103
init_hr = 0
regn = 'south'

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--init", help=f"init date: YYYYMMDD", required=True, type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--fday", help=f"forecast day to plot: 1,...,14, =0 - init. cond.", type=int, required=True)
parser.add_argument("--enmb", help="experiment number: 1, 2, ...", type=int, required=True)
args = parser.parse_args()
  
enmb      = args.enmb 
init_date = args.init
init_hr   = args.ihr if args.ihr else init_hr
fday      = args.fday 
regn      = args.regn if args.regn else regn
  
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

jdim, idim = HH.shape

def read_CryoSat(YR,MM,DD,regn,pthnsidc,varnm):
  if regn == 'south':
    flnsidc = f"sic_pss25_{YR}{MM:02d}{DD:02d}_am2_v06r00.nc"  
  elif regn == 'north':
    flnsidc = f"sic_psn25_{YR}{MM:02d}{DD:02d}_am2_v06r00.nc"  

  with xarray.open_dataset(os.path.join(pthnsidc,flnsidc)) as ds_nsidc:
    A = ds_nsidc[varnm].data.squeeze()

  return A

def find_varnm(dflithkn, var_opt):
  with xarray.open_dataset(dflithkn) as ds_ithkn:
    for varnm in var_opt:
      if varnm in ds_ithkn.data_vars:
        #print(f"Using variable {varnm}")
        return varnm 
        
  raise KeyError("No ice thickness variable name found, check file")
  
  return

def ice_clim_files(enmb, regn, pths_ufs):
  pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
  if enmb >= 30:
    fyaml_rest = 'cice6rest_files.yaml'
    with open(fyaml_rest) as fy:
      pths_clim = safe_load(fy)


    pthice = pths_clim["target_paths"][f"ithkn_{regn}"]["path"]
    fhice = pths_clim["target_paths"][f"ithkn_{regn}"]["file"]

    var_opt = ['ice_thkn', 'ithkn', 'hi', 'ice_thickness']
    varnm = find_varnm(os.path.join(pthice, fhice), var_opt)


  else:
    if regn == 'south':
      pthice = os.path.join(pthdata,'CryoSat2_antarctic_ice_snow_thkn','clim')
      fhice  = f'CryoSat_hice_mnthclim_2011_2020_mesh025_1440x1080_{regn}.nc'
      varnm = 'ice_thkn'
    elif regn == 'north':
      pthice = os.path.join(pthdata,'CryoSat_arctic_ice_snow_thkn','clim')
      fhice = f'ithkn_CryoSat_arcticAWI_mnthclim_2015-2024_1080x1440.nc'
      varnm = 'ithkn'

  return pthice, fhice, varnm

# Get date:
plot_init = fday == 0  # initial conditions

dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
if plot_init:
  dnmb0 = dnmbI
else:
  dnmb0 = dnmbI + fday-1                              # day to plot

yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
YR,MM,DD = mtime.datevec(dnmb0)[:3] 
nsec0 = hr0*3600


# Get ithkn monthly clim interpolated to mesh025
pthice, fhice, varnm = ice_clim_files(enmb, regn, pths_ufs)

dfhice = os.path.join(pthice, fhice)

print(f'Reading ice thickn climatology {dfhice}')
with xarray.open_dataset(dfhice) as ds_hice:
  AI = ds_hice[varnm].isel(time=MM-1).squeeze()

AI = np.where(np.isnan(AI), 0., AI)
AI = np.where(HH>=0, np.nan, AI)

RMsk = np.where(HH>=0, 0, 1)
#if regn == 'south':
# RMsk = np.where(hlat > -60., 0, RMsk)

if plot_init:
  flinp = f"iceh_ic.{yr0}-{mm0:02d}-{dd0:02d}-{nsec0:05d}.nc"
else:
  flinp = f"iceh.{yr0}-{mm0:02d}-{dd0:02d}.nc"

pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/ufs_datm_mx025/expt{enmb:02d}/cice6"
dflice = os.path.join(pthoutp,flinp)

print(f"Processing {YR}/{MM}/{DD} expt{enmb:02d} {dflice}...")
with xarray.open_dataset(dflice) as dcice:
  AICE = dcice['aice_d'].data.squeeze()
  AA = dcice['hi_d'].data.squeeze()    # grid cell mean ice thickness


AA = np.where(RMsk == 0, np.nan, AA)
AI = np.where(RMsk == 0, np.nan, AI)

# Ignore open ocean grid points, only where sea ice is present in the f/cast:
IMask = AICE>0.001    # boolean mask
sqerr = (AA - AI)**2 * IMask
abserr = np.abs(AA - AI) * IMask
diff  = (AA - AI) * IMask

plt.ion()


#clrmp = mclrmps.colormap_conc()
clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.2, 0.2, 0.2])

clrmp_dlt = mclrmps.colormap_uv()
dmin = -2
dmax = 2
clrmp_dlt.set_bad(color=[0.2, 0.2, 0.2])

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=50,lon_0=-10,resolution='l')
  parallels = np.arange(40,89,10.)
  meridians = np.arange(-360,359.,45.)


xh, yh = m(hlon,hlat) # GFS coords
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.05, 0.55, 0.4, 0.4])
m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
img1 = ax1.pcolormesh(xh,yh,AA, cmap=clrmp, vmin=rmin, vmax=rmax)
#img1 = ax1.pcolormesh(xh,yh,sqerr, cmap=clrmp, vmin=rmin, vmax=rmax)
#img1 = ax1.pcolormesh(xh,yh,np.abs(AA-AI), cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.set_title(f'UFS expt{enmb:02d} ithkn {YR}/{MM:02d}/{DD:02d}')

# Interpolated ithkn
ax2 = plt.axes([0.55, 0.55, 0.4, 0.4])
m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
img2 = ax2.pcolormesh(xh, yh, AI, cmap=clrmp, vmin=rmin, vmax=rmax)
ax2.set_title('CryoSat ithkn interp to mesh025')

ax21 = plt.axes([0.05, 0.1, 0.4, 0.4])
m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
ax21.pcolormesh(xh, yh, abserr, cmap=clrmp, vmin=rmin, vmax=rmax)
ax21.set_title(f'|err| ithkn UFS vs  CryoSat {YR}/{MM:02d}/{DD:02d}')

ax22 = plt.axes([0.55, 0.1, 0.4, 0.4])
m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
img2 = ax22.pcolormesh(xh, yh, diff, cmap=clrmp_dlt, vmin=dmin, vmax=dmax)
ax22.set_title(f'diff ithkn UFS vs CryoSat {YR}/{MM:02d}/{DD:02d}')


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
clb.ax.set_xticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'maps_ithkn_err_datmUFSvsCryoSat.py'
bottom_text(btx, pos=[0.1,0.01])



